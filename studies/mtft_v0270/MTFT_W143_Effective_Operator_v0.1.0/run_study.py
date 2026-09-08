#!/usr/bin/env python3
"""Replay the frozen W143 experiment and check its documented outcomes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parent
ARCHIVE_SHA256 = '35cfc00c32877bd6c12b16464ddc8190bcca9e2f054474f4ab0a6c605fe9da38'
PROTOCOL_SHA256 = 'd6e08abc6b41e1aab4099653d1a8eb71df1479b103c6b751a007d3812475f18d'
INPUT_SHA256 = {
    'coupling_matrices.npz': '1b011327896525ffca8d53bbd0279871e102016b2035507575f13614a2543f13',
    'kernel_check_projectors.npz': 'eaeef3f1b6c765b1bb3d676e219bbc2e2b41ae6f295ef4a30d12177a6329c36b',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_check(archive, source):
    """Verify every regular file from the supplied source distribution."""
    checked = 0
    with tarfile.open(archive, 'r:gz') as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            relative = Path(member.name).relative_to('mtft-0.26.2')
            actual = source / relative
            expected = tf.extractfile(member).read()
            if not actual.is_file() or actual.read_bytes() != expected:
                raise RuntimeError(f'Source differs from supplied archive: {relative}')
            checked += 1
    return checked


def execute(source, archive):
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    import numpy as np
    import scipy
    import sympy
    import mpmath
    import matplotlib

    previous = ROOT / 'inputs' / 'previous'
    if digest(ROOT / 'PROTOCOL.md') != PROTOCOL_SHA256:
        raise RuntimeError('Frozen protocol hash mismatch')
    for name, expected in INPUT_SHA256.items():
        if digest(previous / name) != expected:
            raise RuntimeError(f'Frozen input hash mismatch: {name}')
    checked = source_check(archive, source)
    env = dict(os.environ)
    env['PYTHONPATH'] = str(source / 'src') + os.pathsep + env.get('PYTHONPATH', '')
    commands = [
        ['effective_operator_study.py', '--source', str(source), '--previous', str(previous)],
        ['canonical_check.py', '--source', str(source), '--previous', str(previous)],
        ['feedback_identities.py'],
        ['render_effective_operator.py'],
    ]
    with (ROOT / 'run_log.txt').open('w') as log:
        for command in commands:
            result = subprocess.run([sys.executable, *command], cwd=ROOT, env=env,
                                    text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log.write(json.dumps(command) + '\n' + result.stdout + '\n')
            log.flush()
            if result.returncode:
                raise RuntimeError(f'{command[0]} failed; see run_log.txt')

    primary = json.loads((ROOT / 'effective_operator_results.json').read_text())
    independent = json.loads((ROOT / 'canonical_check_results.json').read_text())
    symbolic = json.loads((ROOT / 'feedback_identity_results.json').read_text())
    summary = primary['summary']
    tol = 1e-9
    modes = sorted(primary['modes'], key=lambda x: -x['b'])
    other = sorted(independent['canonical_modes'], key=lambda x: -x['b'])
    agreement = {
        name: max(abs(x[name] - y[name]) for x, y in zip(modes, other))
        for name in ['a', 'b', 'principal_angle_degrees']
    }
    dynamics = primary['dimensionless_unitary_control']['matrix_exponential_checks']
    dynamic_errors = {
        name: max(x[name] for x in dynamics)
        for name in ['full_unitarity', 'projected_trig_identity', 'loss_identity']
    }
    dynamic_errors['mode_fraction_error'] = max(
        abs(v) for x in dynamics for v in x['mode_fraction_errors'])
    matrices = np.load(ROOT / 'effective_operator_matrices.npz')
    external = np.load(ROOT / 'canonical_check_matrices.npz')
    projector_agreement = float(np.linalg.norm(matrices['P'] - external['P']))
    checks = {
        'supplied_mtft_version': primary['mtft_version'] == independent['mtft_version'] == '0.26.2',
        'raw_input_geometry_gate': max(primary['input_checks'].values()) <= tol,
        'block_identity_and_reconstruction_gate': max(primary['block_checks'].values()) <= tol,
        'two_unambiguous_complex_couplings': primary['coupling_rank_complex'] == 2 and not primary['rank_ambiguous'],
        'two_independent_canonical_modes': len(modes) == len(other) == 2,
        'independent_coefficients_agree': max(agreement.values()) <= tol,
        'independent_active_projectors_agree': projector_agreement <= tol,
        'independent_route_checks': independent['all_checks_passed'],
        'uncoupled_active_six_negative_directions': len(primary['uncoupled_active_eigenvalues']) == 6
            and max(abs(x + 1) for x in primary['uncoupled_active_eigenvalues']) <= tol,
        'uncoupled_complement_three_negative_directions': len(primary['uncoupled_fixed_eigenvalues']) == 3
            and max(abs(x + 1) for x in primary['uncoupled_fixed_eigenvalues']) <= tol,
        'registered_regular_grid_count': len(primary['regular_response']) == 907,
        'registered_regular_Schur_gate': summary['max_regular_effective_relative_error'] <= tol,
        'registered_regular_involution_gate': summary['max_regular_involution_relative_error'] <= tol,
        'registered_stress_points_recorded': len(primary['pole_stress']) == 22,
        'singular_complement_inverses_skipped': len(primary['real_compression_poles']) == 2
            and all(x['Schur_evaluation'].startswith('SKIPPED') for x in primary['real_compression_poles']),
        'six_exponential_controls_agree': len(dynamics) == 6 and max(dynamic_errors.values()) <= tol,
        'all_symbolic_controls_pass': symbolic['checks_passed'] == symbolic['checks_total'] == 27
            and all(symbolic['checks'].values()),
        'finite_response_diagnostics': all(np.isfinite(row[key])
            for row in primary['regular_response'] + primary['pole_stress']
            for key in ['effective_relative_error', 'involution_relative_error',
                        'naive_relative_error', 'Schur_denominator_condition']),
    }
    validation = {
        'status': 'PASS' if all(checks.values()) else 'FAIL',
        'scope': 'Replay validation. Registered response gate is the regular grid only; raw Schur stress loss is retained.',
        'checks': checks, 'checks_passed': sum(checks.values()), 'checks_total': len(checks),
        'source_archive_sha256': digest(archive), 'archive_files_verified_unchanged': checked,
        'protocol_sha256': digest(ROOT / 'PROTOCOL.md'),
        'frozen_input_sha256': {name: digest(previous / name)
                              for name in ['coupling_matrices.npz', 'kernel_check_projectors.npz']},
        'maximum_independent_coefficient_differences': agreement,
        'independent_projector_difference': projector_agreement,
        'dynamic_control_maximum_errors': dynamic_errors,
        'response_summary': summary,
        'stress_warning': 'Raw Schur evaluation exceeds 1e-9 in the registered pole-conditioning stress test. It is not certified by the regular-grid pass.',
        'runtime': {'python': platform.python_version(), 'numpy': np.__version__,
                    'scipy': scipy.__version__, 'sympy': sympy.__version__,
                    'mpmath': mpmath.__version__, 'matplotlib': matplotlib.__version__},
    }
    (ROOT / 'VALIDATION.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(validation, indent=2))
    if not all(checks.values()):
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, help='Optional matching extracted MTFT source tree')
    args = parser.parse_args()
    archive = ROOT / 'inputs' / 'mtft-0.26.2.tar.gz'
    if digest(archive) != ARCHIVE_SHA256:
        raise RuntimeError('Supplied source archive hash mismatch')
    if args.source:
        execute(args.source.resolve(), archive)
    else:
        with tempfile.TemporaryDirectory(prefix='mtft_w143_replay_') as temp:
            with tarfile.open(archive, 'r:gz') as tf:
                tf.extractall(temp, filter='data')
            execute(Path(temp) / 'mtft-0.26.2', archive)


if __name__ == '__main__':
    main()
