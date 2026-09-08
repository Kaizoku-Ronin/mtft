#!/usr/bin/env python3
"""Replay the frozen T2 study, including independent real-frame comparison."""
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

ROOT=Path(__file__).resolve().parent
ARCHIVE_SHA='35cfc00c32877bd6c12b16464ddc8190bcca9e2f054474f4ab0a6c605fe9da38'
PROTOCOL_SHA='82a766f60ef267b7246bb09469274f12e67b69c25ee0c3b05f9b14353c61f2cc'
INPUT_SHA={
    'coupling_matrices.npz':'1b011327896525ffca8d53bbd0279871e102016b2035507575f13614a2543f13',
    'kernel_check_projectors.npz':'eaeef3f1b6c765b1bb3d676e219bbc2e2b41ae6f295ef4a30d12177a6329c36b',
}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify_source(archive,source):
    count=0
    with tarfile.open(archive,'r:gz') as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            relative=Path(member.name).relative_to('mtft-0.26.2')
            actual=source/relative
            if not actual.is_file() or actual.read_bytes()!=tf.extractfile(member).read():
                raise RuntimeError(f'Source does not match supplied archive: {relative}')
            count+=1
    return count


def run(source,archive):
    os.environ['OPENBLAS_NUM_THREADS']='1'
    import numpy as np
    import scipy
    import sympy
    import mpmath
    import matplotlib
    previous=ROOT/'inputs'/'previous'
    if sha(ROOT/'PROTOCOL.md')!=PROTOCOL_SHA:
        raise RuntimeError('Frozen protocol hash mismatch')
    for name,expected in INPUT_SHA.items():
        if sha(previous/name)!=expected:
            raise RuntimeError(f'Frozen input mismatch: {name}')
    source_count=verify_source(archive,source)
    env=dict(os.environ)
    env['PYTHONPATH']=str(source/'src')+os.pathsep+env.get('PYTHONPATH','')
    commands=[
        ['t2_effective_operator_study.py','--source',str(source),'--previous',str(previous)],
        ['independent_T2_check.py','--source',str(source),'--previous',str(previous)],
        ['general_feedback_identities.py'],
        ['render_T2_effective_operator.py'],
    ]
    with (ROOT/'run_log.txt').open('w') as log:
        for command in commands:
            result=subprocess.run([sys.executable,*command],cwd=ROOT,env=env,
                                  text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            log.write(json.dumps(command)+'\n'+result.stdout+'\n');log.flush()
            if result.returncode:
                raise RuntimeError(f'{command[0]} failed; see run_log.txt')
    primary=json.loads((ROOT/'t2_effective_operator_results.json').read_text())
    independent=json.loads((ROOT/'independent_T2_results.json').read_text())
    symbolic=json.loads((ROOT/'general_feedback_identity_results.json').read_text())
    pm=np.load(ROOT/'t2_effective_operator_matrices.npz')
    im=np.load(ROOT/'independent_T2_matrices.npz')
    u=pm['U']
    ambient=np.asarray([u@g@u.conj().T for g in pm['regular_projected_responses']])
    response_errors=np.linalg.norm(ambient-im['ambient_projected_responses'],axis=(1,2))/np.linalg.norm(ambient,axis=(1,2))
    agreement={
        'ambient_projected_response_maximum_relative_error':float(np.max(response_errors)),
        'spectral_grid_maximum_difference':float(np.max(abs(pm['regular_z']-im['grid']))),
        'coupling_singular_values_maximum_difference':float(np.max(abs(
            np.asarray(primary['structure']['coupling_singular_values'])-
            independent['coupling_complex_singular_values_from_real_pair_means']))),
        'complement_eigenvalues_maximum_difference':float(np.max(abs(
            np.asarray(primary['structure']['complement_eigenvalues'])-
            independent['compression_complex5_eigenvalues_from_real_pair_means']))),
        'full_eigenvalues_maximum_difference':float(np.max(abs(
            pm['eigenvalues']-independent['full_complex13_eigenvalues_from_real_pair_means']))),
        'active_projector_Frobenius_difference':float(np.linalg.norm(pm['P']-im['P_complex'])),
    }
    exp_keys=['spectral_full_relative_error','full_unitarity','projected_spectral_relative_error',
              'norm_balance_residual','memory_spectral_relative_error']
    exponential_errors={key:max(row[key] for row in primary['exponential_controls']) for key in exp_keys}
    st=primary['structure'];su=primary['summary'];tol=1e-9
    checks={
        'pinned_mtft_version':primary['mtft_version']==independent['version']=='0.26.2',
        'raw_primary_geometry_gates':max(primary['raw_checks'].values())<=tol,
        'independent_real_frame_checks':independent['all_checks_passed'],
        'forty_symbolic_controls':symbolic['checks_passed']==symbolic['checks_total']==40
            and all(symbolic['checks'].values()),
        'regular_grid_1807_without_skips':su['regular_count']==1807 and su['skipped_regular_count']==0,
        'registered_regular_Schur_gate':su['max_regular_effective_relative_error']<=tol,
        'registered_regular_spectral_gate':su['max_regular_spectral_relative_error']<=tol,
        'independent_response_agrees':agreement['ambient_projected_response_maximum_relative_error']<=tol,
        'independent_grid_agrees':agreement['spectral_grid_maximum_difference']==0,
        'independent_coupling_singular_values_agree':agreement['coupling_singular_values_maximum_difference']<=tol,
        'independent_complement_spectrum_agrees':agreement['complement_eigenvalues_maximum_difference']<=tol,
        'independent_full_spectrum_agrees':agreement['full_eigenvalues_maximum_difference']<=tol,
        'independent_active_projectors_agree':agreement['active_projector_Frobenius_difference']<=tol,
        'coupling_rank_routes_agree':2*st['coupling_rank']['resolved_rank']==independent['coupling_rank_real_at_frozen_bands']
            and not st['coupling_rank']['ambiguous'] and not independent['coupling_rank_ambiguous'],
        'six_exponential_controls_agree':len(primary['exponential_controls'])==6 and max(exponential_errors.values())<=tol,
        'memory_residues_sum_correctly':st['residue_sum_relative_residual']<=tol,
        'memory_residues_positive_to_roundoff':all(min(x['residue_eigenvalues'])>=-tol for x in primary['compression_poles']),
        'every_selected_pole_has_eleven_stress_points':len(primary['pole_stress'])==11*len(primary['real_compression_poles']),
        'singular_complement_inverses_skipped':all(x['Schur_evaluation'].startswith('SKIPPED') for x in primary['real_compression_poles']),
        'projected_pole_null_identity':all(x['null_relation_relative_residual']<=tol for x in primary['real_compression_poles']),
        'instantaneous_dark_directions_annihilated':st['instantaneous_dark_coupling_relative']<=tol,
        'registered_dynamics_and_short_time_counts':len(primary['dimensionless_dynamics'])==257 and len(primary['short_time_controls'])==6,
        'sampled_transfer_respects_unitarity_bounds':all(-tol<=row[key]<=1+tol
            for row in primary['dimensionless_dynamics']
            for key in ['average_transfer','worst_input_transfer','dark_average_transfer','dark_worst_input_transfer']
            if row[key] is not None),
        'finite_response_diagnostics':all(np.isfinite(row[key])
            for row in primary['regular_response']+primary['pole_stress']
            for key in ['effective_relative_error','spectral_relative_error','naive_relative_error','Schur_condition']),
    }
    validation={
        'status':'PASS' if all(checks.values()) else 'FAIL',
        'scope':'Replay and identity consistency. Structural hypotheses and raw Schur stress accuracy are separately reported outcomes, not assumed to pass.',
        'checks':checks,'checks_passed':sum(checks.values()),'checks_total':len(checks),
        'independent_agreement':agreement,'exponential_control_maximum_errors':exponential_errors,
        'response_summary':su,
        'structural_outcomes':{
            'coupling_complex_rank':st['coupling_rank'],
            'pair_decomposition_relative_residual':st['independent_pair_decomposition_relative_residual'],
            'instantaneous_dark_dimension':st['dark_complex_dimension'],
            'delayed_dark_rank':st['delayed_dark_rank'],
            'complement_reachable_rank':st['complement_reachable_rank'],
        },
        'short_time_observations':{
            'ordinary_derived_coefficient':st['ordinary_short_time_coefficient'],
            'dark_derived_coefficient':st['dark_short_time_coefficient'],
            'smallest_time_record':primary['short_time_controls'][-1],
            'interpretation':'Taylor coefficients are derived. Finite-time convergence is reported without fitting an exponent or treating its nonzero remainder as an identity failure.'},
        'stress_note':'Raw Schur evaluation loses the regular-grid accuracy near compression poles. All stress points are retained; no stress accuracy gate was substituted.',
        'source_archive_sha256':sha(archive),'source_files_verified_unchanged':source_count,
        'protocol_sha256':sha(ROOT/'PROTOCOL.md'),
        'frozen_input_sha256':{name:sha(previous/name) for name in INPUT_SHA},
        'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
                   'sympy':sympy.__version__,'mpmath':mpmath.__version__,'matplotlib':matplotlib.__version__},
    }
    (ROOT/'VALIDATION.json').write_text(json.dumps(validation,indent=2)+'\n')
    np.savez_compressed(ROOT/'independent_comparison.npz',grid=pm['regular_z'],
                        ambient_response_relative_errors=response_errors)
    print(json.dumps(validation,indent=2))
    if not all(checks.values()):
        raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,help='Optional matching extracted MTFT source')
    args=parser.parse_args()
    archive=ROOT/'inputs'/'mtft-0.26.2.tar.gz'
    if sha(archive)!=ARCHIVE_SHA:
        raise RuntimeError('Source archive hash mismatch')
    if args.source:
        run(args.source.resolve(),archive)
    else:
        with tempfile.TemporaryDirectory(prefix='mtft_t2_replay_') as temporary:
            with tarfile.open(archive,'r:gz') as tf:
                tf.extractall(temporary,filter='data')
            run(Path(temporary)/'mtft-0.26.2',archive)


if __name__=='__main__':
    main()
