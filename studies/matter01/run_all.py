"""Reproduce the Matter-01 milestone, print checks, and retain every result.

Run from any directory: python path/to/studies/matter01/run_all.py.
This is the focused milestone gate, not all archived MTFT long experiments.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
DEFAULT_REPO = ROOT / 'repo_candidate' if (ROOT / 'repo_candidate').exists() else ROOT.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=DEFAULT_REPO)
    args = parser.parse_args()
    repo = args.repo.resolve()
    if not (repo / 'src/mtft/matter/qed.py').exists():
        raise SystemExit('Matter-01 source not found. Pass --repo path/to/patched/mtft.')
    results = ROOT / 'results'
    results.mkdir(exist_ok=True)
    env = os.environ.copy()
    env['PYTHONPATH'] = str(repo / 'src') + os.pathsep + env.get('PYTHONPATH', '')
    log = []

    def run(name, *command, cwd=ROOT):
        print(f'\n{name}', flush=True)
        start = time.perf_counter()
        finished = subprocess.run(list(map(str, command)), cwd=cwd, env=env,
                                  text=True, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT)
        output = finished.stdout
        (results / f'{name}.log').write_text(output, encoding='utf-8')
        log.append({'stage': name, 'returncode': finished.returncode,
                    'seconds': round(time.perf_counter()-start, 3),
                    'log': f'results/{name}.log'})
        # Full logs remain on disk; concise terminal output includes all verdicts.
        lines = output.strip().splitlines()
        print('\n'.join(lines[-10:]) if len(lines)>10 else output, flush=True)
        if finished.returncode:
            (results / 'run_summary.json').write_text(json.dumps({'stages': log, 'passed': False}, indent=2)+'\n')
            raise SystemExit(f'{name} failed; read results/{name}.log')

    python = sys.executable
    run('vacuum_audit', python, ROOT/'studies/vacuum_audit.py', '--repo', repo)
    run('u1_reduction', python, ROOT/'studies/verify_u1_reduction.py')
    run('qed_reference', python, '-m', 'mtft.matter', 'report', '--output', results/'qed_reference_report.json')
    run('independent_qed', python, ROOT/'studies/verify_qed_oracle.py')
    run('atlas_ledger', python, ROOT/'studies/build_evidence_ledger.py')
    run('focused_tests', python, '-m', 'pytest',
        repo/'tests/test_matter01_qed.py', repo/'tests/test_matter01_normalization.py',
        repo/'tests/test_interactions.py', repo/'tests/test_pv02.py',
        repo/'tests/test_repository_docs.py', '-q', cwd=repo)
    run('compose_lab', python, ROOT/'viz/compose_evidence.py', '--engine', results/'qed_reference_report.json')
    run('build_lab', python, ROOT/'viz/build_lab.py', ROOT/'viz/lab_evidence.json')
    node = shutil.which('node')
    if node:
        run('javascript_checks', node, ROOT/'viz/verify_lab.cjs', ROOT/'viz/Matter01_U1_Lab.html',
            ROOT/'viz/lab_runtime_checks.json')
    else:
        print('Node.js not available: JavaScript runtime checks SKIPPED; Python checks and prebuilt lab remain available.')
        log.append({'stage': 'javascript_checks', 'status': 'SKIPPED: Node.js unavailable'})
    run('scientific_figure', python, ROOT/'studies/plot_matter01.py', '--output', ROOT/'viz/Matter01_U1_Overview.png')

    vacuum = json.loads((results/'vacuum_audit.json').read_text())
    qed = json.loads((results/'qed_reference_report.json').read_text())
    oracle = json.loads((results/'qed_oracle_checks.json').read_text())
    reduction = json.loads((results/'u1_reduction_checks.json').read_text())
    passed = (vacuum['summary']['n_failed']==0 and qed['verification']['passed']
              and oracle['all_passed'])
    dependencies = {name: importlib.metadata.version(name)
                    for name in ('numpy','sympy','scipy','pytest','matplotlib')}
    summary = {
        'milestone': 'Matter-01', 'passed': passed,
        'source_version': '0.34.0 + unreleased Matter-01 candidate',
        'python': platform.python_version(), 'platform': platform.platform(),
        'dependencies': dependencies, 'stages': log,
        'exact_vacuum_checks': vacuum['summary']['n_exact_checks'],
        'exact_reduction_identities': len(reduction['exact_gates']),
        'numerical_breathing_samples': len(reduction['breathing_implementation_samples']),
        'qed_structural_gates': qed['verification']['count'],
        'independent_qed_gates': len(oracle['gates']),
        'browser_layout': 'UNVERIFIED; Node mock-DOM checks are not a browser layout test',
        'physical_electron_photon_derivation': 'OPEN',
        'full_MTFT_release_suite': 'NOT RUN by this focused milestone runner',
    }
    (results/'run_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('\nMatter-01:', 'PASS' if passed else 'FAIL')
    print(f"Exact checks: {summary['exact_vacuum_checks']} vacuum + {summary['exact_reduction_identities']} reduction identities.")
    print(f"QED numerical gates: {summary['qed_structural_gates']} structural + {summary['independent_qed_gates']} independent.")
    print('Results: '+str(results))
    raise SystemExit(0 if passed else 1)


if __name__ == '__main__':
    main()
