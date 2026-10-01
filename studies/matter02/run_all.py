#!/usr/bin/env python3
"""Run the focused Matter-02 milestone and retain printed verdicts and logs.

This runs exact arithmetic certificates, regression tests, lab logic and a
scientific figure. It does not run the entire historical MTFT release suite.
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

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=HERE.parents[1])
    parser.add_argument('--skip-tests', action='store_true',
                        help='Run certificates and visuals only; label tests SKIPPED.')
    args = parser.parse_args()
    repo = args.repo.resolve()
    if not (repo/'src/mtft/matter/collinear.py').exists():
        raise SystemExit('Matter-02 source missing; pass --repo path/to/patched/repo.')
    results = HERE/'results'
    results.mkdir(exist_ok=True)
    env = os.environ.copy()
    env['PYTHONPATH'] = str(repo/'src') + os.pathsep + env.get('PYTHONPATH','')
    stages = []

    def run(name, *command):
        print('\n'+name, flush=True)
        start = time.perf_counter()
        process = subprocess.run(list(map(str,command)), cwd=repo, env=env,
                                 text=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
        (results/(name+'.log')).write_text(process.stdout, encoding='utf-8')
        stages.append({'name':name,'returncode':process.returncode,
                       'seconds':round(time.perf_counter()-start,3),
                       'log':'results/'+name+'.log'})
        lines = process.stdout.strip().splitlines()
        print('\n'.join(lines[-15:]), flush=True)
        if process.returncode:
            (results/'run_summary.json').write_text(json.dumps(
                {'passed':False,'stages':stages},indent=2)+'\n')
            raise SystemExit(name+' failed; full log retained.')

    py = sys.executable
    run('cohomology',py,HERE/'cohomology/cohomology_certificate.py',
        '--output',results/'cohomology_exact.json')
    run('yukawa',py,HERE/'yukawa/yukawa_certificate.py',
        '--output',results/'yukawa_exact.json')
    run('photon',py,HERE/'photon/check_photon.py','--repo',repo,
        '--output',results/'photon_checks.json')
    run('public_api',py,'-m','mtft.matter.collinear',
        '--output',results/'collinear_report.json')
    run('capability_index',py,repo/'scripts/build_capability_index.py','--check')
    if args.skip_tests:
        stages.append({'name':'focused_tests','status':'SKIPPED: requested'})
    else:
        run('focused_tests',py,'-m','pytest',
            repo/'tests/test_matter02_collinear.py',
            repo/'tests/test_matter01_qed.py',repo/'tests/test_matter01_normalization.py',
            repo/'tests/test_interactions.py',repo/'tests/test_pv02.py',
            repo/'tests/test_repository_docs.py','-q')
    node = shutil.which('node')
    if node:
        run('lab_logic',node,HERE/'viz/verify_matter02.cjs',
            HERE/'viz/MTFT_Matter02_Mode_Lab.html',results/'lab_runtime_checks.json')
    else:
        stages.append({'name':'lab_logic','status':'SKIPPED: Node.js unavailable'})
        print('Lab logic: SKIPPED; Node.js unavailable. Prebuilt HTML remains usable.')
    run('scientific_figure',py,HERE/'plot_matter02.py','--output',
        HERE/'viz/MTFT_Matter02_Overview.png')
    c=json.loads((results/'cohomology_exact.json').read_text())
    y=json.loads((results/'yukawa_exact.json').read_text())
    p=json.loads((results/'photon_checks.json').read_text())
    versions={}
    for name in ('numpy','sympy','pytest','matplotlib'):
        try: versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: versions[name]='unavailable'
    summary={'milestone':'Matter-02','source_version':'0.34.0 + unreleased Matter-01/02',
             'passed':True,'python':platform.python_version(),'dependencies':versions,
             'exact_certificates':{'cohomology':c['gate_count'],
                                   'yukawa':y['gates_passed'],
                                   'photon':p['exact_gates_passed']},
             'holomorphic_rank_lower_bound':3,
             'physical_mass_and_photon':'OPEN',
             'browser_layout':'UNVERIFIED: mock-DOM logic is not browser layout QA',
             'full_release_suite':'NOT RUN by this focused runner','stages':stages}
    (results/'run_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('\nMatter-02 focused gates: PASS')
    print('Exact certificates:',summary['exact_certificates'])
    print('Holomorphic rank >= 3; physical mass and photon remain OPEN.')
    print('Full results:',results)


if __name__=='__main__':
    main()
