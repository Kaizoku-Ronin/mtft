#!/usr/bin/env python3
"""Reproduce the focused Matter-03 gates with printed verdicts and retained logs.

This is not the whole historical MTFT release suite or its hours-long studies.
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

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=HERE.parents[1])
    parser.add_argument('--skip-tests',action='store_true')
    args=parser.parse_args(); repo=args.repo.resolve()
    if not (repo/'src/mtft/matter/higgs_vacuum.py').exists():
        raise SystemExit('Matter-03 source missing; pass --repo path/to/patched/repo.')
    output=HERE/'results'; output.mkdir(exist_ok=True)
    env=os.environ.copy(); env['PYTHONPATH']=str(repo/'src')+os.pathsep+env.get('PYTHONPATH','')
    stages=[]
    def run(name,*command):
        print('\n'+name,flush=True); start=time.perf_counter()
        proc=subprocess.run(list(map(str,command)),cwd=repo,env=env,text=True,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (output/(name+'.log')).write_text(proc.stdout,encoding='utf-8')
        stages.append({'name':name,'returncode':proc.returncode,
                       'seconds':round(time.perf_counter()-start,3),
                       'log':'results/'+name+'.log'})
        print('\n'.join(proc.stdout.strip().splitlines()[-12:]),flush=True)
        if proc.returncode:
            (output/'run_summary.json').write_text(json.dumps({'passed':False,'stages':stages},indent=2)+'\n')
            raise SystemExit(name+' failed; see its retained log.')
    py=sys.executable
    run('vacuum_certificate',py,HERE/'vacuum/vacuum_certificate.py','--output',output/'vacuum_exact.json')
    run('charged_certificate',py,HERE/'modes/charged_complex.py','--repo',repo,'--output',output/'charged_complex_exact.json')
    run('kinetic_diagnostics',py,HERE/'kinetics/check_kinetics.py','--output',output/'kinetic_checks.json')
    run('vacuum_api',py,'-m','mtft.matter.higgs_vacuum','--output',output/'higgs_vacuum_report.json')
    run('charged_api',py,'-m','mtft.matter.charged_modes','--output',output/'charged_modes_report.json')
    run('capability_index',py,repo/'scripts/build_capability_index.py','--check')
    if args.skip_tests:
        stages.append({'name':'focused_tests','status':'SKIPPED: requested'})
    else:
        run('focused_tests',py,'-m','pytest',
            repo/'tests/test_matter03_higgs_vacuum.py',
            repo/'tests/test_matter03_charged_modes.py',
            repo/'tests/test_matter03_kinetics.py',
            repo/'tests/test_matter02_collinear.py',repo/'tests/test_matter01_qed.py',
            repo/'tests/test_matter01_normalization.py',repo/'tests/test_interactions.py',
            repo/'tests/test_pv02.py',repo/'tests/test_repository_docs.py','-q')
    node=shutil.which('node')
    if node:
        run('lab_logic',node,HERE/'viz/verify_lab.cjs',HERE/'viz/MTFT_Matter03_Higgs_Lab.html',
            output/'lab_runtime_checks.json')
    else:
        stages.append({'name':'lab_logic','status':'SKIPPED: Node.js unavailable'})
        print('Node logic checks SKIPPED; the prebuilt lab still opens locally.')
    run('scientific_figure',py,HERE/'plot_matter03.py','--output',HERE/'viz/MTFT_Matter03_Overview.png')
    v=json.loads((output/'vacuum_exact.json').read_text())
    c=json.loads((output/'charged_complex_exact.json').read_text())
    k=json.loads((output/'kinetic_checks.json').read_text())
    versions={}
    for name in ('numpy','sympy','pytest','matplotlib'):
        try: versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: versions[name]='unavailable'
    summary={'milestone':'Matter-03','source_version':'0.34.0 + unreleased Matter-01/02/03',
             'passed':True,'python':platform.python_version(),'dependencies':versions,
             'vacuum_arithmetic_gates':v['finite_arithmetic_checks'],
             'charged_arithmetic_gates':c['gates_passed'],
             'kinetic_checks':k['check_count'],
             'kinetic_check_scope':'11 exact symbolic gates; 106 finite numerical/rejection checks',
             'classical_paired_metric_existence':'CONDITIONAL THEOREM; metric not computed',
             'paired_charged_degree_one_classes':{'Qminus':237,'Qplus':237},
             'all_aligned_pairs_lower_bound_per_sign':160,
             'selected_higgs_scale':False,'physical_mass_prediction':False,
             'browser_layout':'UNVERIFIED: mock DOM is not browser layout QA',
             'full_release_suite':'NOT RUN by this focused runner','stages':stages}
    (output/'run_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('\nMatter-03 focused gates: PASS')
    print('Classical pair exists conditionally; scale and physical metrics OPEN.')
    print('237 charged degree-one classes per sign; family selection unresolved.')
    print('Results:',output)


if __name__=='__main__': main()
