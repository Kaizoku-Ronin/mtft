"""Replay the bounded audit against a supplied v0.26.2 source tree."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--source',type=Path,required=True)
args=parser.parse_args()
source=args.source.resolve()
env=os.environ.copy()
env['PYTHONPATH']=str(source/'src')+os.pathsep+env.get('PYTHONPATH','')
env.setdefault('OPENBLAS_NUM_THREADS','1')
jobs=[['release_mellin_accounting.py'],['d4_check.py','--source',str(source)],
      ['bimodule_control.py','--source',str(source),'--verify-theta'],
      ['finite_graph_check.py'],['render_finite_controls.py']]
log=[]
for job in jobs:
    result=subprocess.run([sys.executable,str(ROOT/job[0]),*job[1:]],cwd=ROOT,env=env,text=True,capture_output=True)
    log.append('Command: '+' '.join(job)+'\n'+result.stdout+result.stderr+'\nExit code: '+str(result.returncode)+'\n')
    (ROOT/'run_log.txt').write_text('\n'.join(log))
    print(job[0]+': '+('passed' if result.returncode==0 else 'FAILED'),flush=True)
    if result.returncode:
        raise SystemExit(result.returncode)
tests=['tests/test_cc21_cc22.py','tests/test_release_pin.py']
result=subprocess.run([sys.executable,'-m','pytest','-q',*tests],cwd=source,env=env,text=True,capture_output=True)
(ROOT/'targeted_tests.txt').write_text('Command: python -m pytest -q '+' '.join(tests)+'\n'+result.stdout+result.stderr)
print(result.stdout,end='')
raise SystemExit(result.returncode)
