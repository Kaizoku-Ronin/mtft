"""Run companion experiments with an importable patched MTFT source tree."""
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
environment = os.environ.copy()
environment.setdefault('OPENBLAS_NUM_THREADS','1')
for script in ['check_corrections.py','charge_study.py','local_qed_covariance.py',
               'photon_atomic.py','magnetic_and_coupling.py','render_results.py']:
    print('Running',script,flush=True)
    subprocess.run([sys.executable,str(ROOT/script)],cwd=ROOT,env=environment,check=True)
