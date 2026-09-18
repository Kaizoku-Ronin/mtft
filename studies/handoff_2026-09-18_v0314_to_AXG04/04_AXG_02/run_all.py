#!/usr/bin/env python3
"""Run the exact AXG-02 calculations; record reproducibility metadata."""
from pathlib import Path
import hashlib
import json
import platform
import subprocess
import sys
import numpy
import sympy
root=Path(__file__).resolve().parent
reports=[]
for name in ('parent_anomaly','flux_and_scalar'):
    subprocess.run([sys.executable,str(root/(name+'.py'))],check=True,cwd=root)
    result=json.loads((root/(name+'_results.json')).read_text())
    assert result['summary']['failed']==0
    reports.append({'study':name,**result['summary']})
out={'investigation':'AXG-02','classification':'EXACT conditional algebra and necessary constraints',
     'subaudits':reports,'total_passed':sum(r['passed'] for r in reports),'total_failed':0,
     'runtime':{'python':platform.python_version(),'sympy':sympy.__version__,'numpy':numpy.__version__},
     'input_sha256':hashlib.sha256((root/'input/smflux_v0314.py').read_bytes()).hexdigest(),
     'archive_sha256':'2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5',
     'interpretation':'Checks establish the displayed conditional algebra, not a completed parent, vacuum, or measured prediction.'}
(root/'results_summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
