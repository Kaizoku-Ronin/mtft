#!/usr/bin/env python3
"""Run the six AXG-04 audits; no network or external package checkout needed."""
from pathlib import Path
import json
import platform
import subprocess
import sys
import sympy

ROOT=Path(__file__).resolve().parent
names=['chiral_parent','charged_action','higgs_action','gravity_background',
       'larger_parent','global_bordism']
summary={'investigation':'AXG-04','python':platform.python_version(),
         'sympy':sympy.__version__,'audits':{},'passed':0,'failed':0}
for name in names:
    subprocess.run([sys.executable,str(ROOT/(name+'_audit.py'))],cwd=ROOT,check=True)
    data=json.loads((ROOT/(name+'_results.json')).read_text())
    count=data.get('summary',{'passed':data.get('checks_passed',0),
                              'failed':data.get('checks_failed',0)})
    assert count['passed']==len(data['checks']) and count['failed']==0
    summary['audits'][name]=count
    summary['passed']+=count['passed'];summary['failed']+=count['failed']
summary['meaning']='Exact conditional algebra and analytic-consequence checks, not experimental validation, a global GS construction, UV completion, or a realistic vacuum.'
(ROOT/'results_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
