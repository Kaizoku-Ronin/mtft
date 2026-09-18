#!/usr/bin/env python3
"""Run the four bounded AXG-03 audits and summarize exact checks."""
from pathlib import Path
import json
import platform
import subprocess
import sys
import sympy

ROOT = Path(__file__).resolve().parent
audits = ['tensor_factor', 'tensor_flux', 'higgs_parent', 'discrete_operator']
summary = {'investigation': 'AXG-03', 'python': platform.python_version(),
           'sympy': sympy.__version__, 'audits': {}, 'passed': 0, 'failed': 0}
for name in audits:
    subprocess.run([sys.executable, str(ROOT / (name+'_audit.py'))],
                   cwd=ROOT, check=True)
    result = json.loads((ROOT / (name+'_results.json')).read_text())
    counts = result.get('summary', {'passed': result.get('checks_passed', 0),
                                   'failed': result.get('checks_failed', 0)})
    assert counts['passed'] == len(result['checks']) and counts['failed'] == 0
    summary['audits'][name] = counts
    summary['passed'] += counts['passed']
    summary['failed'] += counts['failed']
summary['meaning'] = ('Exact arithmetic and analytic-consequence checks under the stated '
                      'assumptions; not experimental validation or a complete parent theory.')
(ROOT/'results_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
