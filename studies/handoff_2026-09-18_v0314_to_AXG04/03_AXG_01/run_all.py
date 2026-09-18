#!/usr/bin/env python3
"""Run the three AXG-01 exact subaudits and record their combined status."""
from pathlib import Path
import hashlib
import json
import platform
import subprocess
import sys
import numpy
import sympy

root = Path(__file__).resolve().parent
reports = []
for name in ("anomaly", "operator", "vacuum"):
    subprocess.run([sys.executable, str(root / (name + "_audit.py"))], check=True, cwd=root)
    result = json.loads((root / (name + "_results.json")).read_text())
    count = result["check_count"] if name == "anomaly" else result["summary"]["passed"]
    passed = result["all_passed"] if name == "anomaly" else result["summary"]["failed"] == 0
    if not passed:
        raise RuntimeError(name + " failed")
    reports.append({"subaudit": name, "passed": count, "failed": 0})
result = {
    "investigation": "AXG-01",
    "status": "EXACT symbolic checks; conditional effective-theory interpretation",
    "subaudits": reports,
    "total_passed": sum(r["passed"] for r in reports),
    "total_failed": 0,
    "runtime": {"python": platform.python_version(), "sympy": sympy.__version__, "numpy": numpy.__version__},
    "copied_source_sha256": hashlib.sha256((root / "input/smflux_v0314.py").read_bytes()).hexdigest(),
    "meaning": "Algebraic assertions for specified inputs; not experimental validation or a completed parent theory.",
}
(root / "results_summary.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
