"""Build the portable Matter-01 U(1) laboratory from audited JSON evidence.

Only constants and numerical checks supplied by the computation are embedded.
No network libraries or external resources are required by the output HTML.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def build(data_path: Path, output: Path | None = None) -> Path:
    here = Path(__file__).resolve().parent
    output = output or here / "Matter01_U1_Lab.html"
    evidence = json.loads(data_path.read_text())
    embedded = json.dumps(evidence, ensure_ascii=False, separators=(",", ":"))
    embedded = embedded.replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    template = (here / "lab_template.html").read_text()
    output.write_text(template.replace("__MATTER01_EVIDENCE__", embedded))
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="Audited visualization evidence JSON")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    print(build(args.data, args.output))
