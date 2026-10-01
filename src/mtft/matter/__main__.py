"""Run Matter-01's measured-input QED benchmark from a terminal."""

import argparse
import json
from pathlib import Path

from .qed import ALPHA_REFERENCE, ELECTRON_MASS_GEV
from .report import benchmark, benchmark_report


def main(argv=None):
    parser = argparse.ArgumentParser(description="Matter-01 massive QED reference; not MTFT predictions")
    sub = parser.add_subparsers(dest="command", required=True)
    point = sub.add_parser("benchmark", help="one CM point and its two-diagram amplitudes")
    point.add_argument("--sqrt-s-gev", type=float, default=0.01)
    point.add_argument("--cos-theta", type=float, default=0.37)
    point.add_argument("--phi", type=float, default=0.0)
    point.add_argument("--mass-gev", type=float, default=ELECTRON_MASS_GEV)
    point.add_argument("--alpha", type=float, default=ALPHA_REFERENCE)
    point.add_argument("--output", type=Path)
    report = sub.add_parser("report", help="fixed reference grid and numerical structural gates")
    report.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        result = (benchmark(args.sqrt_s_gev, args.cos_theta, mass_gev=args.mass_gev,
                            alpha=args.alpha, phi=args.phi) if args.command == "benchmark"
                  else benchmark_report())
    except (ValueError, TypeError) as exc:
        parser.error(str(exc))
    serialized = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
        print(str(args.output))
    else:
        print(serialized, end="")
    if args.command == "report" and not result["verification"]["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
