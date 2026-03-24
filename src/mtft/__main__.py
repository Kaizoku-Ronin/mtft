"""
MTFT Command-Line Interface
============================

Usage:
    python -m mtft verify          Run all predictions, print pass/fail
    python -m mtft report          Full falsification report
    python -m mtft tower [N]       Multi-N tower analysis (default N=15)
    python -m mtft screen          Materials screening table
    python -m mtft info            Version, constants, module count

Examples:
    $ python -m mtft verify
    MTFT v0.6.0 — 21/23 predictions pass (σ < 2)

    $ python -m mtft tower 20
    [prints full SU(2)→SU(20) landscape]

    $ python -m mtft screen
    [prints Tano metric T_c predictions vs experiment]
"""

from __future__ import annotations

import sys
import argparse


def cmd_verify(args):
    from mtft.falsify import falsification_test
    result = falsification_test()
    import mtft
    print(f"MTFT v{mtft.__version__} — {result['passed']}/{result['total']} "
          f"predictions pass (σ < 2)")
    if result['tension']:
        print(f"  {result['tension']} in tension (2-3σ)")
    if result['failed']:
        print(f"  {result['failed']} beyond 3σ")
    for p in result['predictions']:
        sym = "✓" if p.status == "PASS" else ("~" if p.status == "TENSION" else "✗")
        print(f"  {sym} #{p.number:2d} {p.relation:42s} err={p.error_percent:.4f}%")


def cmd_report(args):
    from mtft.falsify import report
    report(verbose=True)


def cmd_tower(args):
    from mtft.tower import tower_report
    tower_report(N_max=args.N, y=args.y, verbose=True)


def cmd_screen(args):
    from mtft.tano_metric import materials_screening
    results = materials_screening()
    print("MTFT MATERIALS SCREENING (Tano Metric + Geometry)")
    print("=" * 70)
    print(f"  {'Material':>10s} {'ΔT':>8s} {'G':>5s} {'λ_eff':>7s} "
          f"{'T_c pred':>9s} {'T_c obs':>8s} {'Error':>7s}")
    print(f"  {'-' * 58}")
    names = ["LaH₁₀", "H₃S", "YBCO", "MgB₂", "Al", "Pb"]
    for r, name in zip(results, names):
        obs = r.get('observed_Tc_K', 0)
        err = r.get('error_percent', 0)
        print(f"  {name:>10s} {r['delta_T']:8.1f} {r['G']:5.1f} "
              f"{r['lambda_eff']:7.3f} {r['T_c_K']:9.1f} K "
              f"{obs:8.1f} K {err:6.1f}%")


def cmd_info(args):
    import mtft
    from mtft.constants import (
        FEIGENBAUM_DELTA, T_INF, EULER_GAMMA, LAMBERT_OMEGA,
        XI, GAUGE,
    )
    mods = [m for m in dir(mtft) if not m.startswith('_')]
    print(f"MTFT v{mtft.__version__}")
    print(f"  Modules exported: {len(mods)}")
    print(f"  δ  = {FEIGENBAUM_DELTA:.15f}")
    print(f"  T∞ = {T_INF:.15f}")
    print(f"  γ  = {EULER_GAMMA:.15f}")
    print(f"  Ω  = {LAMBERT_OMEGA:.15f}")
    print(f"  ξ  = {XI:.15f}")
    print(f"  α⁻¹ = {1/GAUGE.alpha:.10f}")
    print(f"  sin²θ_W = {GAUGE.sin2_theta_W:.10f} = 3/13")


def main():
    parser = argparse.ArgumentParser(
        prog="mtft",
        description="MTFT — Modular Time Field Theory toolkit",
    )
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("verify", help="Run all predictions")
    sub.add_parser("report", help="Full falsification report")

    p_tower = sub.add_parser("tower", help="Multi-N tower analysis")
    p_tower.add_argument("N", type=int, nargs="?", default=15)
    p_tower.add_argument("--y", type=float, default=0.10)

    sub.add_parser("screen", help="Materials screening")
    sub.add_parser("info", help="Package info")

    args = parser.parse_args()

    commands = {
        "verify": cmd_verify,
        "report": cmd_report,
        "tower": cmd_tower,
        "screen": cmd_screen,
        "info": cmd_info,
    }

    if args.command in commands:
        commands[args.command](args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
