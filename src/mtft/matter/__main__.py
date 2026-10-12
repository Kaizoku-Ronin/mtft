"""Matter command line: Matter-01's measured-input QED benchmark, and the matter-programme toolkit (presets, background, ledger, epsilon, puzzle, levels, timerev, prereg, forces, atlas)."""

import argparse
import json
from pathlib import Path

from .qed import ALPHA_REFERENCE, ELECTRON_MASS_GEV
from .report import benchmark, benchmark_report


def _toolkit(args):
    if args.command == "presets":
        from .background import PRESETS
        for k, b in PRESETS.items(): print(f"{k:20s} {b.note}")
        return 0
    if args.command == "background":
        from .background import PRESETS
        if args.name not in PRESETS:
            print(f"unknown preset {args.name!r}; choose from {', '.join(PRESETS)}"); return 2
        b = PRESETS[args.name]
        print(json.dumps(b.analyse(), indent=2, default=str) if args.json else b.text()); return 0
    if args.command == "ledger":
        from . import ledger as L
        if args.markdown: print(L.to_markdown(), end=""); return 0
        for e in L.entries(status=args.status):
            print(f"{e['id']:7s} {e['wave']:10s} {e['class']:10s} {e['status']:10s} {e['statement']}" + (f"  [{e['note']}]" if e["note"] else ""))
        return 0
    if args.command == "puzzle":
        from . import puzzle as Z
        if args.json: print(Z.json_text()); return 0
        if args.markdown: print(Z.to_markdown(), end=""); return 0
        if args.dot: print(Z.to_dot(), end=""); return 0
        if args.piece and args.piece not in Z.pieces():
            print(f"unknown piece {args.piece!r}; choose from {', '.join(Z.ORDER)}"); return 2
        print(Z.text(args.piece, only_open=args.open)); return 0
    if args.command == "timerev":
        from . import time_reversal as T
        if args.candidate:
            kind, value = args.candidate
            try:
                if kind == "cm":
                    cand = ("cm", tuple(int(x) for x in value.split(",")))
                elif kind == "cusp":
                    from fractions import Fraction
                    cand = ("cusp", None if value in ("inf", "oo") else Fraction(value))
                elif kind == "numeric":
                    if args.tolerance is None:
                        print("a numerical candidate needs --tolerance"); return 2
                    cand = ("numeric", complex(value))
                else:
                    print("candidate kind is cm, cusp or numeric"); return 2
                r = T.strong_cp_test(cand, det_phase=args.det_phase, tolerance=args.tolerance)
            except ValueError as exc:
                print(exc); return 2
            print(json.dumps(r)); return 0
        if args.json: print(json.dumps(T.report(), default=str)); return 0
        print(T.text()); return 0
    if args.command == "atlas":
        from . import vacuum_atlas as VA
        if args.json: print(json.dumps(VA.report(), default=str)); return 0
        print(VA.text()); return 0
    if args.command == "forces":
        from . import forces as F
        if args.touch:
            print(json.dumps(F.touch_scan(args.touch), default=str)); return 0
        if args.column:
            print(json.dumps(F.exact_sm_solutions(args.column, box=args.box), default=str)); return 0
        if args.json: print(json.dumps(F.report(box=args.box), default=str)); return 0
        print(F.text(box=args.box)); return 0
    if args.command == "prereg":
        from . import prereg as P
        if args.verify:
            v = P.verify()
            for k, x in v.items(): print(k, x)
            return 0 if all(x is True or x == "superseded" for x in v.values()) else 1
        if args.json: print(json.dumps(P.report(), default=str)); return 0
        print(P.text()); return 0
    if args.command == "levels":
        from . import star_quotients as S
        if args.json: print(json.dumps(S.report(), default=str)); return 0
        if args.level is not None:
            try:
                print(S.text(args.level))
            except ValueError as exc:
                print(exc); return 2
            return 0
        print(S.text()); return 0
    if args.command == "phases":
        from .condensate_phases import report as cr
        r = cr(); i = r["invariants"]; g = r["gj"]
        print(f"spin-1 quartic: E|s|^4 = {i['a']:.4f} |c|^4 {i['b']:+.4f} |theta|^2  (coherent maximal, polar minimal: tree level prefers polar)")
        print(f"down-quark rho = D2/D1^2 = {g['rho_down']:.3f} (orthonormal-coherent = 0.5)")
        print(f"lepton sign flip: m_mu/m_s {g['pred_mmu_over_ms']:.2f} (data {g['data_mmu_over_ms']:.2f}), m_e/m_d {g['pred_me_over_md']:.3f} (data {g['data_me_over_md']:.3f})")
        return 0
    if args.command == "leptons":
        from .leptons import report as lr
        r = lr(); g = r["gj"]
        print(f"at {g['mu_star_GeV']:.3g} GeV: m_mu/m_s = {g['mmu_over_ms']:.2f}, m_e/m_d = {g['me_over_md']:.3f}; eps_e/eps_d needed {g['eps_ratio_second_family']:.2f} (2nd) vs {g['eps_ratio_first_family']:.2f} (1st)")
        print("Koide Q:", {k: round(x, 6) for k, x in r["koide"].items()})
        print("equal-overlap alignments:", {k: tuple(round(x, 4) for x in s) for k, s in r["fits"].items()})
        return 0
    if args.command == "ratio3":
        from .ratio_three import report as rt
        r = rt(); print(f"leading eps_d/eps_u = {r['leading_ratio']:.0f}  (data {r['targets']['eps_d'] / r['targets']['eps_u']:.2f})")
        for K, w in r["window"].items():
            print(f"  K = {K}: r = {w['r_from_eps_u']:.0f} (from eps_u) -> ratio {w['ratio_at_r_u']:.3f};  r = {w['r_from_eps_d']:.0f} (from eps_d) -> ratio {w['ratio_at_r_d']:.3f}")
        print(f"jet-sector light neutrino (third family): {1e3 * r['nu_light_GeV'][1.0]:.1f} MeV, independent of the D-term suppression")
        return 0
    if args.command == "scan":
        from .flux_scan import families, higgs_charge_rule
        for k, s in families().items():
            ch = s["charges"]
            print(f"{k}:\n  a = {s['a']}  gamma = {s['g']}  eps = {s['eps']} A_E  kappa = {s['kappa']}  condensate {s['pair']}  nu^c {ch['nu_c']}  Higgs (up, down) ({ch['delta_u']}, {ch['delta_d']})  zero modes {s['fermions']}")
        print("eps_d/eps_u admissible with SU(4) kept (smallest charges):", ", ".join(f"{k} {v}" for k, v in sorted(higgs_charge_rule().items()) if 1 < k <= 5))
        return 0
    from .epsilon import report as eps_report
    r = eps_report()
    print("rank-one entries for a torus-constant D-term:", r["rank_one_entries"])
    print(f"eps_d/eps_u needed by the data: {r['eps_ratio_needed']:.3f}")
    print(f"Arakelov density at the cusp: {r['rho_Ar_cusp']:.5f};  phi' = {r['phi_prime_coeff']:.4f}/sqrt(r)")
    print(f"eps_gamma at r = 176: {r['eps_gamma_r176'][0]:.3f} - {r['eps_gamma_r176'][1]:.3f};  eps_d = 0.14 needs r = {r['r_for_eps_d'][0]:.0f} - {r['r_for_eps_d'][1]:.0f}")
    return 0


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
    sub.add_parser("presets", help="list the toolkit's background presets")
    bg = sub.add_parser("background", help="analyse a background preset (Matter-04 ... EPS-01)")
    bg.add_argument("name", nargs="?", default="matter04")
    bg.add_argument("--json", action="store_true")
    led = sub.add_parser("ledger", help="the matter-programme ledger")
    led.add_argument("--status", default=None, choices=["active", "caveat", "superseded", "retracted"])
    led.add_argument("--markdown", action="store_true")
    sub.add_parser("epsilon", help="EPS-01: what sets the family-hierarchy parameter")
    sub.add_parser("scan", help="FS-01: stack-level flux scan (curve degrees and gamma-charges)")
    sub.add_parser("ratio3", help="RT-01: the ratio-3 background — eps_d/eps_u and the jet-sector neutrinos")
    sub.add_parser("leptons", help="LEP-01: the charged-lepton texture (Georgi–Jarlskog, Koide)")
    sub.add_parser("phases", help="CP-01: spin-1 phases of the condensate alignment and the sign-flip test")
    pz = sub.add_parser("puzzle", help="the programme map: every piece, its status and its own open puzzles")
    pz.add_argument("--piece", default=None, help="one piece (I ... X) with its question, verdict and sub-puzzles")
    pz.add_argument("--open", action="store_true", help="only the open sub-puzzles, with what would settle each")
    fmt = pz.add_mutually_exclusive_group()
    fmt.add_argument("--json", action="store_true", help="machine-readable export")
    fmt.add_argument("--markdown", action="store_true", help="the generated docs/SM/MATTER_PUZZLE_MAP.md")
    fmt.add_argument("--dot", action="store_true", help="the dependency graph in Graphviz DOT")
    lv = sub.add_parser("levels", help="SQ-01: the star quotients X0(pq)* for supersingular p < q — is X0(143) special?")
    lv.add_argument("--level", type=int, default=None, help="one level N = pq: its shape, what is unique to it, its siblings")
    lv.add_argument("--json", action="store_true", help="machine-readable census summary")
    tr = sub.add_parser("timerev", help="TR-01: time reversal — the two mirror classes, the thirteen star-quotient curves, the sealed strong-CP test")
    tr.add_argument("--json", action="store_true", help="machine-readable report")
    tr.add_argument("--candidate", nargs=2, metavar=("KIND", "VALUE"), help="run the sealed I.4 test on a candidate: cm A,B,C | cusp p/q | numeric x+yj")
    tr.add_argument("--det-phase", type=float, default=None, help="arg det(M_u M_d) in radians, if known")
    tr.add_argument("--tolerance", type=float, default=None, help="tolerance for a numerical candidate")
    pr = sub.add_parser("prereg", help="PR-01: sealed pre-registrations and the trials counter for 'is 143 special?'")
    pr.add_argument("--json", action="store_true", help="machine-readable report")
    pr.add_argument("--verify", action="store_true", help="recompute every seal from the code as it stands; exit 1 if any changed")
    fo = sub.add_parser("forces", help="GF-01: forces first — the eight dots, the touch scan, charge quantisation, hypercharge, the torus, the parent columns, the clock of twelve")
    fo.add_argument("--json", action="store_true", help="machine-readable report")
    fo.add_argument("--touch", choices=("U8", "SO16", "E8"), default=None, help="the touch scan of one parent's menu (sealed rule PR-S3)")
    fo.add_argument("--column", choices=("U8", "SO16", "E8"), default=None, help="the exact-SM chirality system of one parent column")
    fo.add_argument("--box", type=int, default=4, help="search box for the integer twist indices phi (default 4)")
    at = sub.add_parser("atlas", help="PG-01: the vacuum as a dynamical system — the critical-point atlas, Hubble friction, the toy run, the simulation gate")
    at.add_argument("--json", action="store_true", help="machine-readable report")
    args = parser.parse_args(argv)
    if args.command in ("presets", "background", "ledger", "epsilon", "scan", "ratio3", "leptons", "phases", "puzzle", "levels", "timerev", "prereg", "forces", "atlas"):
        try:
            return _toolkit(args)
        except BrokenPipeError:                       # output piped into head/less and closed early
            import os, sys
            sys.stdout = open(os.devnull, "w"); return 0
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
