import argparse, json, math
from . import gates, riemann_matrix, cp_channel_report, bergman_density

def main():
    ap = argparse.ArgumentParser(prog="python -m mtft.periods")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("verify"); sub.add_parser("tau"); sub.add_parser("physics")
    sub.add_parser("involutions"); sub.add_parser("oldtorus")
    h = sub.add_parser("hamiltonian"); h.add_argument("--potential", default="width")
    sub.add_parser("crossover")
    b = sub.add_parser("bergman"); b.add_argument("--x", type=float, default=0.0)
    b.add_argument("--y", type=float, default=1/math.sqrt(143))
    a = ap.parse_args()
    if a.cmd == "verify":
        out = gates.run_all()
        for k, v in out.items():
            print(f"[PASS] {k}")
        print(f"{len(out)}/{len(out)} period/Hodge gates PASS")
    elif a.cmd == "tau":
        t = riemann_matrix(50)
        print(t)
    elif a.cmd == "physics":
        print(json.dumps(cp_channel_report("width"), indent=1, default=str))
    elif a.cmd == "bergman":
        print(bergman_density(complex(a.x, a.y)))
    elif a.cmd == "involutions":
        from . import involutions as IV
        print(json.dumps({"eps": IV.al_signs()["eps"],
                          "census": IV.sector_census(),
                          "route2": IV.route2_fixed_intersections()},
                         default=str, indent=1))
    elif a.cmd == "oldtorus":
        from . import oldtorus as OT
        print(json.dumps({"polarization": OT.polarization_type(),
                          "l9_index": OT.l9_index(),
                          "product": OT.product_charpoly(),
                          "entropy": OT.entropy()}, default=str, indent=1))
    elif a.cmd == "hamiltonian":
        from . import hamiltonian as HM
        print(json.dumps({"report": HM.channel_report(a.potential),
                          "rho": HM.pairing_stability(a.potential),
                          "oldspace": HM.oldspace_routing(a.potential)},
                         default=str, indent=1))
    elif a.cmd == "crossover":
        from . import channels as CH
        print(CH.mode_crossover(4, 1))

main()
