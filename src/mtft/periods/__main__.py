import argparse, json, math
from . import gates, riemann_matrix, cp_channel_report, bergman_density

def main():
    ap = argparse.ArgumentParser(prog="python -m mtft.periods")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("verify"); sub.add_parser("tau"); sub.add_parser("physics")
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

main()
