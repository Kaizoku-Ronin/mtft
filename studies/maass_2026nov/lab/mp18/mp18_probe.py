"""MP18: test given (subsector, R) pairs directly on several horocycle pairs, and locate missing forms with Route T.
Usage: python3 -I mp18_probe.py at                 the oldform partners and the forms only pair A found
       python3 -I mp18_probe.py locate out.json L1.json ...   Route T's residual against the list, scanned in r0"""
import json, math, os, sys, importlib.util
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("tw", os.path.join(HERE, "mp18_twisted.py")); tw = importlib.util.module_from_spec(spec); spec.loader.exec_module(tw)

def at():
    pairs = {"A": (0.96, 0.88), "B": (0.93, 0.85), "C": (0.90, 0.82)}
    H = {k: tw.setup(*v) for k, v in pairs.items()}
    tests = [((0, (1, 1)), 1.0548437), ((1, (1, 1)), 1.1532812), ((0, (1, -1)), 1.1651519), ((1, (1, -1)), 1.3763584),
             ((0, (-1, 1)), 1.5243501), ((1, (1, 1)), 1.7632399)]
    out = []
    for (par, eps), R0 in tests:
        row = {"subsector": f"{'even' if par == 0 else 'odd'} e11={eps[0]:+d} e13={eps[1]:+d}", "R0": R0}
        for k, (Y0, M, hcs) in H.items():
            fr = np.arange(1, M + 1) - 0.5
            F = lambda r: (lambda cs: cs[0][1] - cs[1][1])(tw.coeffs(r, hcs, fr, par, eps))
            try:
                R = brentq(F, R0 - 3e-4, R0 + 3e-4, xtol=1e-12, rtol=1e-14)
                c1, c2 = tw.coeffs(R, hcs, fr, par, eps)
                ag, he, b11, b13 = tw.tests(c1, c2)
                row[k] = {"R": R, "agree": ag, "hecke": he, "b11_sqrt11": b11, "b13_sqrt13": b13}
            except ValueError:
                c1, c2 = tw.coeffs(R0, hcs, fr, par, eps)
                row[k] = {"no sign change": True, "agree at R0": tw.tests(c1, c2)[0]}
        out.append(row)
        print(json.dumps(row), flush=True)
    json.dump(out, open(os.path.join(HERE, "probe_at.json"), "w"), indent=1)

def locate(outf, lists, sig=0.3):
    t16 = tw._t16()
    W = np.load(os.path.join(HERE, "weights_tw.npz"))["hyp"]
    Rl = sorted({(e["subsector"], round(e["R"], 7)): e["R"] for f in lists for e in json.load(open(f))["found"]}.values())   # keyed by subsector: oldform pairs share R
    rows = []
    for r0 in np.arange(0.0, 1.81, 0.02):
        h, g = t16.testfun(float(r0), sig)
        geo = tw.geometric(h, g, W)["total"]
        spec_ = float(sum(h(R) for R in Rl))
        rows.append((float(r0), geo, spec_, geo - spec_))
    json.dump({"sigma": sig, "forms": len(Rl), "rows": rows}, open(outf, "w"))
    for r0, geo, sp, d in rows:
        if r0 <= 1.6:
            print(f"  r0 {r0:.2f}  geometry {geo:9.5f}  list {sp:9.5f}  missing {d:+.4f}")

if __name__ == "__main__":
    if sys.argv[1] == "at":
        at()
    else:
        locate(sys.argv[2], sys.argv[3:])
