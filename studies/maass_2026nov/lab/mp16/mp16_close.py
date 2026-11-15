#!/usr/bin/env python3
"""MP16: the three odd-sector identities at level 143 re-evaluated with the forms found by the targeted searches added to MP15's
list, for the test functions the plan registered (A1: r0 = 0 with sigma 0.35 and 0.45, r0 = 0.374 with sigma 0.35; A2: r0 = 0.556
and 0.845 with sigma 0.35).  The searches covered R up to 1.65 in every odd sector, so the tolerance is now the Weyl bound for the
forms above 1.65: tau' = (1/4) int_{1.65}^inf 14 r h(r) dr.
Usage: python3 -I mp16_close.py extra.json"""
import importlib.util
import json
import os
import sys

import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("t", os.path.join(HERE, "mp16_trace.py"))
t = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t)
c_ell, c_refl = t.constants()
cos = {Q: t.Coset(143, Q, 23.5, c_ell, c_refl, cache=os.path.join(HERE, f"weights_143_{Q}.npz")) for Q in (11, 13, 143)}
arr = {Q: tuple(np.array([r[i] for r in lst]) for lst in (c.hyp, c.gli) for i in (1, 2)) for Q, c in cos.items()}
spec_ = t.odd_spectrum()
extra = json.load(open(sys.argv[1]))
for e in extra:
    spec_.setdefault((e["e11"], e["e13"]), []).append(e["R"])
EPS = {11: lambda a, b: a, 13: lambda a, b: b, 143: lambda a, b: a * b}
rows = []
for sig, r0, which in ((0.35, 0.0, "A1"), (0.45, 0.0, "A1"), (0.35, 0.374, "A1"), (0.35, 0.556, "A2"), (0.35, 0.845, "A2")):
    h, g = t.testfun(r0, sig)
    tau = 0.25 * quad(lambda r: 14 * r * h(r), 1.65, 60, limit=200)[0]
    for Q in (11, 13, 143):
        lh, wh, lg, wg = arr[Q]
        c = cos[Q]
        rhs = (c.E * c_ell * t.int_h_sech(h) - c.L * c_refl * float(g(0.0)) + float(np.sum(g(lh) * wh)) - float(np.sum(g(lg) * wg))) / 8
        lhs = 0.25 * sum(EPS[Q](a, b) * sum(float(h(R)) for R in Rs) for (a, b), Rs in spec_.items())
        rows.append({"test": which, "sigma": sig, "r0": r0, "Q": Q, "LHS": lhs, "RHS": rhs, "diff": lhs - rhs, "tau'": tau, "pass": abs(lhs - rhs) <= tau + 1e-6})
        print(f"{which} sigma {sig} r0 {r0:5.3f} Q {Q:3d}: LHS {lhs:+.6f} RHS {rhs:+.6f} diff {lhs - rhs:+.2e} tau' {tau:.1e} {'PASS' if rows[-1]['pass'] else 'FAIL'}")
out = {"extra forms": extra, "list used": {str(k): sorted(v) for k, v in spec_.items()}, "rows": rows}
json.dump(out, open(os.path.join(HERE, "mp16_close.json"), "w"), indent=1)
