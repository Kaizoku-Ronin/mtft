#!/usr/bin/env python3
"""MP16: the odd Fricke identity at the prime levels 11 and 13, whose odd spectra below R = 2 are known completely
(level 11: none, by mp16_old.py and LMFDB; level 13: only 13.1, R = 0.9708154170, Fricke -1).  A test of this code's
geometric side on spectra that MP15 did not compute:  S_odd,F+ - S_odd,F- = (1/2)(T_{W_N} - T_{R W_N}).
Usage: python3 -I mp16_prime.py"""
import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("t", os.path.join(HERE, "mp16_trace.py"))
t = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t)
c_ell, c_refl = t.constants()
known = {11: [], 13: [(0.9708154170, -1)]}
out = {}
for N in (11, 13):
    c = t.Coset(N, N, 23.5, c_ell, c_refl, cache=os.path.join(HERE, f"weights_{N}_{N}.npz"))
    lh, wh = np.array([r[1] for r in c.hyp]), np.array([r[2] for r in c.hyp])
    lg, wg = np.array([r[1] for r in c.gli]), np.array([r[2] for r in c.gli])
    rows = []
    for sig in (0.35, 0.45):
        for r0 in (0.0, 0.3, 0.6, 0.97):
            h, g = t.testfun(r0, sig)
            ell = c.E * c_ell * t.int_h_sech(h)
            mir = c.L * c_refl * float(g(0.0))
            rhs = 0.5 * (ell - mir + float(np.sum(g(lh) * wh)) - float(np.sum(g(lg) * wg)))
            lhs = sum(e * float(h(R)) for R, e in known[N])
            rows.append({"sigma": sig, "r0": r0, "LHS": lhs, "RHS": rhs, "diff": lhs - rhs})
            print(f"level {N}: sigma {sig} r0 {r0:4.2f}: LHS {lhs:+.8f} RHS {rhs:+.8f} diff {lhs - rhs:+.2e}", flush=True)
    out[str(N)] = {"fixed points": c.E, "mirror length": c.L, "rows": rows}
json.dump(out, open(os.path.join(HERE, "mp16_prime.json"), "w"), indent=1)
