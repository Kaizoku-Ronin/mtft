#!/usr/bin/env python3
"""MP16: odd Maass cusp forms at a prime level p (11 or 13) in a window of R, by MP15's solver, for the oldforms they give at 143.
An odd form of level p with Fricke sign e appears at level 143 = p q twice, in the odd sectors (e_p = e, e_q = +1) and (e, -1).
Usage: python3 -I mp16_old.py p lo hi dR out.json"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("x", os.path.join(HERE, "..", "mp15", "mp15_extra.py"))
x = importlib.util.module_from_spec(spec)
spec.loader.exec_module(x)
m = x.m
T0 = time.time()
p, lo, hi, dR, out = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
Y0, hcs, M = x.horocycles(p, 0.01)
print(f"level {p}: lowest point met {Y0:.5f}, M = {M}", flush=True)
res = {"level": p, "Y0": Y0, "M": M, "window": [lo, hi], "step": dR, "found": [], "rejected": []}
grid = np.arange(lo, hi + 1e-12, dR)
for e in (1, -1):
    d = lambda r, k=1: (lambda c: c[0][0][k] - c[1][0][k])(x.coeffs(hcs, M, r, 1, [(e,)]))
    vals = np.array([[d(R, k) for k in (1, 2, 3)] for R in grid])
    for i in range(len(grid) - 1):
        if not any(np.sign(vals[i, k]) != np.sign(vals[i + 1, k]) for k in range(3)):
            continue
        Rr = None
        for k in (1, 2, 3, 4):
            try:
                Rr = brentq(lambda r: d(r, k), grid[i], grid[i + 1], xtol=1e-12, rtol=1e-14)
                break
            except ValueError:
                continue
        if Rr is None or any(abs(Rr - f["R"]) < 1e-7 for f in res["found"]):
            continue
        c = x.coeffs(hcs, M, Rr, 1, [(e,)])
        agree, hecke, cp = x.checks_at(c[0][0], c[1][0], (p,))
        rec = {"fricke": e, "R": Rr, "agree": agree, "hecke": hecke, "cp_sqrt_p": cp[p]}
        if agree > 1e-6 or hecke > 1e-4:
            res["rejected"].append(rec)
            continue
        res["found"].append(rec)
        print(f"  odd, Fricke {e:+d}: R = {Rr:.10f}  agree {agree:.1e}  hecke {hecke:.1e}  c_p sqrt p {cp[p]:+.6f}", flush=True)
json.dump(res, open(out, "w"), indent=1)
print(f"done in {time.time() - T0:.0f} s: {[(f['fricke'], round(f['R'], 6)) for f in res['found']]}")
