#!/usr/bin/env python3
"""MP16: targeted searches for odd Maass cusp forms of Gamma0(143) that MP15's list may lack, in one sector and one window of R,
with MP15's solver unchanged (mp15_maass.setup143 and mp15_extra.coeffs / checks_at).

On a grid of step dR the two horocycles' coefficients c_2..c_9 are compared; candidates are sign changes of any difference and
local minima of their norm.  Each candidate is refined by Brent's method on one difference and kept only if the two horocycles
agree (<= 1e-6) and the Hecke relations hold (<= 1e-4), MP15's own acceptance rule; rejected brackets are poles of the system.

Usage: python3 -I mp16_find.py e11 e13 lo hi dR out.json [f1 f2]   (f1, f2: horocycle heights as fractions of Y0; MP15 used 0.96 and 0.88)"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
MP15 = os.path.join(HERE, "..", "mp15")
sys.path.insert(0, MP15)
spec = importlib.util.spec_from_file_location("x", os.path.join(MP15, "mp15_extra.py"))
x = importlib.util.module_from_spec(spec)
spec.loader.exec_module(x)
m = x.m
T0 = time.time()
e11, e13, lo, hi, dR, out = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
N, Y0, hcs, M = m.setup143()
F1, F2 = (float(sys.argv[7]), float(sys.argv[8])) if len(sys.argv) > 8 else (0.96, 0.88)
if (F1, F2) != (0.96, 0.88):                       # moved horocycles: same M and Qs, other heights; poles move, eigenvalues do not
    Qs = hcs[0].Qs if hasattr(hcs[0], "Qs") else int(math.ceil(1.08 * M / 0.88)) + 12
    hcs = (m.Horocycle(N, F1 * Y0, Qs), m.Horocycle(N, F2 * Y0, Qs))
print(f"horocycles at {F1} Y0 and {F2} Y0 (Y0 = {Y0:.6f}), M = {M}", flush=True)
sec = [(e11, e13)]


def diffs(R):
    c = x.coeffs(hcs, M, R, 1, sec)
    return c[0][0][1:9] - c[1][0][1:9]


grid = np.arange(lo, hi + 1e-12, dR)
Dg = np.array([diffs(R) for R in grid])
nrm = np.linalg.norm(Dg, axis=1)
print(f"odd e11={e11:+d} e13={e13:+d}: {len(grid)} points on [{lo}, {hi}] in {time.time() - T0:.0f} s", flush=True)
cand = set()
for i in range(len(grid) - 1):
    if any(np.sign(Dg[i, k]) != np.sign(Dg[i + 1, k]) for k in range(3)):
        cand.add(i)
for i in range(1, len(grid) - 1):
    if nrm[i] < nrm[i - 1] and nrm[i] < nrm[i + 1]:
        cand.add(i)
found, rejected = [], []
for i in sorted(cand):
    Rr = None
    for k in range(5):
        f = lambda r, k=k: diffs(r)[k]
        a, b = grid[max(i - 1, 0)], grid[min(i + 2, len(grid) - 1)]
        try:
            Rr = brentq(f, a, b, xtol=1e-12, rtol=1e-14)
            break
        except ValueError:
            continue
    if Rr is None or any(abs(Rr - e["R"]) < 1e-7 for e in found):
        continue
    c = x.coeffs(hcs, M, Rr, 1, sec)
    agree, hecke, cp = x.checks_at(c[0][0], c[1][0], (11, 13))
    if agree > 1e-6 or hecke > 1e-4:
        rejected.append({"R": Rr, "agree": agree, "hecke": hecke})
        continue
    new = all(abs(abs(v) - 1) < 1e-4 for v in cp.values())
    found.append({"sector": f"odd e11={e11:+d} e13={e13:+d}", "R": Rr, "lambda": 0.25 + Rr * Rr, "new": bool(new), "agree": agree, "hecke": hecke,
                  "c11_sqrt11": cp[11], "c13_sqrt13": cp[13], "c2": float(c[0][0][1]), "c3": float(c[0][0][2])})
    print(f"  FOUND R = {Rr:.10f}  {'new' if new else 'OLD'}  agree {agree:.1e}  hecke {hecke:.1e}  c11*sqrt11 {cp[11]:+.6f}  c13*sqrt13 {cp[13]:+.6f}", flush=True)
json.dump({"sector": [e11, e13], "window": [lo, hi], "step": dR, "M": M, "Y0": Y0, "horocycle fractions": [F1, F2], "found": found, "rejected (poles)": rejected,
           "seconds": time.time() - T0}, open(out, "w"), indent=1)
print(f"done: {len(found)} forms, {len(rejected)} rejected brackets, {time.time() - T0:.0f} s")
