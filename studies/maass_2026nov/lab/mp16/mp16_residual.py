#!/usr/bin/env python3
"""MP16, after A1 missed: the residuals LHS - RHS of the three odd-sector identities at level 143, as functions of the centre r0
of the test function (sigma = 0.35), and where they say odd forms are missing from MP15's list.

A form missing from the list, at R in sector chi = (e11, e13), changes identity Q by -(1/4) eps_Q(chi) h_{r0}(R).  The residual
curves are fitted with non-negative amounts of such forms (atoms on a grid of R, one family per odd sector) by NNLS.  Adding one
form to each of the four sectors at the same R changes no identity (the four sign vectors sum to zero), so the fit can only
locate forms up to that blind direction; the Hejhal searches that follow decide.

Usage: python3 -I mp16_residual.py [extra.json ...]   (extra.json: forms to add to MP15's list, [{"R":..,"e11":..,"e13":..}, ...])"""
import importlib.util
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import nnls

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("t", os.path.join(HERE, "mp16_trace.py"))
t = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t)

c_ell, c_refl = t.constants()
cos = {Q: t.Coset(143, Q, 23.5, c_ell, c_refl, cache=os.path.join(HERE, f"weights_143_{Q}.npz")) for Q in (11, 13, 143)}
arr = {Q: (np.array([r[1] for r in c.hyp]), np.array([r[2] for r in c.hyp]), np.array([r[1] for r in c.gli]), np.array([r[2] for r in c.gli]))
       for Q, c in cos.items()}
spec_ = t.odd_spectrum()
for fn in sys.argv[1:]:
    for e in json.load(open(fn)):
        spec_.setdefault((e["e11"], e["e13"]), []).append(e["R"])
EPS = {11: lambda a, b: a, 13: lambda a, b: b, 143: lambda a, b: a * b}
SIG = 0.35
r0s = np.round(np.arange(0.0, 1.2001, 0.025), 4)


def rhs(Q, r0):
    h, g = t.testfun(r0, SIG)
    lh, wh, lg, wg = arr[Q]
    c = cos[Q]
    tot = c.E * c_ell * t.int_h_sech(h) - c.L * c_refl * float(g(0.0)) + float(np.sum(g(lh) * wh)) - float(np.sum(g(lg) * wg))
    return tot / 8


def lhs(Q, r0):
    h, _ = t.testfun(r0, SIG)
    return 0.25 * sum(EPS[Q](a, b) * sum(float(h(R)) for R in Rs) for (a, b), Rs in spec_.items())


D = {Q: np.array([lhs(Q, r0) - rhs(Q, r0) for r0 in r0s]) for Q in (11, 13, 143)}
out = {"sigma": SIG, "r0": r0s.tolist(), "residual": {str(Q): D[Q].tolist() for Q in D}, "list used": {str(k): v for k, v in spec_.items()}}
for i, r0 in enumerate(r0s):
    if i % 4 == 0:
        print(f"r0 {r0:5.3f}: " + "  ".join(f"Q{Q} {D[Q][i]:+.5f}" for Q in (11, 13, 143)), flush=True)

# NNLS on r0 <= 0.9 (forms above the grid's top weigh little there)
use = r0s <= 0.9
Rgrid = np.round(np.arange(0.90, 2.3001, 0.005), 4)
cols, labels = [], []
for (a, b) in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
    for R in Rgrid:
        col = []
        for Q in (11, 13, 143):
            for r0 in r0s[use]:
                h, _ = t.testfun(float(r0), SIG)
                col.append(-0.25 * EPS[Q](a, b) * float(h(R)))
        cols.append(col)
        labels.append(((a, b), float(R)))
A = np.array(cols).T
y = np.concatenate([D[Q][use] for Q in (11, 13, 143)])
x, rn = nnls(A, y, maxiter=20000)
found = {}
for (sec, R), v in zip(labels, x):
    if v > 0.02:
        found.setdefault(str(sec), []).append((R, round(float(v), 3)))
clusters = {}
for sec, lst in found.items():
    groups, cur = [], [lst[0]]
    for R, v in lst[1:]:
        if R - cur[-1][0] <= 0.03:
            cur.append((R, v))
        else:
            groups.append(cur)
            cur = [(R, v)]
    groups.append(cur)
    clusters[sec] = [{"R mean": round(sum(R * v for R, v in g) / sum(v for _, v in g), 3), "amount": round(sum(v for _, v in g), 3)} for g in groups]
out["nnls"] = {"residual norm": float(rn), "data norm": float(np.linalg.norm(y)), "clusters (sector (e11,e13): forms)": clusters}
print(json.dumps(out["nnls"], indent=1))
tag = "_".join(os.path.splitext(os.path.basename(f))[0] for f in sys.argv[1:]) or "mp15"
with open(os.path.join(HERE, f"mp16_residual_{tag}.json"), "w") as fh:
    json.dump(out, fh, indent=1)
