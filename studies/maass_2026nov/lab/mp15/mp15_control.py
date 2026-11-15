#!/usr/bin/env python3
"""MP15 control (the plan's amendment of 13:58): the lowest odd Maass cusp forms with eps_N = +1 at the squarefree levels
133 = 7*19 and 145 = 5*29, by the same solver as mp15_maass.py.  Usage: python3 mp15_control.py N R0 R1 dR out.json"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq

spec = importlib.util.spec_from_file_location("m", os.path.join(os.path.dirname(os.path.abspath(__file__)), "mp15_maass.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
T0 = time.time()
N, R0, R1, dR, out = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
ps = m.primes_of(N)
probe = [m.reduce_point(float(x), 0.002, N)[1] for x in np.linspace(0.0001, 0.4999, 4000)]
Y0 = min(probe)
while True:                                           # the probe can miss the lowest part of the domain: lower until every point is raised
    M = int(math.ceil(40.0 / (2 * math.pi * Y0)))
    Qs = int(math.ceil(1.08 * M / 0.88)) + 12
    try:
        hcs = (m.Horocycle(N, 0.96 * Y0, Qs), m.Horocycle(N, 0.88 * Y0, Qs))
        break
    except AssertionError:
        Y0 *= 0.9
SECT = [(1, 1, 1), (1, -1, -1)]                          # odd, eps_p, eps_q with eps_N = +1
l = np.arange(1, M + 1)
print(f"level {N} = {ps[0]}*{ps[1]}: lowest point met {Y0:.6f}, M = {M}, Qs = {Qs}", flush=True)


def coeffs(R):
    res = []
    for hc in hcs:
        Kc = (m.Kir(R, 2 * np.pi * np.outer(hc.ys, l)), m.Kir(R, 2 * np.pi * l * hc.Y))
        res.append([c for c, _ in m.solve_coeffs(R, hc, M, 1, [hc.chi((a, b)) for (_, a, b) in SECT], Kc)])
    return res                                            # res[h][sector] = c vector


Rs = np.arange(R0, R1 + 1e-12, dR)
D = np.zeros((len(Rs), len(SECT), 8))
for i, R in enumerate(Rs):
    cs = coeffs(R)
    for s in range(len(SECT)):
        D[i, s] = (cs[0][s] - cs[1][s])[1:9]
    if i % 80 == 0:
        print(f"  R = {R:.4f} ({i + 1}/{len(Rs)}, {time.time() - T0:.0f} s)", flush=True)

found = []
for s, (par, a, b) in enumerate(SECT):
    nrm = np.linalg.norm(D[:, s], axis=1)
    med = np.median(nrm)
    cand = {i for i in range(1, len(Rs) - 1) if nrm[i] < nrm[i - 1] and nrm[i] < nrm[i + 1] and nrm[i] < 0.2 * med}
    cand |= {i for i in range(len(Rs) - 1) if np.sign(D[i, s, 0]) != np.sign(D[i + 1, s, 0])}
    seen = []
    for i in sorted(cand):
        Rr = None
        for k in range(8):
            f = lambda r: (lambda c: c[0][s][k + 1] - c[1][s][k + 1])(coeffs(r))
            lo, hi = Rs[max(i - 1, 0)], Rs[min(i + 2, len(Rs) - 1)]
            try:
                Rr = brentq(f, lo, hi, xtol=1e-12, rtol=1e-14)
                break
            except ValueError:
                continue
        if Rr is None or any(abs(Rr - x) < 1e-7 for x in seen):
            continue
        c1, c2 = (c[s] for c in coeffs(Rr))
        agree = float(np.max(np.abs(c1[1:20] - c2[1:20])) / max(1.0, np.max(np.abs(c1[1:20]))))
        c = c1
        hecke = max(abs(c[1] * c[2] - c[5]), abs(c[1] * c[4] - c[9]), abs(c[3] - (c[1] ** 2 - 1)), abs(c[2] * c[4] - c[14]), abs(c[8] - (c[2] ** 2 - 1)))
        if agree > 1e-6 or hecke > 1e-4:
            continue
        seen.append(Rr)
        new = all(abs(abs(c[p - 1]) * math.sqrt(p) - 1) < 1e-4 for p in ps)
        found.append({"sector": f"odd e{ps[0]}={a:+d} e{ps[1]}={b:+d}", "R": Rr, "lambda": 0.25 + Rr * Rr, "ratio_to_quarter": (0.25 + Rr * Rr) / 0.25,
                      "new": bool(new), "agree": agree, "hecke": float(hecke)})
        print(f"  odd e{ps[0]}={a:+d} e{ps[1]}={b:+d}: R = {Rr:.10f}  lambda/(1/4) = {(0.25 + Rr * Rr) / 0.25:.5f}  {'new' if new else 'OLD'}  agree {agree:.1e}  hecke {hecke:.1e}", flush=True)
found.sort(key=lambda e: e["R"])
json.dump({"N": N, "Y0": Y0, "M": M, "Qs": Qs, "R_range": [R0, R1], "step": dR, "found": found}, open(out, "w"), indent=1)
print(f"level {N}: lowest odd eps_N = +1 forms {[round(e['R'], 6) for e in found[:3]]}; {time.time() - T0:.0f} s")
