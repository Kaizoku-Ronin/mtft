#!/usr/bin/env python3
"""MP15, the checks added during the run (same solver as mp15_maass.py):

  python3 mp15_extra.py dips             the three eigenvalues at level 143 that the sign-change filter missed, found as dips of the
                                         horocycle disagreement (near R = 0.8993, 1.4143, 1.4640), refined and checked; writes roots_dips.json
  python3 mp15_extra.py lmfdb 21 1 0.85 1.0
                                         a level and parity scanned for sign changes of c_2 between the horocycles, each refined and
                                         checked: at level 21, odd, it finds LMFDB's form 21.1 (R = 0.92400721, Fricke -1)
  python3 mp15_extra.py sectors 133      all eight sectors at a level, R from 0.05 to 1.5 in steps of 0.025: where c_2 changes sign
"""
import importlib.util
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

spec = importlib.util.spec_from_file_location("m", os.path.join(os.path.dirname(os.path.abspath(__file__)), "mp15_maass.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def horocycles(N, y_probe):
    probe = [m.reduce_point(float(x), y_probe, N)[1] for x in np.linspace(0.0001, 0.4999, 4000)]
    Y0 = min(probe)
    while True:
        M = int(math.ceil(40.0 / (2 * math.pi * Y0)))
        Qs = int(math.ceil(1.08 * M / 0.88)) + 12
        try:
            return Y0, (m.Horocycle(N, 0.96 * Y0, Qs), m.Horocycle(N, 0.88 * Y0, Qs)), M
        except AssertionError:
            Y0 *= 0.9


def coeffs(hcs, M, R, par, chis_eps):
    l = np.arange(1, M + 1)
    out = []
    for hc in hcs:
        Kc = (m.Kir(R, 2 * np.pi * np.outer(hc.ys, l)), m.Kir(R, 2 * np.pi * l * hc.Y))
        out.append([c for c, _ in m.solve_coeffs(R, hc, M, par, [hc.chi(e) for e in chis_eps], Kc)])
    return out


def checks_at(c1, c2, ps):
    """Agreement of c_2 .. c_20 between the horocycles, and Hecke relations: c_mn = c_m c_n for coprime m, n; c_{p^2} = c_p^2 - 1
    for p not dividing the level, c_p^2 for p dividing it."""
    c = c1
    agree = float(np.max(np.abs(c1[1:20] - c2[1:20])) / max(1.0, np.max(np.abs(c1[1:20]))))
    sq = lambda p: c[p - 1] ** 2 - (0 if p in ps else 1)
    hecke = max(abs(c[1] * c[2] - c[5]), abs(c[1] * c[4] - c[9]), abs(c[3] - sq(2)), abs(c[2] * c[4] - c[14]), abs(c[8] - sq(3)))
    return agree, float(hecke), {p: float(c[p - 1] * math.sqrt(p)) for p in ps}


def dips():
    N, Y0, hcs, M = m.setup143()
    res = []
    for s, R0 in ((0, 0.89925523), (2, 1.41428254), (7, 1.46398671)):
        par, e11, e13 = m.SECTORS[s]
        R = None
        for k in (2, 3, 4, 5, 1):
            f = lambda r: (lambda c: c[0][0][k] - c[1][0][k])(coeffs(hcs, M, r, par, [(e11, e13)]))
            try:
                R = brentq(f, R0 - 2e-4, R0 + 2e-4, xtol=1e-12, rtol=1e-14)
                break
            except ValueError:
                continue
        c = coeffs(hcs, M, R, par, [(e11, e13)])
        agree, hecke, cp = checks_at(c[0][0], c[1][0], (11, 13))
        new = all(abs(abs(v) - 1) < 1e-4 for v in cp.values())
        print(f"{m.sector_name(m.SECTORS[s])}: R = {R:.10f}  lambda = {0.25 + R * R:.8f}  {'new' if new else 'OLD'}  agree {agree:.1e}  hecke {hecke:.1e}  {cp}")
        res.append({"sector": m.sector_name(m.SECTORS[s]), "R": R, "lambda": 0.25 + R * R, "new": bool(new), "agree": agree, "hecke": hecke,
                    "c11_sqrt11": cp[11], "c13_sqrt13": cp[13], "found_by": "dip of the horocycle disagreement (missed by the sign-change filter)"})
    json.dump(res, open("roots_dips.json", "w"), indent=1)


def lmfdb(N, par, lo, hi):
    Y0, hcs, M = horocycles(N, 0.01)
    ps = m.primes_of(N)
    print(f"level {N} {ps}: lowest point met {Y0:.5f}, M = {M}")
    epss = [(a, b) for a in (1, -1) for b in (1, -1)]
    grid = np.arange(lo, hi + 1e-9, 0.005)
    for j, eps in enumerate(epss):
        d = [(lambda c: c[0][0][1] - c[1][0][1])(coeffs(hcs, M, R, par, [eps])) for R in grid]
        for i in range(len(grid) - 1):
            if np.sign(d[i]) == np.sign(d[i + 1]):
                continue
            try:
                R = brentq(lambda r: (lambda c: c[0][0][1] - c[1][0][1])(coeffs(hcs, M, r, par, [eps])), grid[i], grid[i + 1], xtol=1e-12)
            except ValueError:
                continue
            c = coeffs(hcs, M, R, par, [eps])
            agree, hecke, cp = checks_at(c[0][0], c[1][0], ps)
            print(f"  eps{ps} = {eps}: R = {R:.8f}  agree {agree:.1e}  hecke {hecke:.1e}  {'(a pole)' if agree > 1e-6 else ''}  {cp}")


def sectors(N):
    Y0, hcs, M = horocycles(N, 0.002)
    grid = np.arange(0.05, 1.501, 0.025)
    epss = [(a, b) for a in (1, -1) for b in (1, -1)]
    D = np.zeros((len(grid), 8))
    for i, R in enumerate(grid):
        for par in (0, 1):
            c = coeffs(hcs, M, R, par, epss)
            for j in range(4):
                D[i, 4 * par + j] = c[0][j][1] - c[1][j][1]
    ps = m.primes_of(N)
    for s in range(8):
        par, a, b = s // 4, epss[s % 4][0], epss[s % 4][1]
        sc = [round(float(grid[i]), 3) for i in range(len(grid) - 1) if np.sign(D[i, s]) != np.sign(D[i + 1, s])]
        print(f"level {N} {'even' if par == 0 else 'odd'} eps{ps} = ({a:+d}, {b:+d}): c_2 changes sign at {sc}")


if __name__ == "__main__":
    if sys.argv[1] == "dips":
        dips()
    elif sys.argv[1] == "lmfdb":
        lmfdb(int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
    elif sys.argv[1] == "sectors":
        sectors(int(sys.argv[2]))
