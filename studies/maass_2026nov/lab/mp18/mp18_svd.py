#!/usr/bin/env python3
"""MP18, a pole-free detector for Route H.  Hejhal's normalized system (c_1 = 1, first equation dropped) has poles wherever the
reduced matrix is singular, and the twisted sectors are crowded with them (find2 and find3 rejected about 160 candidates).
The full homogeneous system V c = 0 on one horocycle has no normalization and so no poles: its smallest singular value,
relative to its largest, dips toward 0 only where a form exists.

  python3 -I mp18_svd.py scan f R0 R1 dR out.npz      sigma_min / sigma_max for all 8 subsectors on the horocycle f Y0
                                                      (a 300-coefficient system: a detector, not a solver)
  python3 -I mp18_svd.py find out.json scan1.npz ...  dips -> refined by Brent on the two-horocycle difference (pair A, then
                                                      B, then C) -> accepted only under P3's rules
  python3 -I mp18_svd.py dips found.json out.json [0|1]   every unconfirmed dip of find, tested with the full system on the
                                                      horocycles 0.96 Y0 and 0.88 Y0 (the optional 0 or 1 takes every other dip)
  python3 -I mp18_svd.py fine which R0 R1 dR s1,s2 out.json   the full system on one horocycle (which = 0 or 1) on a fine grid,
                                                      for subsectors numbered 0 to 7 (even, then odd; ++, +-, -+, --)"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("tw", os.path.join(HERE, "mp18_twisted.py"))
tw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tw)
T0 = time.time()


def vmat(R, hc, fr, parity, chi, Kc):
    Kst, Kdiag = Kc
    M = fr.size
    trig = np.cos if parity == 0 else np.sin
    P = np.sqrt(hc.ys)[:, None] * Kst * trig(2 * np.pi * np.outer(hc.xs, fr))
    Tn = trig(2 * np.pi * np.outer(fr, hc.xm)) * (2.0 / hc.Qs)
    V = (Tn * chi[None, :]) @ P
    V[np.arange(M), np.arange(M)] -= np.sqrt(hc.Y) * Kdiag
    sc = np.abs(V).max(axis=0)
    sc[sc == 0] = 1.0
    return V / sc[None, :]


def smin(R, hc, fr, s, Kc=None):
    par, eps = tw.SUB[s]
    if Kc is None:
        Kc = tw.kcache(R, hc, fr)
    sv = np.linalg.svd(vmat(R, hc, fr, par, hc.chi(eps), Kc), compute_uv=False)
    return sv[-1] / sv[0]


MD = 300                                                # detector size: K_iR(2 pi 300 Y) ~ e^-17, enough to see a dip


def det_setup(f):
    Y0 = tw.setup(0.96, 0.88)[0]
    Qs = int(math.ceil(1.08 * MD / 0.88)) + 12
    return Y0, tw.Horo(f * Y0, Qs), np.arange(1, MD + 1) - 0.5


def scan(f, R0, R1, dR, out):
    Y0, hc, fr = det_setup(f)
    M = MD
    Rs = np.arange(R0, R1 + 1e-12, dR)
    S = np.zeros((len(Rs), 8))
    for i, R in enumerate(Rs):
        Kc = tw.kcache(R, hc, fr)
        for s in range(8):
            S[i, s] = smin(R, hc, fr, s, Kc)
        if i % 40 == 0:
            print(f"  R = {R:.4f} ({i + 1}/{len(Rs)}, {time.time() - T0:.0f} s)", flush=True)
    np.savez(out, Rs=Rs, S=S, f=f, Y0=Y0, M=M)
    print(f"svd scan {R0}..{R1} on {f} Y0: {time.time() - T0:.0f} s")


def find(out, files):
    Z = [np.load(f) for f in files]
    Rs = np.concatenate([z["Rs"] for z in Z])
    S = np.concatenate([z["S"] for z in Z])
    Rs, first = np.unique(np.round(Rs, 9), return_index=True)        # the two halves overlap at their seam
    S = S[first]
    f = float(Z[0]["f"])
    Y0, M, hcsA = tw.setup(0.96, 0.88)
    Qs = hcsA[0].Qs
    pairs = {"A": hcsA, "B": (tw.Horo(0.93 * Y0, Qs), tw.Horo(0.85 * Y0, Qs)), "C": (tw.Horo(0.90 * Y0, Qs), tw.Horo(0.82 * Y0, Qs))}
    _, hc, frd = det_setup(f)
    fr = np.arange(1, M + 1) - 0.5
    found, rejected = [], []
    L = np.log(S)
    for s in range(8):
        par, eps = tw.SUB[s]
        dips = [i for i in range(1, len(Rs) - 1) if L[i, s] < L[i - 1, s] and L[i, s] <= L[i + 1, s]]
        for i in dips:
            res = minimize_scalar(lambda r: math.log(smin(r, hc, frd, s)), bounds=(Rs[i - 1], Rs[i + 1]), method="bounded",
                                  options={"xatol": 1e-10})
            r1 = float(res.x)
            if res.fun > -13:                               # a shallow minimum of the background, not a form
                continue
            got = None
            for name, hcs in pairs.items():
                F = lambda r: (lambda cs: cs[0][1] - cs[1][1])(tw.coeffs(r, hcs, fr, par, eps))
                try:
                    Rr = brentq(F, r1 - 2e-4, r1 + 2e-4, xtol=1e-12, rtol=1e-14)
                except ValueError:
                    continue
                c1, c2 = tw.coeffs(Rr, hcs, fr, par, eps)
                agree, hecke, b11, b13 = tw.tests(c1, c2)
                if agree <= 1e-6 and hecke <= 1e-4:
                    got = (name, Rr, agree, hecke, b11, b13, c1)
                    break
            if got is None:
                rejected.append({"subsector": tw.sub_name(s), "dip": r1, "log smin": float(res.fun)})
                print(f"  {tw.sub_name(s)}: dip at {r1:.6f} (log smin {res.fun:.1f}) not confirmed", flush=True)
                continue
            name, Rr, agree, hecke, b11, b13, c1 = got
            new = abs(abs(b11) - 1) < 1e-4 and abs(abs(b13) - 1) < 1e-4
            found.append({"subsector": tw.sub_name(s), "parity": "even" if par == 0 else "odd", "eps11": eps[0], "eps13": eps[1], "R": Rr,
                          "lambda": 0.25 + Rr * Rr, "new": bool(new), "agree": agree, "hecke": hecke, "b3": float(c1[1]), "b5": float(c1[2]),
                          "b11_sqrt11": b11, "b13_sqrt13": b13, "confirmed on pair": name, "log smin at dip": float(res.fun)})
            print(f"  {tw.sub_name(s)}: R = {Rr:.10f}  lambda/(1/4) = {(0.25 + Rr * Rr) * 4:.6f}  {'new' if new else 'OLD'}  pair {name}  "
                  f"agree {agree:.1e}  hecke {hecke:.1e}  (dip log smin {res.fun:.1f})", flush=True)
    json.dump({"found": sorted(found, key=lambda e: e["R"]), "unconfirmed dips": rejected, "seconds": time.time() - T0}, open(out, "w"), indent=1)
    print(f"svd find: {len(found)} forms, {len(rejected)} unconfirmed dips, {time.time() - T0:.0f} s")


def dips(found_file, out, part=None):
    """Every dip that find left unconfirmed, tested with the full system (all M coefficients) on the two horocycles of pair A:
    the minimum of log(sigma_min / sigma_max) within 0.003 of the dip.  A form dips to about -26 at the same R on both
    horocycles; a dip that is shallow, at the window's edge or at a different R on the second horocycle is not a form."""
    d = json.load(open(found_file))["unconfirmed dips"]
    Y0, M, hcs = tw.setup(0.96, 0.88)
    fr = np.arange(1, M + 1) - 0.5
    names = [tw.sub_name(s) for s in range(8)]
    rows = []
    for i, x in enumerate(d):
        if part is not None and i % 2 != part:
            continue
        s = names.index(x["subsector"])
        lo, hi = x["dip"] - 0.003, x["dip"] + 0.003
        row = {"index": i, "subsector": x["subsector"], "dip (300-coefficient detector, 0.96 Y0)": x["dip"], "its log smin": x["log smin"]}
        for name, hc in (("0.96 Y0", hcs[0]), ("0.88 Y0", hcs[1])):
            res = minimize_scalar(lambda r: math.log(smin(r, hc, fr, s)), bounds=(lo, hi), method="bounded", options={"xatol": 1e-8})
            row[name] = {"R at minimum": float(res.x), "log smin": float(res.fun), "at the window's edge": bool(min(res.x - lo, hi - res.x) < 2e-5)}
        rows.append(row)
        a, b = row["0.96 Y0"], row["0.88 Y0"]
        print(f"  {x['subsector']}  dip {x['dip']:.6f}: full system on 0.96 Y0 -> {a['R at minimum']:.6f} ({a['log smin']:.1f})   "
              f"on 0.88 Y0 -> {b['R at minimum']:.6f} ({b['log smin']:.1f}){' [edge]' if b["at the window's edge"] else ''}", flush=True)
    json.dump({"full-system coefficients": M, "rows": rows, "seconds": time.time() - T0}, open(out, "w"), indent=1)
    print(f"dips: {len(rows)} tested, {time.time() - T0:.0f} s")


def fine(which, R0, R1, dR, subs, out):
    """log(sigma_min / sigma_max) of the full system on one horocycle of pair A (which = 0: 0.96 Y0, 1: 0.88 Y0), on a fine
    grid, for the listed subsectors: the test that separates a form (the same dip on both horocycles) from a dip of one."""
    Y0, M, hcs = tw.setup(0.96, 0.88)
    fr = np.arange(1, M + 1) - 0.5
    hc = hcs[which]
    Rs = np.arange(R0, R1 + 1e-9, dR)
    res = {tw.sub_name(s): [float(math.log(smin(float(R), hc, fr, s))) for R in Rs] for s in subs}
    json.dump({"horocycle": ["0.96 Y0", "0.88 Y0"][which], "full-system coefficients": M, "Rs": [float(r) for r in Rs], "log smin": res,
               "seconds": time.time() - T0}, open(out, "w"))
    for k, v in res.items():
        i = int(np.argmin(v))
        print(f"  {['0.96 Y0', '0.88 Y0'][which]}  {k}: minimum {v[i]:.2f} at R = {Rs[i]:.4f}; median {np.median(v):.2f}; maximum {max(v):.2f}", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "fine":
        fine(int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), [int(v) for v in sys.argv[6].split(",")], sys.argv[7])
    elif sys.argv[1] == "scan":
        scan(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), sys.argv[6])
    elif sys.argv[1] == "dips":
        dips(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else None)
    else:
        find(sys.argv[2], sys.argv[3:])
