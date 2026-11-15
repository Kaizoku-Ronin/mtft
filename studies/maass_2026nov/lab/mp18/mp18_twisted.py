#!/usr/bin/env python3
"""MP18, the eta-twist: Maass forms on Gamma_0(143) with the character chi = nu_eta^12, the sign of the permutation that gamma mod 2
makes of P^1(F_2).  All four cusp widths (1, 11, 13, 143) are odd, so chi = -1 around every cusp: no Eisenstein series, no continuum.
The twisted space splits by parity and by the characters X_eps = eps * chi~ of Gamma_0*(143), where chi~(M) is the sign of the
integer matrix M (odd determinant) mod 2 and eps = (eps_11, eps_13) is the Atkin-Lehner part.

Predictions P1 to P4 and F1 are in MP18_Eta_Twist_Plan.md, saved at 08:52:21 New York time before this file existed.

Route H (Hejhal's method): MP15's solver with two changes.  The frequencies are n - 1/2, since f(z + 1) = -f(z), so
    f = sum_{n >= 1} c_n sqrt(y) K_iR(2 pi (n - 1/2) y) trig(2 pi (n - 1/2) x),  c_1 = 1,
and the half-range points x_m = (m - 1/2)/(2 Qs) stay orthogonal (DCT-IV and DST-IV).  The reduction into the Ford domain of
Gamma_0*(143) carries, besides the Atkin-Lehner parities, the sign mod 2 of every step's integer matrix (T^k gives (-1)^k).
Route T (the twisted trace formula, both parities and all eps together):
    sum_j h(r_j) = (vol/4 pi) int h(r) r tanh(pi r) dr + sum_t w_t^chi g(l_t) - 4 g(0) log 2,
with 2 cosh(l_t/2) = t and w_t^chi = 2 W^chi(t^2 - 4): MP16's W_of with each term of index u signed by chi, +1 for odd t and
(-1)^u for even t (V2 checks the rule).  No elliptic terms (Gamma_0(143) is torsion-free), no continuum, no scattering terms.

Usage: python3 -I mp18_twisted.py checks                         V1, V2, V3 -> mp18_checks.json
       python3 -I mp18_twisted.py scan f1 f2 R0 R1 dR out.npz    horocycles at f1 Y0 and f2 Y0
       python3 -I mp18_twisted.py find out.json lo hi scan.npz   candidates refined and checked (P3's rules)
       python3 -I mp18_twisted.py find2 out.json lo hi scan.npz  the same, faster (interpolation, then secant steps)
       python3 -I mp18_twisted.py find3 out.json lo hi scan.npz known1.json ...   find2 with every candidate, known forms skipped
       python3 -I mp18_twisted.py weights                        Route T's signed weights to length 23.5 -> weights_tw.npz
       python3 -I mp18_twisted.py trace out.json list1.json [list2.json ...]
Needs numpy, scipy, mpmath, cypari2, and mp15_maass.py, mp16_trace.py in ../mp15, ../mp16 (MP17's weights_143_1.npz in ../mp17)."""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m = _load("mp15_maass", os.path.join(HERE, "..", "mp15", "mp15_maass.py"))
T0 = time.time()
N = 143
VOL = 56 * math.pi                                     # index 168 times pi/3
EPS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
SUB = [(par, e) for par in (0, 1) for e in EPS]        # 8 subsectors: (parity, (eps11, eps13))


def sub_name(s):
    par, (e11, e13) = SUB[s]
    return f"{'even' if par == 0 else 'odd'} e11={e11:+d} e13={e13:+d}"


# ------------------------------------------------------------------------------------------- chi and the reduction
def sign_mod2(a, b, c, d):
    """chi~ of an integer matrix of odd determinant: the sign of its permutation of P^1(F_2).  In GL(2, F_2) = S_3 the identity
    and the two elements of trace 1 (3-cycles) give +1; the three other elements (trace 0, order 2) give -1."""
    a, b, c, d = a & 1, b & 1, c & 1, d & 1
    if (a, b, c, d) == (1, 0, 0, 1):
        return 1
    return 1 if (a + d) & 1 else -1


def reduce_tw(x, y):
    """MP15's reduction (mp15_maass.reduce_point) with the mod-2 sign carried along.  Returns x*, y*, the parities of the
    Atkin-Lehner steps at p = 11, 13, and s: a form with character X_eps satisfies f(z) = eps^parity * s * f(z*)."""
    Qs = (1, 11, 13, 143)
    par = [0, 0]
    s = 1
    for _ in range(10000):
        k = math.floor(x + 0.5)
        x -= k
        if k & 1:
            s = -s
        best, arg = 1.0 - 1e-13, None
        for Q in Qs:
            Nq = N // Q
            cmax = int(math.sqrt(Q) / (N * y)) + 1
            for c in range(1, cmax + 1):
                fl = math.floor(-N * c * x / Q)
                for d in range(fl - 1, fl + 3):
                    if math.gcd(Q * d, Nq * c) != 1:
                        continue
                    v = ((N * c * x + Q * d) ** 2 + (N * c * y) ** 2) / Q
                    if v < best:
                        best, arg = v, (Q, c, d)
        if arg is None:
            return x, y, tuple(par), s
        Q, c, d = arg
        _, u, w = m.egcd(Q * d, (N // Q) * c)
        a, b = u, -w                                      # det [[Q a, b], [N c, Q d]] = Q
        den = (N * c * x + Q * d) ** 2 + (N * c * y) ** 2
        xn = ((Q * a * x + b) * (N * c * x + Q * d) + Q * a * N * c * y * y) / den
        y = Q * y / den
        x = xn
        s *= sign_mod2(Q * a, b, N * c, Q * d)
        if Q % 11 == 0:
            par[0] ^= 1
        if Q % 13 == 0:
            par[1] ^= 1
    raise RuntimeError("no convergence")


class Horo:
    def __init__(self, Y, Qs):
        self.Y, self.Qs = Y, Qs
        self.xm = (np.arange(1, Qs + 1) - 0.5) / (2 * Qs)
        red = [reduce_tw(float(x), Y) for x in self.xm]
        self.xs = np.array([r[0] for r in red])
        self.ys = np.array([r[1] for r in red])
        self.par = np.array([r[2] for r in red], dtype=int)
        self.sg = np.array([r[3] for r in red], dtype=float)
        assert self.ys.min() > Y * 1.0001, (self.ys.min(), Y)

    def chi(self, eps, twist=True):
        e = np.array(eps, dtype=float)
        al = np.prod(np.where(self.par == 1, e[None, :], 1.0), axis=1)
        return al * self.sg if twist else al


def setup(f1=0.96, f2=0.88):
    probe = [m.reduce_point(float(x), 0.002, N)[1] for x in np.linspace(0.0001, 0.4999, 4000)]
    Y0 = min(probe)
    M = int(math.ceil(40.0 / (2 * math.pi * Y0)))
    Qs = int(math.ceil(1.08 * M / 0.88)) + 12
    return Y0, M, (Horo(f1 * Y0, Qs), Horo(f2 * Y0, Qs))


def kcache(R, hc, fr):
    return m.Kir(R, 2 * np.pi * np.outer(hc.ys, fr)), m.Kir(R, 2 * np.pi * fr * hc.Y)


def solve(R, hc, fr, parity, chis, Kc):
    """MP15's solve_coeffs with the frequencies fr (n - 1/2 for the twist, n without it)."""
    Kst, Kdiag = Kc
    M = fr.size
    trig = np.cos if parity == 0 else np.sin
    P = np.sqrt(hc.ys)[:, None] * Kst * trig(2 * np.pi * np.outer(hc.xs, fr))
    Tn = trig(2 * np.pi * np.outer(fr, hc.xm)) * (2.0 / hc.Qs)
    out = []
    for chi in chis:
        V = (Tn * chi[None, :]) @ P
        V[np.arange(M), np.arange(M)] -= np.sqrt(hc.Y) * Kdiag
        A = V[1:, 1:]
        sc = np.abs(A).max(axis=0)
        sc[sc == 0] = 1.0
        cr = np.linalg.solve(A / sc[None, :], -V[1:, 0]) / sc
        out.append(np.concatenate([[1.0], cr]))
    return out


def coeffs(R, hcs, fr, parity, eps, twist=True):
    return [solve(R, hc, fr, parity, [hc.chi(eps, twist)], kcache(R, hc, fr))[0] for hc in hcs]


def tests(c1, c2):
    """P3's rules.  b_m is the coefficient of frequency m/2, m odd: b_m = c[(m - 1)/2]."""
    agree = float(np.max(np.abs(c1[1:20] - c2[1:20])) / max(1.0, np.max(np.abs(c1[1:20]))))
    c = c1
    hecke = float(max(abs(c[1] * c[2] - c[7]),             # b3 b5 = b15
                      abs(c[1] * c[3] - c[10]),            # b3 b7 = b21
                      abs(c[4] - (c[1] ** 2 - 1)),         # b9 = b3^2 - 1
                      abs(c[12] - (c[2] ** 2 - 1)),        # b25 = b5^2 - 1
                      abs(c[2] * c[3] - c[17])))           # b5 b7 = b35
    return agree, hecke, float(c[5] * math.sqrt(11)), float(c[6] * math.sqrt(13))


# ------------------------------------------------------------------------------------------- checks V1, V2, V3
def checks():
    rng = np.random.default_rng(18)
    out = {}
    # V1: the carried sign is a character
    n_ok = n_bad = n_skip = 0
    while n_ok + n_bad < 200:
        k = int(rng.choice([-2, -1, 1, 2]))
        c = 143 * k
        d = int(rng.integers(-40, 41))
        if math.gcd(c, d) != 1:
            continue
        _, u, v = m.egcd(d, c)                                  # u d + v c = 1
        a, b = u, -v                                            # a d - b c = 1
        z = complex(-d / c + rng.uniform(-0.02, 0.02), rng.uniform(0.01, 0.04))
        gz = (a * z + b) / (c * z + d)
        r1 = reduce_tw(z.real, z.imag)
        r2 = reduce_tw(gz.real, gz.imag)
        if abs(r1[0] - r2[0]) + abs(r1[1] - r2[1]) > 1e-9 * max(1.0, r1[1]):
            n_skip += 1                                         # met the domain at two boundary points
            continue
        if r1[2] == r2[2] and r2[3] == r1[3] * sign_mod2(a, b, c, d):
            n_ok += 1
        else:
            n_bad += 1
    out["V1"] = {"pairs": n_ok + n_bad, "consistent": n_ok, "inconsistent": n_bad, "skipped (boundary)": n_skip, "pass": n_bad == 0}
    print(f"V1: {n_ok} of {n_ok + n_bad} pairs consistent ({n_skip} skipped at the boundary)", flush=True)
    # V2: the sign of a hyperbolic class from (t, u)
    cnt = {"t odd": 0, "t even, u odd": 0, "t even, u even": 0}
    bad = 0
    for i in range(3000):
        k = int(rng.integers(1, 12)) * int(rng.choice([-1, 1]))
        c = 143 * k
        d = int(rng.integers(-400, 401))
        if math.gcd(c, d) != 1:
            continue
        _, u_, v_ = m.egcd(d, c)
        a, b = u_, -v_
        mats = [(a, b, c, d)]
        if (a + d) % 2 == 0:                                    # squares of even-trace elements reach u even
            mats.append((a * a + b * c, a * b + b * d, c * a + d * c, c * b + d * d))
        for (A, B, C, Dd) in mats:
            t = abs(A + Dd)
            if t <= 2:
                continue
            uu = math.gcd(math.gcd(C, Dd - A), B)
            rule = 1 if t % 2 else (-1 if uu % 2 else 1)
            key = "t odd" if t % 2 else ("t even, u even" if uu % 2 == 0 else "t even, u odd")
            cnt[key] += 1
            if rule != sign_mod2(A, B, C, Dd):
                bad += 1
    out["V2"] = {"counts": cnt, "disagreements": bad, "pass": bad == 0 and min(cnt.values()) > 20}
    print(f"V2: {cnt}, disagreements {bad}", flush=True)
    # V3: with the twist switched off, MP15's lowest odd form (eps = (-1, -1)) is found again
    Y0, M, hcs = setup(0.96, 0.88)
    fr = np.arange(1, M + 1, dtype=float)
    g = lambda R: (lambda cs: cs[0][1] - cs[1][1])(coeffs(R, hcs, fr, 1, (-1, -1), twist=False))
    R = brentq(g, 0.3730, 0.3746, xtol=1e-13, rtol=1e-14)
    out["V3"] = {"R": R, "MP15": 0.3738518157494396, "difference": R - 0.3738518157494396, "pass": abs(R - 0.3738518157494396) < 1e-8,
                 "Y0": Y0, "M": M, "Qs": hcs[0].Qs}
    print(f"V3: R = {R:.12f} against MP15's 0.373851815749, difference {R - 0.3738518157494396:.1e}", flush=True)
    json.dump(out, open(os.path.join(HERE, "mp18_checks.json"), "w"), indent=1)
    print(f"checks done, {time.time() - T0:.0f} s")


# ------------------------------------------------------------------------------------------- Route H: scan and find
def scan(f1, f2, R0, R1, dR, out):
    Y0, M, hcs = setup(f1, f2)
    fr = np.arange(1, M + 1) - 0.5
    chis = [[hc.chi(e) for e in EPS] for hc in hcs]
    Rs = np.arange(R0, R1 + 1e-12, dR)
    D = np.zeros((len(Rs), 8, 8))
    t = time.time()
    for i, R in enumerate(Rs):
        sols = []
        for h, hc in enumerate(hcs):
            Kc = kcache(R, hc, fr)
            row = []
            for parity in (0, 1):
                row += solve(R, hc, fr, parity, chis[h], Kc)
            sols.append(row)
        for s in range(8):
            D[i, s, :] = (sols[0][s] - sols[1][s])[1:9]
        if i % 40 == 0:
            print(f"  R = {R:.4f}  ({i + 1}/{len(Rs)}, {time.time() - t:.0f} s)", flush=True)
    np.savez(out, Rs=Rs, D=D, Y0=Y0, M=M, Qs=hcs[0].Qs, f=[f1, f2])
    print(f"scan {R0}..{R1} step {dR}, horocycles {f1} and {f2} Y0: Y0 = {Y0:.6f}, M = {M}, Qs = {hcs[0].Qs}, {time.time() - T0:.0f} s")


def find(out, lo, hi, scanfile):
    """mp16_find's search on every subsector of one scan: candidates are sign changes of the first three coefficient differences
    and local minima of their norm; each is refined by Brent's method and kept only under P3's rules (agreement <= 1e-6,
    Hecke <= 1e-4).  Rejected brackets are poles of the linear system."""
    z = np.load(scanfile)
    f1, f2 = [float(v) for v in z["f"]]
    Y0, M, hcs = setup(f1, f2)
    assert abs(Y0 - float(z["Y0"])) < 1e-12 and M == int(z["M"])
    fr = np.arange(1, M + 1) - 0.5
    Rs, D = z["Rs"], z["D"]
    keep = (Rs >= lo - 1e-12) & (Rs <= hi + 1e-12)
    Rs, D = Rs[keep], D[keep]
    found, rejected = [], []
    for s in range(8):
        par, eps = SUB[s]
        Dg = D[:, s, :]
        nrm = np.linalg.norm(Dg, axis=1)
        cand = set()
        for i in range(len(Rs) - 1):
            if any(np.sign(Dg[i, k]) != np.sign(Dg[i + 1, k]) for k in range(3)):
                cand.add(i)
        for i in range(1, len(Rs) - 1):
            if nrm[i] < nrm[i - 1] and nrm[i] < nrm[i + 1]:
                cand.add(i)
        mine = []
        for i in sorted(cand):
            a, b = Rs[max(i - 1, 0)], Rs[min(i + 2, len(Rs) - 1)]
            if any(a - 1e-9 <= e["R"] <= b + 1e-9 for e in mine):
                continue
            Rr = None
            for k in range(5):
                fk = lambda r, k=k: (lambda cs: cs[0][k + 1] - cs[1][k + 1])(coeffs(r, hcs, fr, par, eps))
                try:
                    Rr = brentq(fk, a, b, xtol=1e-12, rtol=1e-14, maxiter=60)
                    break
                except ValueError:
                    continue
            if Rr is None:
                continue
            c1, c2 = coeffs(Rr, hcs, fr, par, eps)
            agree, hecke, b11, b13 = tests(c1, c2)
            if agree > 1e-6 or hecke > 1e-4:
                rejected.append({"subsector": sub_name(s), "R": Rr, "agree": agree, "hecke": hecke})
                continue
            new = abs(abs(b11) - 1) < 1e-4 and abs(abs(b13) - 1) < 1e-4
            e = {"subsector": sub_name(s), "parity": "even" if par == 0 else "odd", "eps11": eps[0], "eps13": eps[1], "R": Rr,
                 "lambda": 0.25 + Rr * Rr, "new": bool(new), "agree": agree, "hecke": hecke, "b3": float(c1[1]), "b5": float(c1[2]),
                 "b11_sqrt11": b11, "b13_sqrt13": b13}
            mine.append(e)
            print(f"  {sub_name(s)}: R = {Rr:.10f}  lambda/(1/4) = {(0.25 + Rr * Rr) * 4:.6f}  {'new' if new else 'OLD'}  agree {agree:.1e}  "
                  f"hecke {hecke:.1e}  b11 sqrt11 {b11:+.6f}  b13 sqrt13 {b13:+.6f}", flush=True)
        found += mine
    json.dump({"scan": os.path.basename(scanfile), "horocycle fractions": [f1, f2], "window": [lo, hi], "Y0": Y0, "M": M,
               "found": sorted(found, key=lambda e: e["R"]), "rejected (poles)": rejected, "seconds": time.time() - T0},
              open(out, "w"), indent=1)
    print(f"find: {len(found)} forms, {len(rejected)} rejected, {time.time() - T0:.0f} s")


def find2(out, lo, hi, scanfile, exhaustive=False, known=()):
    """A faster find, used for the record (find, mp16_find's search, spent most of its time on false candidates).
    Candidates: sign changes of the c_2 or c_3 difference that pass MP15's pole test (|difference| must shrink into the
    bracket), and deep dips of the difference norm (below a quarter of its value two steps either side), which is where
    MP16's hidden forms were.  Each starts at the scan's linear interpolation and is refined by secant steps on the c_2
    difference (at most 10 evaluations, and it must stay within 0.01); accepted only under P3's rules."""
    z = np.load(scanfile)
    f1, f2 = [float(v) for v in z["f"]]
    Y0, M, hcs = setup(f1, f2)
    assert abs(Y0 - float(z["Y0"])) < 1e-12 and M == int(z["M"])
    fr = np.arange(1, M + 1) - 0.5
    keep = (z["Rs"] >= lo - 0.0101) & (z["Rs"] <= hi + 0.0101)
    Rs, D = z["Rs"][keep], z["D"][keep]
    found, rejected, nev = [], [], 0
    for s in range(8):
        par, eps = SUB[s]
        d = D[:, s, :]
        nrm = np.linalg.norm(d, axis=1)
        cands = []
        for k in ((0, 1, 2) if exhaustive else (0, 1)):
            for i in range(len(Rs) - 1):
                if np.sign(d[i, k]) == np.sign(d[i + 1, k]):
                    continue
                if not exhaustive and 2 <= i < len(Rs) - 3 and max(abs(d[i, k]), abs(d[i + 1, k])) >= min(abs(d[i - 2, k]), abs(d[i + 3, k])):
                    continue
                cands.append(Rs[i] + d[i, k] / (d[i, k] - d[i + 1, k]) * (Rs[i + 1] - Rs[i]))
        for i in range(2, len(Rs) - 2):
            if nrm[i] <= nrm[i - 1] and nrm[i] <= nrm[i + 1] and (exhaustive or nrm[i] < 0.25 * min(nrm[i - 2], nrm[i + 2])):
                cands.append(float(Rs[i]))
        cands = sorted(c for c in cands if lo - 0.005 <= c <= hi + 0.005)
        merged = []
        for c in cands:
            if not merged or c - merged[-1] > (0.002 if exhaustive else 0.004):
                merged.append(c)
        kn = [R for (nm_, R) in known if nm_ == sub_name(s)]
        merged = [c for c in merged if all(abs(c - R) > 0.003 for R in kn)]
        mine = []
        for R0 in merged:
            F = lambda r: (lambda cs: cs[0][1] - cs[1][1])(coeffs(r, hcs, fr, par, eps))
            a, b = R0, R0 + 2e-4
            fa, fb = F(a), F(b)
            nev += 2
            ok = False
            for _ in range(10):
                if fb == fa:
                    break
                c = b - fb * (b - a) / (fb - fa)
                if abs(c - R0) > 0.01:
                    break
                a, fa, b = b, fb, c
                fb = F(b)
                nev += 1
                if abs(b - a) < 1e-11:
                    ok = True
                    break
            if not ok or not (lo <= b < hi) or any(abs(b - e["R"]) < 1e-7 for e in mine):
                if ok and lo <= b < hi:
                    continue
                rejected.append({"subsector": sub_name(s), "start": R0, "reason": "secant did not settle in the window"})
                continue
            c1, c2 = coeffs(b, hcs, fr, par, eps)
            nev += 1
            agree, hecke, b11, b13 = tests(c1, c2)
            if agree > 1e-6 or hecke > 1e-4:
                rejected.append({"subsector": sub_name(s), "R": b, "agree": agree, "hecke": hecke})
                continue
            new = abs(abs(b11) - 1) < 1e-4 and abs(abs(b13) - 1) < 1e-4
            e = {"subsector": sub_name(s), "parity": "even" if par == 0 else "odd", "eps11": eps[0], "eps13": eps[1], "R": b,
                 "lambda": 0.25 + b * b, "new": bool(new), "agree": agree, "hecke": hecke, "b3": float(c1[1]), "b5": float(c1[2]),
                 "b11_sqrt11": b11, "b13_sqrt13": b13}
            mine.append(e)
            print(f"  {sub_name(s)}: R = {b:.10f}  lambda/(1/4) = {(0.25 + b * b) * 4:.6f}  {'new' if new else 'OLD'}  agree {agree:.1e}  "
                  f"hecke {hecke:.1e}  b11 sqrt11 {b11:+.6f}  b13 sqrt13 {b13:+.6f}", flush=True)
        found += mine
        print(f"  [{sub_name(s)} done: {len(merged)} candidates, {time.time() - T0:.0f} s]", flush=True)
    json.dump({"scan": os.path.basename(scanfile), "method": "find2", "horocycle fractions": [f1, f2], "window": [lo, hi], "Y0": Y0, "M": M,
               "found": sorted(found, key=lambda e: e["R"]), "rejected": rejected, "evaluations": nev, "seconds": time.time() - T0},
              open(out, "w"), indent=1)
    print(f"find2: {len(found)} forms, {len(rejected)} rejected, {nev} evaluations, {time.time() - T0:.0f} s")


# ------------------------------------------------------------------------------------------- Route T: the twisted trace formula
def _t16():
    return _load("mp16_trace", os.path.join(HERE, "..", "mp16", "mp16_trace.py"))


def W_tw(t16, Delta, t):
    """MP16's W_of(Delta, 143, 1), each term of index u = f signed by chi: +1 for odd t, (-1)^u for even t."""
    cd = t16.pari.coredisc(Delta, 1)
    D, c0 = int(cd[0]), int(cd[1])
    L = t16.L1(D)
    fac = [(int(p), int(e)) for p, e in zip(*t16.pari.factor(c0))] if c0 > 1 else []
    tot = tot_plain = 0.0
    for f in t16.divisors_from(fac):
        c = c0 // f
        term = L / f
        for p, _ in fac:
            if c % p == 0:
                term *= 1 - t16.kron(D, p) / p
        for p in (11, 13):
            term *= (p + 1) if f % p == 0 else (1 + (0 if c % p == 0 else t16.kron(D, p)))
        tot_plain += term
        tot += -term if (t % 2 == 0 and f % 2 == 1) else term
    return tot, tot_plain


def weights():
    t16 = _t16()
    z = np.load(os.path.join(HERE, "..", "mp17", "weights_143_1.npz"))
    hyp = z["hyp"].copy()                                   # rows (t, l, w = 2 W)
    out = hyp.copy()
    worst = 0.0
    for i in range(hyp.shape[0]):
        t = int(round(hyp[i, 0]))
        if t % 2:
            continue
        if hyp[i, 2] == 0.0:
            out[i, 2] = 0.0
            continue
        Wt, Wp = W_tw(t16, t * t - 4, t)
        worst = max(worst, abs(2 * Wp - hyp[i, 2]) / max(1.0, abs(hyp[i, 2])))
        out[i, 2] = 2 * Wt
        if i % 20000 == 0:
            print(f"  t = {t}  ({i}/{hyp.shape[0]}, {time.time() - T0:.0f} s)", flush=True)
    np.savez(os.path.join(HERE, "weights_tw.npz"), lmax=float(z["lmax"]), hyp=out)
    print(f"weights: {out.shape[0]} traces to length {float(z['lmax'])}; unsigned recomputation against MP17's cache, worst "
          f"relative difference {worst:.1e}; {time.time() - T0:.0f} s")


def geometric(h, g, W):
    ident = (VOL / (2 * math.pi)) * quad(lambda r: h(r) * r * math.tanh(math.pi * r), 0, 30, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    hyp = float(np.sum(W[:, 2] * g(W[:, 1])))
    par = -4 * float(g(0.0)) * math.log(2)
    return {"identity": ident, "hyperbolic": hyp, "cusps": par, "total": ident + hyp + par}


def smooth_count(R):
    return (VOL / (2 * math.pi)) * quad(lambda r: r * math.tanh(math.pi * r), 0, R)[0] - 4 * math.log(2) / math.pi * R


def trace(out, lists, Rtop=1.75):
    t16 = _t16()
    W = np.load(os.path.join(HERE, "weights_tw.npz"))["hyp"]
    forms = {}
    for f in lists:
        for e in json.load(open(f))["found"]:
            key = (e["subsector"], round(e["R"], 7))
            forms.setdefault(key, e)
    Rlist = sorted(e["R"] for e in forms.values())
    rows = []
    for r0, sig in ((0.0, 0.35), (0.0, 0.45), (0.374, 0.35), (0.6325, 0.35)):
        h, g = t16.testfun(r0, sig)
        geo = geometric(h, g, W)
        spec = float(sum(h(R) for R in Rlist))
        tau = quad(lambda r: 28 * r * h(r), Rtop, 40)[0]
        rows.append({"r0": r0, "sigma": sig, **{f"geometric {k}": v for k, v in geo.items()}, "spectral (list)": spec,
                     "list - geometry": spec - geo["total"], "tau_chi": tau, "pass": abs(spec - geo["total"]) <= tau})
        print(f"  r0 {r0:<6} sigma {sig}: geometry {geo['total']:+.9f} (identity {geo['identity']:.6f}, hyperbolic {geo['hyperbolic']:+.6f}, "
              f"cusps {geo['cusps']:+.6f})  list {spec:+.9f}  diff {spec - geo['total']:+.2e}  tau {tau:.1e}", flush=True)
    counts = {str(R): {"found": sum(1 for x in Rlist if x < R), "smooth": smooth_count(R)} for R in (0.6325, 1.0, 1.5)}
    json.dump({"forms": len(Rlist), "rows": rows, "counts": counts}, open(out, "w"), indent=1)
    print(f"trace: {len(Rlist)} forms; counts {counts}; {time.time() - T0:.0f} s")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "checks":
        checks()
    elif cmd == "scan":
        scan(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7])
    elif cmd == "find":
        find(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
    elif cmd == "find2":
        find2(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
    elif cmd == "find3":                                   # exhaustive candidates, forms already in the given lists skipped
        known = [(e["subsector"], e["R"]) for f in sys.argv[6:] for e in json.load(open(f))["found"]]
        find2(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5], exhaustive=True, known=known)
    elif cmd == "weights":
        weights()
    elif cmd == "trace":
        trace(sys.argv[2], sys.argv[3:])
