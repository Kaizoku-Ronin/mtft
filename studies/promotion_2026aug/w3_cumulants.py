#!/usr/bin/env python3
"""
w3_cumulants.py — third cumulants of the Tano weight ensemble
============================================================================

MIT License — Copyright (c) 2026 Roger Tano

Sequel to w2_susceptibility.py.  The ensemble is P_beta(n) = n^-beta /
zeta(beta), beta > 1, with w_n = sum_{d|n} (log d)/d and l = log n.
Certified there: <w> = -zeta'(beta+1), <w^2> = T(beta) (closed),
Cov(w,l) = zeta''(beta+1), Fisher metric g of the exponential family
p ~ exp(lambda*w - beta*l).  This study delivers the full THIRD
cumulant tensor at lambda = 0 — the Amari–Chentsov cubic tensor of the
(beta, lambda) manifold, equal in natural coordinates to the complete
gradient of the Fisher metric (d_k g_ij = T_ijk).

THE NEW OBJECT.  kappa_www requires the triple-lcm Dirichlet series
   U(s) = sum_{d,e,f} (log d)(log e)(log f) / (d e f [d,e,f]^s),
with sum_n w_n^3 n^-s = zeta(s) U(s), hence <w^3>_beta = U(beta).

ROUTE 1 (Euler-product engine, exact local factors).  The local factor
sum_{a,b,c>=0} A^a B^b C^c S^max(a,b,c) has an exact closed rational
form by max-layer telescoping (Abel): with partial sums h_m, g_m of
q^a and a q^a, the layer products reduce to finitely many closed
geometric sums S_r(Qq^j) = sum m^r z^m, r <= 3.  Per-prime mixed
log-derivatives are then exact; tails over p > P are resummed exactly
through prime-zeta derivatives (Mobius over log zeta) applied to the
exact (n_q, n_Q) Fraction series of the log-derivative ratios.
U(s) = -M0 * (F3 + 3 F1 F2 + F1^3)   (Faa di Bruno on log M3).

ROUTE 2 (sieve).  Divisor sieve of w_n to N = 4e6; every moment gated
against its closed form within a proven tail bound using
w_n <= (log^2 n)/2 + 0.7 (validated with margin on the sieve).

PRE-REGISTRATIONS AND VERDICTS
------------------------------
P1 CONFIRMED: the pair specialization of the engine reproduces the
   certified T(s) of w2_susceptibility to ~1e-51 at s = 2.5, 3.
P2 CONFIRMED: the triple Euler-log decomposition does NOT terminate —
   the triple-lcm series is not a finite product of zetas (the pair
   case was; the triple case leaves the class).  The skeleton opens
     prod_i zeta(sigma+x_i) * prod_{i<j} zeta(v+x_i+x_j)/zeta(u+..)
     * zeta(s+3+X) * zeta(2s+3+X)^-3 * zeta(3s+3+X)^2 * ...  (X = x+y+z)
   with nonzero exponent counts at every degree 10, 11, 12.
P3 CONFIRMED: HAGEDORN TENSOR SPLIT.  kappa_lll = -(log zeta)'''(beta)
   diverges as exactly 2/(beta-1)^3 at the wall, while EVERY
   w-containing component stays finite:
     kappa_wll(beta) = -zeta'''(beta+1)            [EXACT, one line:
        d/d(-beta) of the certified Cov identity]
     kappa_wwl(beta) = -chi_w'(beta)
                     = -T'(beta) + 2 zeta'(beta+1) zeta''(beta+1)
     kappa_www(beta) = U + 3 T zeta'(beta+1) - 2 zeta'(beta+1)^3.
P4 CONFIRMED: cold skewness positive: gamma1_cold = 1.7232... > 0.

COLD CONSTANTS (beta -> 1+, engine dps 50, robustness-gated)
------------------------------------------------------------
  <w^3>_cold = U(1)      = 4.42947284842615649140232922679   (Pr)
  kappa3_cold            = 1.28839000968718081092964935357   (Pr)
  gamma1_cold (skewness) = 1.72320112367975760472641003784   (Pr)
  kappa_wwl_cold         = 2.41326066271876888229315220338   (Pr)
  kappa_wll_cold         = -zeta'''(2)
                         = 6.00014580284304486564394121754   (EXACT)

AG-D5 FILED AND DISMISSED: zeta'''(2) = -6.000145803 sits 1.46e-4
from -6.  Pole-echo heuristic: zeta ~ 1/(s-1) gives zeta''' ~
-6/(s-1)^4, which at s = 2 is exactly -6; the offset is the regular
part.  Proximity noted, no identity claimed.

THREE-ROUTE HIGHLIGHTS: kappa_wll agrees between (a) -zeta'''(beta+1)
via the tensor identity, (b) the closed moment combination
<w l^2> - 2<l><w l> - <w><l^2> + 2<w><l>^2 (Leibniz cancellation,
gap ~1e-52), and (c) the sieve.  kappa from the metric: curvature of
the (beta, lambda) manifold needs kappa_4 — queued, not claimed.

GATES G1-G8 as in main().  Runtime ~6 min (seven U evaluations).
Writes w3_cumulants_ledger.json next to itself.
"""

from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction as Fr
from itertools import product as iproduct

import numpy as np
from mpmath import mp, mpf, zeta, log, exp, diff, quad, inf

mp.dps = 50


# ── primes / mobius / prime-zeta log sums ───────────────────────────

def sieve_primes(P):
    s = np.ones(P + 1, bool)
    s[:2] = False
    for i in range(2, int(P ** .5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(x) for x in np.nonzero(s)[0]]


def _mobius(J):
    mu = np.ones(J + 1, np.int64)
    for p in sieve_primes(J):
        mu[p::p] *= -1
        mu[p * p::p * p] = 0
    return [int(x) for x in mu]


MU = _mobius(400)
zz = lambda t: zeta(t, derivative=1) / zeta(t)
_DK = {}


def primesum_logk(k, a):
    """sum_p (log p)^k p^-a via Mobius over log zeta (exact analytic)."""
    key = (k, str(a))
    if key in _DK:
        return _DK[key]
    a = mpf(a)
    J = max(1, int(300 / float(a)))
    tot = mpf(0)
    for j in range(1, J + 1):
        m = MU[j]
        if m == 0:
            continue
        x = mpf(j) * a
        if k == 0:
            tot += mpf(m) / j * log(zeta(x))
        elif k == 1:
            tot += m * zz(x)
        elif k == 2:
            tot += m * j * diff(zz, x, n=1)
        elif k == 3:
            tot += m * j * j * diff(zz, x, n=2)
    r = tot * (-1) ** k if k > 0 else tot
    _DK[key] = r
    return r


# ── exact per-prime W_k (max-layer telescoping) ─────────────────────

def rep_mul(X, Y):
    Z = {}
    for jx, cx in X.items():
        for jy, cy in Y.items():
            j = jx + jy
            cc = Z.setdefault(j, [mpf(0)] * 4)
            for a, ca in enumerate(cx):
                if ca == 0:
                    continue
                for b, cb in enumerate(cy):
                    if cb == 0 or a + b > 3:
                        continue
                    cc[a + b] += ca * cb
    return Z


def S_r(z):
    o = 1 - z
    return [1 / o, z / o ** 2, z * (1 + z) / o ** 3,
            z * (1 + 4 * z + z * z) / o ** 4]


def rep_sumQ(X, Q, q):
    tot = mpf(0)
    for j, cc in X.items():
        S = S_r(Q * q ** j)
        for r, c in enumerate(cc):
            if c != 0:
                tot += c * S[r]
    return tot


def W_all(q, Q, r):
    h = {0: [1 / (1 - q), 0, 0, 0], 1: [-q / (1 - q), 0, 0, 0]}
    g = {0: [q / (1 - q) ** 2, 0, 0, 0],
         1: [-q / (1 - q) ** 2, -q / (1 - q), 0, 0]}
    out = []
    for k in range(r + 1):
        X = {0: [mpf(1), 0, 0, 0]}
        for _ in range(k):
            X = rep_mul(X, g)
        for _ in range(r - k):
            X = rep_mul(X, h)
        out.append((1 - Q) * rep_sumQ(X, Q, q))
    return out


def brute_W(q, Q, r, M=100):
    outs = [mpf(0)] * (r + 1)
    for tup in iproduct(range(M), repeat=r):
        base = q ** sum(tup) * Q ** max(tup)
        outs[0] += base
        pr = 1
        for k in range(1, r + 1):
            pr *= tup[k - 1]
            if pr == 0:
                break
            outs[k] += pr * base
    return outs


# ── exact Fraction series in (n_q, n_Q) for tails ───────────────────

def series_W(r, OM):
    Ws = [dict() for _ in range(r + 1)]
    for tup in iproduct(range(OM + 1), repeat=r):
        nq, m = sum(tup), max(tup)
        if nq + m > OM:
            continue
        key = (nq, m)
        Ws[0][key] = Ws[0].get(key, 0) + 1
        pr = 1
        for k in range(1, r + 1):
            pr *= tup[k - 1]
            if pr == 0:
                break
            Ws[k][key] = Ws[k].get(key, 0) + pr
    return [{k: Fr(v) for k, v in W.items()} for W in Ws]


def ser_mul(X, Y, OM):
    Z = {}
    for kx, cx in X.items():
        for ky, cy in Y.items():
            k = (kx[0] + ky[0], kx[1] + ky[1])
            if k[0] + k[1] > OM:
                continue
            Z[k] = Z.get(k, Fr(0)) + cx * cy
    return Z


def ser_inv(X, OM):
    u = {k: v for k, v in X.items() if k != (0, 0)}
    out = {(0, 0): Fr(1)}
    pw = {(0, 0): Fr(1)}
    sgn = 1
    for _ in range(OM):
        pw = ser_mul(pw, u, OM)
        sgn = -sgn
        if not pw:
            break
        for k, v in pw.items():
            out[k] = out.get(k, Fr(0)) + sgn * v
    return out


def ser_log(X, OM):
    u = {k: v for k, v in X.items() if k != (0, 0)}
    out = {}
    pw = {(0, 0): Fr(1)}
    for n in range(1, OM + 1):
        pw = ser_mul(pw, u, OM)
        if not pw:
            break
        for k, v in pw.items():
            out[k] = out.get(k, Fr(0)) + Fr((-1) ** (n + 1), n) * v
    return {k: v for k, v in out.items() if v != 0}


_SER = {}


def build_R_series(r, OM):
    if (r, OM) in _SER:
        return _SER[(r, OM)]
    Ws = series_W(r, OM)
    inv0 = ser_inv(Ws[0], OM)
    R1 = ser_mul(Ws[1], inv0, OM)
    out = {"R1": R1, "logW0": ser_log(Ws[0], OM)}
    if r >= 2:
        R2raw = ser_mul(Ws[2], inv0, OM)
        R11 = ser_mul(R1, R1, OM)
        out["R2"] = {k: R2raw.get(k, Fr(0)) - R11.get(k, Fr(0))
                     for k in set(R2raw) | set(R11)}
    if r >= 3:
        R3raw = ser_mul(Ws[3], inv0, OM)
        R21 = ser_mul(ser_mul(Ws[2], inv0, OM), R1, OM)
        R111 = ser_mul(ser_mul(R1, R1, OM), R1, OM)
        out["R3"] = {k: (R3raw.get(k, Fr(0)) - 3 * R21.get(k, Fr(0))
                         + 2 * R111.get(k, Fr(0)))
                     for k in set(R3raw) | set(R21) | set(R111)}
    _SER[(r, OM)] = out
    return out


def F_sums(s, r, P=500, OM=20):
    s = mpf(s)
    pr = sieve_primes(P)
    ser = build_R_series(r, OM)
    keys = {"R1": (1, -1), "logW0": (0, 1)}
    if r >= 2:
        keys["R2"] = (2, 1)
    if r >= 3:
        keys["R3"] = (3, -1)
    res = {}
    Wcache = {}
    for p in pr:
        Wcache[p] = W_all(mpf(1) / p, mpf(p) ** (-s), r)
    for name, (k, sgn) in keys.items():
        tot = mpf(0)
        for p in pr:
            W = Wcache[p]
            lp = log(p)
            if name == "logW0":
                val = log(W[0])
            elif name == "R1":
                val = W[1] / W[0]
            elif name == "R2":
                val = W[2] / W[0] - (W[1] / W[0]) ** 2
            else:
                val = (W[3] / W[0] - 3 * W[2] * W[1] / W[0] ** 2
                       + 2 * (W[1] / W[0]) ** 3)
            tot += lp ** k * val
        for (nq, nQ), c in ser[name].items():
            if nq == 0 and nQ == 0:
                continue
            a = nq + s * nQ
            full = primesum_logk(k, a)
            part = sum(log(p) ** k * mpf(p) ** (-a) for p in pr)
            tot += mpf(c.numerator) / c.denominator * (full - part)
        res[name] = sgn * tot
    return res


def U3(s, P=500, OM=20):
    F = F_sums(s, 3, P, OM)
    return -exp(F["logW0"]) * (F["R3"] + 3 * F["R1"] * F["R2"]
                               + F["R1"] ** 3)


def T2_engine(s, P=500, OM=20):
    F = F_sums(s, 2, P, OM)
    return exp(F["logW0"]) * (F["R2"] + F["R1"] ** 2)


def T2_closed(s):
    s = mpf(s)
    sig, u, v = s + 1, 2 * s + 2, s + 2
    z = lambda a: zeta(a)
    z1 = lambda a: zeta(a, derivative=1)
    z2 = lambda a: zeta(a, derivative=2)
    C0 = z(sig) ** 2 / z(u)
    C1 = z(sig) ** 2 * z1(u) / z(u) ** 2 - z(sig) * z1(sig) / z(u)
    C2 = (z(sig) ** 2 * (2 * z1(u) ** 2 / z(u) ** 3 - z2(u) / z(u) ** 2)
          - 2 * z(sig) * z1(sig) * z1(u) / z(u) ** 2
          + z1(sig) ** 2 / z(u))
    return z2(v) * C0 - 2 * z1(v) * C1 + z(v) * C2


def tensor_closed(beta, Ub):
    b = mpf(beta)
    T = T2_closed(b)
    zp = zeta(b + 1, derivative=1)
    chi = T - zp ** 2
    return dict(
        T=T, chi=chi,
        k_www=Ub + 3 * T * zp - 2 * zp ** 3,
        k_wwl=-diff(T2_closed, b, n=1) + 2 * zp * zeta(b + 1, derivative=2),
        k_wll=-zeta(b + 1, derivative=3),
        k_lll=-diff(lambda t: log(zeta(t)), b, n=3),
    )


def main() -> int:
    t00 = time.time()
    ledger = {"study": "w3_cumulants", "gates": {}, "constants": {}}
    ok = True

    def gate(name, passed, **info):
        nonlocal ok
        ok &= bool(passed)
        ledger["gates"][name] = {"passed": bool(passed), **info}
        print(f"[{'PASS' if passed else 'FAIL'}] {name}  "
              + "  ".join(f"{k}={v}" for k, v in info.items()))

    # G1 engine soundness: closed W vs brute; pair engine vs certified T
    q, Q = mpf(1) / 3, mpf(1) / 7
    d3 = max(abs(a - b) for a, b in zip(W_all(q, Q, 3), brute_W(q, Q, 3)))
    d2 = max(abs(a - b) for a, b in zip(W_all(q, Q, 2),
                                        brute_W(q, Q, 2, M=160)))
    pg = [abs(T2_engine(s) - T2_closed(s)) for s in ("2.5", "3")]
    gate("G1_engine_vs_certified_pair",
         d3 < mpf("1e-40") and d2 < mpf("1e-40")
         and max(pg) < mpf("1e-40"),
         W_gap=f"{float(max(d2, d3)):.1e}",
         pair_T_gap=f"{float(max(pg)):.1e}")

    # G2 Euler-log decomposition census (exact Fractions, degree 12)
    DMAX = 12
    V = {}
    for a, b, c in iproduct(range(DMAX + 1), repeat=3):
        m = max(a, b, c)
        if a + b + c + m > DMAX:
            continue
        k = (a, b, c, m)
        V[k] = V.get(k, Fr(0)) + 1

    def pmul(X, Y):
        Z = {}
        for kx, cx in X.items():
            for ky, cy in Y.items():
                k = tuple(i + j for i, j in zip(kx, ky))
                if sum(k) > DMAX:
                    continue
                Z[k] = Z.get(k, Fr(0)) + cx * cy
        return {k: v for k, v in Z.items() if v != 0}

    u = {k: v for k, v in V.items() if k != (0, 0, 0, 0)}
    L = {}
    pw = {(0, 0, 0, 0): Fr(1)}
    for nn in range(1, DMAX + 1):
        pw = pmul(pw, u)
        if not pw:
            break
        for k, v in pw.items():
            L[k] = L.get(k, Fr(0)) + Fr((-1) ** (nn + 1), nn) * v
    work = {k: v for k, v in L.items() if v != 0}
    exps = {}
    while work:
        mu0 = min(work, key=lambda k: (sum(k), k))
        e = work[mu0]
        assert e.denominator == 1
        exps[mu0] = int(e)
        for r in range(1, DMAX // max(1, sum(mu0)) + 1):
            k = tuple(r * i for i in mu0)
            if sum(k) > DMAX:
                break
            work[k] = work.get(k, Fr(0)) - e * Fr(1, r)
        work = {k: v for k, v in work.items() if v != 0}
    bydeg = {}
    for k, e in exps.items():
        bydeg.setdefault(sum(k), []).append((k, e))
    skeleton_ok = (
        exps.get((1, 0, 0, 1)) == 1 and exps.get((0, 1, 1, 1)) == 1
        and exps.get((0, 1, 1, 2)) == -1 and exps.get((1, 1, 1, 1)) == 1
        and exps.get((1, 1, 1, 2)) == -3 and exps.get((1, 1, 1, 3)) == 2
        and exps.get((2, 1, 1, 2)) == -1)
    nonterm = all(len(bydeg.get(d, [])) > 0 for d in (10, 11, 12))
    gate("G2_decomposition_census", skeleton_ok and nonterm,
         skeleton="sigma^3 v^3 / u^3 * z(s+3) z(2s+3)^-3 z(3s+3)^2 ...",
         counts={d: len(bydeg.get(d, [])) for d in sorted(bydeg)},
         terminates=False)
    ledger["euler_log_exponents_deg_le_8"] = {
        str(k): e for k, e in sorted(exps.items()) if sum(k) <= 8}

    # sieve (route 2)
    N = 4_000_000
    w = np.zeros(N + 1)
    for d in range(2, N + 1):
        w[d::d] += np.log(d) / d
    n = np.arange(1, N + 1, dtype=np.float64)
    lv = np.log(n)
    wv = w[1:]
    wb_margin = float((wv - 0.5 * lv ** 2 - 0.7).max())

    def tail(a, bpow, s):
        f = lambda x: ((0.5 * mp.log(x) ** 2 + mpf("0.7")) ** a
                       * mp.log(x) ** bpow * x ** (-mpf(s)))
        return f(N) + quad(f, [N, inf])

    # G3 triple E2
    U = {s: U3(s) for s in ("4", "3.5", "3", "2", "1.5")}
    g3ok = wb_margin < 0
    info3 = {"w_bound_margin": f"{wb_margin:.3f}"}
    for s in ("4", "3.5"):
        sv = float(s)
        sieve_v = float((wv ** 3 * n ** (-sv)).sum())
        gap = abs(zeta(mpf(s)) * U[s] - mpf(sieve_v))
        tb = tail(3, 0, s) + mpf("5e-12")
        g3ok &= gap <= tb
        info3[f"s{s}"] = f"gap={float(gap):.1e}<=tail={float(tb):.1e}"
    gate("G3_triple_E2", g3ok, **info3)
    ledger["U_values"] = {s: str(v) for s, v in U.items()}

    # G4 mixed-moment E2 at beta = 3
    b = mpf(3)
    G = lambda s: -zeta(s) * zeta(s + 1, derivative=1)
    closed_m = {
        "w3": zeta(b) * U["3"],
        "w2": zeta(b) * T2_closed(b),
        "w2l": -diff(lambda s: zeta(s) * T2_closed(s), b, n=1),
        "wl": -diff(G, b, n=1),
        "wl2": diff(G, b, n=2),
        "l2": zeta(b, derivative=2),
        "l3": -zeta(b, derivative=3),
    }
    base = n ** (-3.0)
    sieve_m = {
        "w3": (wv ** 3 * base).sum(), "w2": (wv ** 2 * base).sum(),
        "w2l": (wv ** 2 * lv * base).sum(), "wl": (wv * lv * base).sum(),
        "wl2": (wv * lv ** 2 * base).sum(), "l2": (lv ** 2 * base).sum(),
        "l3": (lv ** 3 * base).sum(),
    }
    ab = {"w3": (3, 0), "w2": (2, 0), "w2l": (2, 1), "wl": (1, 1),
          "wl2": (1, 2), "l2": (0, 2), "l3": (0, 3)}
    g4ok = True
    info4 = {}
    for kk in closed_m:
        gap = abs(closed_m[kk] - mpf(float(sieve_m[kk])))
        tb = tail(*ab[kk], 3) + mpf("5e-12")
        g4ok &= gap <= tb
        info4[kk] = f"{float(gap):.1e}<={float(tb):.1e}"
    gate("G4_moment_E2_beta3", g4ok, **info4)

    # G5 multi-route cumulants
    g5ok = True
    info5 = {}
    for bt in ("3", "2.5"):
        bb = mpf(bt)
        mw = -zeta(bb + 1, derivative=1)
        ml = -zeta(bb, derivative=1) / zeta(bb)
        ml2 = zeta(bb, derivative=2) / zeta(bb)
        mwl = -diff(G, bb, n=1) / zeta(bb)
        mwl2 = diff(G, bb, n=2) / zeta(bb)
        r_a = -zeta(bb + 1, derivative=3)
        r_b = mwl2 - 2 * ml * mwl - mw * ml2 + 2 * mw * ml ** 2
        gap = abs(r_a - r_b)
        g5ok &= gap < mpf("1e-38")
        info5[f"k_wll_b{bt}"] = f"{float(gap):.1e}"
        mw2 = T2_closed(bb)
        mw2l = -diff(lambda s: zeta(s) * T2_closed(s), bb, n=1) / zeta(bb)
        r_c = mw2l - 2 * mw * mwl - ml * mw2 + 2 * mw ** 2 * ml
        r_d = (-diff(T2_closed, bb, n=1)
               + 2 * zeta(bb + 1, derivative=1) * zeta(bb + 1, derivative=2))
        gap2 = abs(r_c - r_d)
        g5ok &= gap2 < mpf("1e-30")
        info5[f"k_wwl_b{bt}"] = f"{float(gap2):.1e}"
    gate("G5_multiroute_cumulants", g5ok, **info5)

    # G6 Hagedorn split at the wall
    g6ok = True
    info6 = {}
    for bt in ("1.01", "1.001"):
        bb = mpf(bt)
        klll = -diff(lambda t: log(zeta(t)), bb, n=3)
        ratio = klll * (bb - 1) ** 3 / 2
        g6ok &= abs(ratio - 1) < mpf("1e-3")
        info6[f"lll_ratio_b{bt}"] = f"{float(ratio):.6f}"
    wall_wll = [abs(-zeta(mpf(bt) + 1, derivative=3)
                    - (-zeta(mpf(2), derivative=3)))
                for bt in ("1.01", "1.001")]
    g6ok &= wall_wll[1] < wall_wll[0]
    info6["wll_approach"] = f"{float(wall_wll[0]):.3f}->{float(wall_wll[1]):.3f}"
    gate("G6_hagedorn_split", g6ok, **info6)

    # G7 cold constants + robustness
    U1a = U3("1", P=500, OM=20)
    U1b = U3("1", P=300, OM=16)
    rob = abs(U1a - U1b)
    T1 = T2_closed(1)
    zp2 = zeta(mpf(2), derivative=1)
    chi_c = T1 - zp2 ** 2
    k3_c = U1a + 3 * T1 * zp2 - 2 * zp2 ** 3
    g1_c = k3_c / chi_c ** mpf("1.5")
    kwwl_c = (-diff(T2_closed, mpf(1), n=1)
              + 2 * zp2 * zeta(mpf(2), derivative=2))
    kwll_c = -zeta(mpf(2), derivative=3)
    gate("G7_cold_constants", rob < mpf("1e-30") and g1_c > 0,
         U1_robustness=f"{float(rob):.1e}",
         gamma1_cold=f"{float(g1_c):.12f}")
    for name, v in [("w3_cold", U1a), ("kappa3_cold", k3_c),
                    ("gamma1_cold", g1_c), ("kappa_wwl_cold", kwwl_c),
                    ("kappa_wll_cold", kwll_c), ("chi_cold", chi_c)]:
        ledger["constants"][name] = mp.nstr(v, 30)

    # G8 skewness curve positive on grid
    grid = {}
    g8ok = True
    for bt in ("1", "1.5", "2", "3"):
        Ub = U1a if bt == "1" else U[bt] if bt in U else U3(bt)
        tc = tensor_closed(bt, Ub)
        sk = tc["k_www"] / tc["chi"] ** mpf("1.5")
        grid[bt] = float(sk)
        g8ok &= sk > 0
    gate("G8_skewness_positive", g8ok, skew_grid=grid)
    ledger["tensor_beta3"] = {k: mp.nstr(v, 25)
                              for k, v in tensor_closed(3, U["3"]).items()}

    ledger["all_passed"] = ok
    ledger["runtime_s"] = round(time.time() - t00, 1)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "w3_cumulants_ledger.json")
    with open(out, "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print(f"\nledger -> {out}   [{ledger['runtime_s']}s]")
    print("ALL GATES PASS" if ok else "GATE FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
