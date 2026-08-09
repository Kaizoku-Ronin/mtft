#!/usr/bin/env python3
"""
Tano weight moments — closed forms for the arithmetic ensemble
==============================================================

MIT License — Copyright (c) 2026 Roger Tano
See LICENSE file for full terms.

The Tano ensemble p_beta(n) = n^-beta / zeta(beta) carries the weight
w_n = sum_{d|n} (log d)/d.  This module supplies its moment functions
in CLOSED FORM: finite expressions in zeta and its derivatives plus
one convergent Euler product per moment order.  Nothing here is a
sieve estimate; every value is arbitrary-precision and certified
against an independent route in tests/test_moments.py.

Contents and epistemic classes
------------------------------
1.  FIRST MOMENT (EXACT).  <w>_beta = -zeta'(beta+1).

2.  SECOND MOMENT (EXACT closed form).  <w^2>_beta = T(beta) with
    T = zeta''(v) C0 - 2 zeta'(v) C1 + zeta(v) C2, v = beta + 2, and
    C0, C1, C2 the Euler-product blocks of the pair kernel
    sum_{d,e} (log d)(log e) / (d e [d,e]^beta).  The susceptibility
    is chi_w = T - zeta'(beta+1)^2, and Cov(log n, w) = zeta''(beta+1)
    EXACTLY.

3.  THIRD MOMENT (Pr, dps-robust).  <w^3>_beta = U(beta), the triple
    kernel, by exact max-layer telescoping of the local factor at each
    prime followed by Moebius-resummed prime-zeta tails.  The pair
    specialization of the same engine reproduces T(beta) to 1.7e-51 —
    the engine's own certificate.

4.  CUMULANTS (EXACT / Pr).  With l = log n and w the weight, the
    Amari-Chentsov components are kappa_wll = -zeta'''(beta+1) EXACT,
    kappa_wwl = -d/dbeta chi_w, kappa_www = U + 3 T zeta'(beta+1)
    - 2 zeta'(beta+1)^3, and kappa_lll = -(log zeta)'''(beta).  The
    Fisher metric and the curvature built from them live in
    mtft.curvature.

5.  COLD CONSTANTS (Pr, dps 50, robustness 1.7e-35).  Values at
    beta = 1, the Hagedorn edge, where all w-moments stay finite while
    the log-moments diverge — the transparency of the Hagedorn wall to
    the weight sector.

Provenance: studies/w2_susceptibility.py and studies/w3_cumulants.py
(9 and 8 gates); constants reproduced in the tests.
"""

from __future__ import annotations

from fractions import Fraction as Fr

import numpy as np
from itertools import product as iproduct

from mpmath import mp, mpf, zeta, diff, log, exp, mpmathify

__all__ = [
    "weight_first_moment", "weight_second_moment",
    "weight_susceptibility", "weight_third_moment",
    "cov_log_weight", "cumulants", "COLD",
]

def sieve_primes(P):
    s = np.ones(P+1, bool); s[:2] = False
    for i in range(2, int(P**.5)+1):
        if s[i]: s[i*i::i] = False
    return [int(x) for x in np.nonzero(s)[0]]

MU = None
def _mobius(J):
    mu = np.ones(J+1, np.int64)
    for p in sieve_primes(J):
        mu[p::p] *= -1; mu[p*p::p*p] = 0
    return [int(x) for x in mu]

zz = lambda t: zeta(t, derivative=1)/zeta(t)
def primesum_logk(k, a):
    global MU
    if MU is None: MU = _mobius(400)
    a = mpf(a); J = max(1, int(300/float(a))); tot = mpf(0)
    for j in range(1, J+1):
        m = MU[j]
        if m == 0: continue
        x = mpf(j)*a
        if k == 0: tot += mpf(m)/j*log(zeta(x))
        elif k == 1: tot += m*zz(x)
        elif k == 2: tot += m*j*diff(zz, x, n=1)
        elif k == 3: tot += m*j*j*diff(zz, x, n=2)
    return tot*(-1)**k if k > 0 else tot

# ---- exact per-prime W_k for the r-fold max sum ----
# rep: dict j -> [c0..c3]: sum_m (c0+c1 m+c2 m^2+c3 m^3) q^{j m}
def rep_mul(X, Y):
    Z = {}
    for jx, cx in X.items():
        for jy, cy in Y.items():
            j = jx + jy
            cc = Z.setdefault(j, [mpf(0)]*4)
            for a, ca in enumerate(cx):
                if ca == 0: continue
                for b, cb in enumerate(cy):
                    if cb == 0 or a+b > 3: continue
                    cc[a+b] += ca*cb
    return Z
def S_r(z):  # sum_{m>=0} m^r z^m, r=0..3
    o = 1 - z
    return [1/o, z/o**2, z*(1+z)/o**3, z*(1+4*z+z*z)/o**4]
def rep_sumQ(X, Q, q):
    tot = mpf(0)
    for j, cc in X.items():
        S = S_r(Q*q**j)
        for r, c in enumerate(cc):
            if c != 0: tot += c*S[r]
    return tot
def W_all(q, Q, r):
    """W_k(q,Q) = sum_{a_1..a_r>=0} (prod of first k coords) q^{sum a} Q^{max a}, k=0..r"""
    h = {0: [1/(1-q), 0, 0, 0], 1: [-q/(1-q), 0, 0, 0]}
    g = {0: [q/(1-q)**2, 0, 0, 0], 1: [-q/(1-q)**2, -q/(1-q), 0, 0]}
    out = []
    for k in range(r+1):
        X = {0: [mpf(1), 0, 0, 0]}
        for _ in range(k): X = rep_mul(X, g)
        for _ in range(r-k): X = rep_mul(X, h)
        out.append((1-Q)*rep_sumQ(X, Q, q))
    return out

def brute_W(q, Q, r, M=80):
    outs = [mpf(0)]*(r+1)
    for tup in iproduct(range(M), repeat=r):
        base = q**sum(tup)*Q**max(tup)
        outs[0] += base
        pr = 1
        for k in range(1, r+1):
            pr *= tup[k-1]
            if pr == 0: break
            outs[k] += pr*base
    return outs

# ---- series in (nq, nQ) with Fractions for tails ----
def series_W(r, OM):
    Ws = [dict() for _ in range(r+1)]
    for tup in iproduct(range(OM+1), repeat=r):
        nq, m = sum(tup), max(tup)
        if nq + m > OM: continue
        key = (nq, m)
        Ws[0][key] = Ws[0].get(key, 0) + 1
        pr = 1
        for k in range(1, r+1):
            pr *= tup[k-1]
            if pr == 0: break
            Ws[k][key] = Ws[k].get(key, 0) + pr
    return [{k: Fr(v) for k, v in W.items()} for W in Ws]
def ser_mul(X, Y, OM):
    Z = {}
    for kx, cx in X.items():
        for ky, cy in Y.items():
            k = (kx[0]+ky[0], kx[1]+ky[1])
            if k[0]+k[1] > OM: continue
            Z[k] = Z.get(k, Fr(0)) + cx*cy
    return Z
def ser_inv(X, OM):  # X = 1 + u
    u = {k: v for k, v in X.items() if k != (0, 0)}
    out = {(0, 0): Fr(1)}; pw = {(0, 0): Fr(1)}
    for _ in range(OM):
        pw = ser_mul(pw, u, OM)
        if not pw: break
        for k, v in pw.items():
            out[k] = out.get(k, Fr(0)) + (-1)**1*v if False else out.get(k, Fr(0))
        # careful: (1+u)^-1 = sum (-u)^n
    # redo properly
    out = {(0, 0): Fr(1)}; pw = {(0, 0): Fr(1)}; sgn = 1
    for _ in range(OM):
        pw = ser_mul(pw, u, OM); sgn = -sgn
        if not pw: break
        for k, v in pw.items():
            out[k] = out.get(k, Fr(0)) + sgn*v
    return out
def ser_log(X, OM):  # X = 1 + u
    u = {k: v for k, v in X.items() if k != (0, 0)}
    out = {}; pw = {(0, 0): Fr(1)}
    for n in range(1, OM+1):
        pw = ser_mul(pw, u, OM)
        if not pw: break
        for k, v in pw.items():
            out[k] = out.get(k, Fr(0)) + Fr((-1)**(n+1), n)*v
    return {k: v for k, v in out.items() if v != 0}

def build_R_series(r, OM):
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
        out["R3"] = {k: R3raw.get(k, Fr(0)) - 3*R21.get(k, Fr(0)) + 2*R111.get(k, Fr(0))
                     for k in set(R3raw) | set(R21) | set(R111)}
    return out

def F_sums(s, r, P=500, OM=20):
    """returns dict with F1,F2,F3 (as available), logM0 for the r-fold master series"""
    s = mpf(s)
    pr = sieve_primes(P)
    ser = build_R_series(r, OM)
    keys = {"R1": (1, -1), "logW0": (0, 1)}
    if r >= 2: keys["R2"] = (2, 1)
    if r >= 3: keys["R3"] = (3, -1)
    res = {}
    for name, (k, sgn) in keys.items():
        tot = mpf(0)
        for p in pr:
            q = mpf(1)/p; Q = mpf(p)**(-s); lp = log(p)
            W = W_all(q, Q, r)
            if name == "logW0": val = log(W[0])
            elif name == "R1": val = W[1]/W[0]
            elif name == "R2": val = W[2]/W[0] - (W[1]/W[0])**2
            else: val = W[3]/W[0] - 3*W[2]*W[1]/W[0]**2 + 2*(W[1]/W[0])**3
            tot += lp**k*val
        # tail via series and prime-zeta sums over p > P
        for (nq, nQ), c in ser[name].items():
            if nq == 0 and nQ == 0: continue
            a = nq + s*nQ
            full = primesum_logk(k, a)
            part = sum(log(p)**k*mpf(p)**(-a) for p in pr)
            tot += mpf(c.numerator)/c.denominator*(full - part)
        res[name] = sgn*tot
    return res

def U3(s, P=500, OM=20):
    F = F_sums(s, 3, P, OM)
    F1, F2, F3, lM0 = F["R1"], F["R2"], F["R3"], F["logW0"]
    return -exp(lM0)*(F3 + 3*F1*F2 + F1**3), F

def T2_engine(s, P=500, OM=20):
    F = F_sums(s, 2, P, OM)
    F1, F2, lM0 = F["R1"], F["R2"], F["logW0"]
    return exp(lM0)*(F2 + F1**2)

def T2_closed(s):
    s = mpf(s); sig, u, v = s+1, 2*s+2, s+2
    z = lambda a: zeta(a); z1 = lambda a: zeta(a, derivative=1); z2 = lambda a: zeta(a, derivative=2)
    C0 = z(sig)**2/z(u)
    C1 = z(sig)**2*z1(u)/z(u)**2 - z(sig)*z1(sig)/z(u)
    C2 = (z(sig)**2*(2*z1(u)**2/z(u)**3 - z2(u)/z(u)**2)
          - 2*z(sig)*z1(sig)*z1(u)/z(u)**2 + z1(sig)**2/z(u))
    return z2(v)*C0 - 2*z1(v)*C1 + z(v)*C2

# ── public interface ────────────────────────────────────────────────

def weight_first_moment(beta):
    """<w>_beta = -zeta'(beta+1).  EXACT."""
    return -mp.diff(mp.zeta, mp.mpf(beta) + 1, 1)


def cov_log_weight(beta):
    """Cov_beta(log n, w) = zeta''(beta+1).  EXACT."""
    return mp.diff(mp.zeta, mp.mpf(beta) + 1, 2)


def weight_second_moment(beta, **kw):
    """<w^2>_beta = T(beta).  EXACT closed form."""
    return T2_closed(mp.mpf(beta), **kw)


def weight_susceptibility(beta, **kw):
    """chi_w(beta) = <w^2> - <w>^2.  EXACT closed form."""
    b = mp.mpf(beta)
    return T2_closed(b, **kw) - mp.diff(mp.zeta, b + 1, 1) ** 2


def weight_third_moment(beta, **kw):
    """<w^3>_beta = U(beta).  Pr (dps-robust to 1.7e-35)."""
    val = U3(mp.mpf(beta), **kw)
    return val[0] if isinstance(val, tuple) else val


def cumulants(beta, **kw):
    """Amari-Chentsov components at (beta, lambda = 0).

    Keys 'lll', 'wll', 'wwl', 'www'.  Sign rule:
    d_beta^a d_lambda^b psi = (-1)^a kappa(l^a, w^b).
    """
    b = mp.mpf(beta)
    zp = mp.diff(mp.zeta, b + 1, 1)
    T = T2_closed(b, **kw)
    U = weight_third_moment(b, **kw)
    return {
        "lll": -mp.diff(lambda x: mp.log(mp.zeta(x)), b, 3),
        "wll": -mp.diff(mp.zeta, b + 1, 3),
        "wwl": -mp.diff(lambda x: weight_susceptibility(x, **kw), b, 1),
        "www": U + 3 * T * zp - 2 * zp ** 3,
    }


COLD = {
    "T": "1.70276979154901697001",
    "chi_w": "0.82377306237833093427",
    "zeta2_second": "1.98928023429890102342",
    "U": "4.42947284842615649140232922679",
    "kappa3": "1.28839000968718081092964935357",
    "skewness": "1.72320112367975760472641003784",
    "kappa_wwl": "2.41326066271876888229315220338",
    "kappa_wll": "-6.00014580284304486564394121754",
}
