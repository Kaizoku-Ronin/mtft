#!/usr/bin/env python3
"""
Eisenstein congruences — the congruence primes of X_0(143)
==========================================================

MIT License — Copyright (c) 2026 Roger Tano
See LICENSE file for full terms.

In weight k and level 1, Herbrand-Ribet makes l | num(B_k / 2k)
equivalent to the existence of a weight-k cusp form congruent to the
Eisenstein series modulo l — for k = 12 that is the classical
tau(n) = sigma_11(n) mod 691.  The bottom triangular Faulhaber layer
a_{p,2} = 2p B_{p-1} (see mtft.combinatorial and
studies/triangular_layers.py) IS that Bernoulli numerator, so the
level-1 congruence modulus is computable from figurate geometry.

In weight 2 and level N no Bernoulli number is available, but the SAME
invariant is computable directly from the Manin/Merel machinery of
mtft.hecke as

    C(block) = gcd over good p of det((p+1) I - T_p | block),

together with the bad-prime conditions at p | N, since mtft.hecke
certifies that the Eisenstein complement carries T_p eigenvalue
exactly p + 1 at good p.  This module is that measurement.

Contents and epistemic classes
------------------------------
1.  eisenstein_modulus(block) — the congruence modulus (Pr, stable
    over all good primes to 43).  Results for X_0(143):

        143a1       (dim 2)   C = 1       no Eisenstein congruence
        11a1 ghost  (dim 4)   C = 5^4     modulus 5
        f2 quartic  (dim 8)   C = 7^2     norm-modulus 7
        f3 sextic   (dim 12)  C = 12^2    norm-modulus 12 = 2^2 * 3

    The determinant is the norm of (p + 1 - a_p) raised to the power
    dim/e, where e is the Galois-orbit degree: homology doubles each
    newform, and the level-11 form additionally enters twice through
    the two degeneracy maps.  Hence 5^4 for the ghost (e = 1,
    dim = 4) but 7^2 and 12^2 for the new orbits (e = 4, 6 with
    dim = 8, 12).

2.  E2 on the two elliptic blocks (Cert).  Weierstrass point counting
    over 27-28 primes gives gcd(p + 1 - a_p) = 5 for 11a1 and 1 for
    143a1, matching the Manin/Merel determinants exactly.  The 5 is
    independently Mazur's Eisenstein number for level 11,
    numerator((11 - 1)/12) = 5.

3.  STURM CERTIFICATE (Cert).  sturm_bound(143, 2) = 28.  The
    11a1-ghost congruence holds mod 5 for every n <= 28 — all primes
    2, 3, 5, 7, 11, 13, 17, 19, 23 including both bad primes, with
    prime powers supplied by the Hecke recursions — so the congruence
    is certified, not merely observed on a finite sample.

4.  BAD PRIMES (EXACT).  U_11 and U_13 come from Merel's matrices with
    off-P^1 images dropped (mtft.hecke).  On 143a1 both are the scalar
    -1; on the ghost U_11 is the scalar 1 = a_11(11a1) while U_13 has
    characteristic polynomial exactly (x^2 - 4x + 13)^2, the classical
    oldspace signature at the degeneracy prime.  The Eisenstein
    complement carries U_p eigenvalues in {1, p}.

Provenance: studies/eisenstein_congruences.py (7 gates).
"""

from __future__ import annotations

from fractions import Fraction as Fr
from math import gcd

from . import hecke as _h

__all__ = [
    "sturm_bound", "eisenstein_modulus", "congruence_census",
    "curve_ap", "hecke_on_block", "MODULI", "NORM_MODULI",
    "ORBIT_DEGREE", "GOOD_PRIMES",
]

GOOD_PRIMES = (2, 3, 5, 7, 17, 19, 23, 29, 31, 37, 41, 43)

MODULI = {"143a1": 1, "11a1_ghost": 625, "f2_quartic": 49,
          "f3_sextic": 144}
NORM_MODULI = {"143a1": 1, "11a1_ghost": 5, "f2_quartic": 7,
               "f3_sextic": 12}
ORBIT_DEGREE = {"143a1": 1, "11a1_ghost": 1, "f2_quartic": 4,
                "f3_sextic": 6}
_ALIAS = {"143a1": "ell", "11a1_ghost": "old", "f2_quartic": "q4",
          "f3_sextic": "q6"}


def sturm_bound(N=143, k=2):
    """Sturm bound (k/12)[SL_2(Z):Gamma_0(N)] for squarefree N."""
    idx = N
    n = N
    p = 2
    while p * p <= n:
        if n % p == 0:
            idx = idx * (p + 1) // p
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        idx = idx * (n + 1) // n
    return (k * idx) // 12


def _det(M):
    M = [list(r) for r in M]
    n = len(M)
    d = Fr(1)
    for c in range(n):
        pr = next((i for i in range(c, n) if M[i][c] != 0), None)
        if pr is None:
            return Fr(0)
        if pr != c:
            M[c], M[pr] = M[pr], M[c]
            d = -d
        d *= M[c][c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(c + 1, n):
            if M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    return d


def _restrict_block(Cp, Bv):
    d = len(Bv)
    n = len(Cp)
    img = [[sum(Cp[i][j] * v[j] for j in range(n)) for i in range(n)]
           for v in Bv]
    Aug = [[Bv[a][i] for a in range(d)] + [img[a][i] for a in range(d)]
           for i in range(n)]
    R, piv = _h._rref(Aug)
    if piv[:d] != list(range(d)):
        raise ValueError("supplied vectors are not independent")
    return [[R[a][d + b] for b in range(d)] for a in range(d)]


def hecke_on_block(name, p):
    """Matrix of T_p restricted to a named block.  EXACT."""
    key = _ALIAS.get(name, name)
    Cp = [list(r) for r in _h.cuspidal_hecke(p)]
    return _restrict_block(Cp, [list(v) for v in _h.blocks()[key]])


def eisenstein_modulus(name, primes=GOOD_PRIMES):
    """gcd over good p of det((p+1) I - T_p | block).  Pr.

    The Eisenstein complement has T_p eigenvalue exactly p + 1 at good
    p, so this determinant is the norm of (p + 1 - a_p) over the
    block; its gcd is the congruence modulus.  Divide out the
    homological doubling (take the square root) for the norm-modulus
    of the underlying newform orbit.
    """
    g = 0
    for p in primes:
        X = hecke_on_block(name, p)
        d = len(X)
        Np = _det([[Fr((p + 1) if i == j else 0) - X[i][j]
                    for j in range(d)] for i in range(d)])
        if Np.denominator != 1:
            raise ValueError("non-integral determinant")
        g = gcd(g, abs(int(Np)))
    return g


def congruence_census(primes=GOOD_PRIMES):
    """The congruence modulus of every block.  Pr."""
    return {name: eisenstein_modulus(name, primes) for name in _ALIAS}


def curve_ap(curve, pmax):
    """a_p by Weierstrass point counting — the independent route.

    `curve` is the coefficient tuple (a1, a2, a3, a4, a6); returns a
    dict over odd primes up to pmax.  Used to confirm the Manin/Merel
    determinants without sharing any step with them.
    """
    a1, a2, a3, a4, a6 = curve
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    primes = [p for p in range(3, pmax + 1)
              if all(p % q for q in range(2, int(p ** .5) + 1))]
    out = {}
    for p in primes:
        sq = bytearray(p)
        for k in range(p // 2 + 1):
            sq[(k * k) % p] = 1
        s = 0
        for x in range(p):
            g = (((4 * x + b2) * x + 2 * b4) * x + b6) % p
            if g:
                s += 1 if sq[g] else -1
        out[p] = -s
    return out
