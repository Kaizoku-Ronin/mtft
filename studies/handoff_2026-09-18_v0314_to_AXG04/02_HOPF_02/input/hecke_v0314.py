#!/usr/bin/env python3
"""
Manin symbols and Hecke blocks of X_0(143)
==========================================

MIT License — Copyright (c) 2026 Roger Tano
See LICENSE file for full terms.

The Farey skeleton of X_0(143) and its modular symbols live on the
SAME index set, P^1(Z/143): Manin symbols are the oriented edges of
the tessellation (S-orbits with sign), the 2-term relation is
orientation reversal, and the 3-term relation is exactly the triangle
boundary matrix.  Hence

    modular symbols = Q^84 / im d2 = H_1(X, cusps)   (dim 29),
    cuspidal part   = ker(cusp boundary)             (dim 26 = H_1).

Everything below is exact rational or integer linear algebra;
floating point appears only as a guess generator whose output is
certified exactly.

Contents and epistemic classes
------------------------------
1.  MODEL (EXACT).  rank d2 = 55, quotient 29, cuspidal 26; the cusp
    boundary annihilates im d2.

2.  HECKE (EXACT).  Merel's matrices (det p, a > b >= 0, d > c >= 0)
    give T_p for ALL p once images falling outside P^1(Z/143) are
    dropped (those Manin symbols are zero) — including the bad primes
    11 and 13, where the naive [[1,r],[0,p]] set fails to descend.
    The operators pairwise commute and act on the 3-dimensional
    Eisenstein complement with eigenvalue exactly p + 1 at good p.

3.  THE FOUR PARTICLES (Pr / EXACT).  By exact integer polynomial
    multiplication,

        charpoly(T_2 | H_1) = x^2 (x+2)^4 g4(x)^2 h6(x)^2,
        g4 = x^4 - 3x^3 - x^2 + 5x + 1        (irreducible)
        h6 = x^6 - 10x^4 + 2x^3 + 24x^2 - 7x - 12  (irreducible)

    giving blocks of dimensions [2, 4, 8, 12] = 26:
      ell  (2)  ker T_2        = 143a1, a_p matching point counts
      old  (4)  ker(T_2 + 2)   = the LEVEL-11 OLDSPACE, the ghost of
                                 X_0(11) inside 143 = 11 * 13
      q4   (8)  ker g4(T_2)    = the quartic new orbit f2
      q6  (12)  ker h6(T_2)    = the sextic new orbit f3
    The corpus statement "three Galois orbits [1, 4, 6]" is the
    NEWSPACE and is re-derived here from scratch; the homology
    additionally carries old-11a with multiplicity two.

    CORPUS RECONCILIATION.  mtft.x0_143.FIELD_POLY_F2 and g4 define
    the SAME quartic field: both have polynomial discriminant 1957
    (the stored field discriminant) and identical Frobenius
    factorization patterns at every prime below 400.  The stored
    polynomial is a reduced generator; g4 is the minimal polynomial of
    a_2 itself.  Likewise FIELD_POLY_F3 equals h6(-x) exactly — the
    opposite Hecke sign convention.  Both records are correct and
    complementary.  (Separately: mtft.x0_143.A2_COMPLEX remains
    exposed but was struck by corpus correction CC-01; weight-2 Hecke
    eigenvalues are totally real.)

4.  STAR INVOLUTION (EXACT).  iota*(c:d) = (-c:d) descends, commutes
    with every Hecke operator, and splits the blocks (1,1), (2,2),
    (4,4), (6,6) — the real Hodge split of each particle.

5.  HARMONIC DENSITIES (Pr).  Each block gets its unique harmonic
    representatives (orthogonal to im d2 in R^84); the presence
    density rho = diagonal of the orthogonal projector is
    basis-independent, iota*-invariant, and has trace exactly the
    block dimension.  Every cuspidal density vanishes identically on
    the self-loop edge — the particles avoid the one edge carrying
    the skeleton's graph noncommutativity.

Provenance: studies/x0143_hecke_particles.py (9 gates) and
studies/x0143_ribbon_embedding.py (10 gates).
"""

from __future__ import annotations

from fractions import Fraction as Fr
from functools import lru_cache
from math import gcd

__all__ = [
    "LEVEL", "G4", "H6", "BLOCK_DIMS", "model", "merel",
    "hecke_matrix", "cuspidal_hecke", "blocks", "star_involution",
    "harmonic_density", "charpoly_T2_factorization",
]

LEVEL = 143
G4 = [1, 5, -1, -3, 1]                       # low -> high
H6 = [-12, -7, 24, 2, -10, 0, 1]             # low -> high
BLOCK_DIMS = {"ell": 2, "old": 4, "q4": 8, "q6": 12}
BAD_PRIMES = (11, 13)


# ── exact linear algebra helpers ────────────────────────────────────

def _rref(M):
    M = [[Fr(x) for x in row] for row in M]
    R, C = len(M), len(M[0])
    piv, r = [], 0
    for c in range(C):
        pr = next((i for i in range(r, R) if M[i][c] != 0), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(R):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == R:
            break
    return M, piv


def _nullspace(M):
    R, piv = _rref(M)
    C = len(M[0])
    ps = set(piv)
    out = []
    for fc in [c for c in range(C) if c not in ps]:
        v = [Fr(0)] * C
        v[fc] = Fr(1)
        for r_, c_ in enumerate(piv):
            v[c_] = -R[r_][fc]
        out.append(v)
    return out


def _orbits(s):
    seen, out = set(), []
    for i in range(len(s)):
        if i in seen:
            continue
        o, j = [], i
        while j not in seen:
            seen.add(j)
            o.append(j)
            j = s[j]
        out.append(o)
    return out


def merel(p):
    """Merel's matrices of determinant p with a > b >= 0, d > c >= 0."""
    out = []
    for a in range(1, p + 1):
        for b in range(0, a):
            for c in range(0, p + 1):
                num = p + b * c
                if num % a:
                    continue
                d = num // a
                if d > c and a * d - b * c == p:
                    out.append((a, b, c, d))
    return out


@lru_cache(maxsize=1)
def model():
    """Build the exact P^1(Z/143) Manin model.  Memoized.

    Returns a dict with the flag set, the S/T/R actions, the edge and
    triangle data, the boundary matrices, the 29-dimensional quotient
    map and the 26-dimensional cuspidal basis.
    """
    N = LEVEL
    units = [u for u in range(1, N) if gcd(u, N) == 1]

    def valid(c, d):
        return not (c % 11 == 0 and d % 11 == 0) and \
               not (c % 13 == 0 and d % 13 == 0)

    def canon(c, d):
        return min(((c * u) % N, (d * u) % N) for u in units)

    P1 = sorted({canon(c, d) for c in range(N) for d in range(N)
                 if valid(c, d)})
    idx = {p: i for i, p in enumerate(P1)}
    perm = lambda f: [idx[canon(*f(*p))] for p in P1]
    sS = perm(lambda c, d: (d % N, (-c) % N))
    sT = perm(lambda c, d: (c, (c + d) % N))
    sR = perm(lambda c, d: (d % N, (d - c) % N))
    iota = perm(lambda c, d: ((-c) % N, d))

    tris = _orbits(sR)
    tri_of = {}
    for ti, t in enumerate(tris):
        for f in t:
            tri_of[f] = ti
    fans = sorted(_orbits(sT), key=len)
    cusp_of = {}
    for k, o in enumerate(fans):
        for f in o:
            cusp_of[f] = k
    Eorb = _orbits(sS)
    E = len(Eorb)
    erep, eid, esign = {}, {}, {}
    for k, o in enumerate(Eorb):
        r = min(o)
        erep[k] = r
        for f in o:
            eid[f] = k
            esign[f] = 1 if f == r else -1

    D2 = [[0] * 56 for _ in range(E)]
    for ti, t in enumerate(tris):
        for f in t:
            D2[eid[f]][ti] += esign[f]
    DEL = [[0] * 4 for _ in range(E)]
    for k in range(E):
        x = erep[k]
        DEL[k][cusp_of[x]] += 1
        DEL[k][cusp_of[sS[x]]] -= 1

    _, cols = _rref([[D2[e][t] for t in range(56)] for e in range(E)])
    D2r = [[D2[e][t] for t in cols] for e in range(E)]
    _, pivE = _rref([[D2r[e][j] for e in range(E)]
                     for j in range(len(cols))])
    free = [e for e in range(E) if e not in set(pivE)]

    B = [[Fr(0)] * E for _ in range(E)]
    for e in range(E):
        for j in range(len(cols)):
            B[e][j] = Fr(D2r[e][j])
    for j, e in enumerate(free):
        B[e][len(cols) + j] = Fr(1)
    M = [row[:] + [Fr(1) if i == j else Fr(0) for j in range(E)]
         for i, row in enumerate(B)]
    for c in range(E):
        pr = next(i for i in range(c, E) if M[i][c] != 0)
        M[c], M[pr] = M[pr], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(E):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    Binv = [row[E:] for row in M]
    nq = len(free)

    DB = [[Fr(DEL[free[j]][c]) for j in range(nq)] for c in range(4)]
    Mr, pivK = _rref(DB)
    freeK = [j for j in range(nq) if j not in set(pivK)]
    K = []
    for fj in freeK:
        v = [Fr(0)] * nq
        v[fj] = Fr(1)
        for r_, c_ in enumerate(pivK):
            v[c_] = -Mr[r_][fj]
        K.append(v)

    return dict(N=N, P1=P1, idx=idx, canon=canon, sS=sS, sT=sT, sR=sR,
                iota=iota, tris=tris, tri_of=tri_of, fans=fans, cusp_of=cusp_of,
                E=E, erep=erep, eid=eid, esign=esign, D2=D2, DEL=DEL,
                cols=cols, D2r=D2r, free=free, Binv=Binv, nq=nq,
                K=K, freeK=freeK)


def _qcoords(m, v):
    E, Binv, ncols = m["E"], m["Binv"], len(m["cols"])
    x = [sum(Binv[i][k] * v[k] for k in range(E) if v[k])
         for i in range(E)]
    return x[ncols:]


def _raw(m, mats):
    N, E = m["N"], m["E"]
    P1, idx, canon = m["P1"], m["idx"], m["canon"]
    erep, eid, esign = m["erep"], m["eid"], m["esign"]
    T = [[0] * E for _ in range(E)]
    for k in range(E):
        c, d = P1[erep[k]]
        for (a, b, cc, dd) in mats:
            u, v = (c * a + d * cc) % N, (c * b + d * dd) % N
            if gcd(gcd(u, v), N) != 1:
                continue                  # symbol vanishes off P^1
            y = idx[canon(u, v)]
            T[eid[y]][k] += esign[y]
    return T


def _quotient(m, T):
    nq, free, E = m["nq"], m["free"], m["E"]
    cs = [_qcoords(m, [Fr(T[e][fj]) for e in range(E)]) for fj in free]
    return [[cs[j][i] for j in range(nq)] for i in range(nq)]


def _restrict26(m, Q):
    K, freeK, nq = m["K"], m["freeK"], m["nq"]
    out = []
    for k in K:
        img = [sum(Q[i][j] * k[j] for j in range(nq))
               for i in range(nq)]
        out.append([img[fj] for fj in freeK])
    d = len(K)
    return [[out[j][i] for j in range(d)] for i in range(d)]


def hecke_matrix(p):
    """T_p on the 29-dimensional modular-symbol quotient.  EXACT."""
    m = model()
    return _quotient(m, _raw(m, merel(p)))


@lru_cache(maxsize=32)
def cuspidal_hecke(p):
    """T_p on the 26-dimensional cuspidal homology.  EXACT.

    Valid for every prime, including the bad primes 11 and 13 where
    this is U_p.
    """
    m = model()
    return tuple(tuple(r) for r in
                 _restrict26(m, _quotient(m, _raw(m, merel(p)))))


def _as_list(M):
    return [list(r) for r in M]


def _poly_of_op(coefs, C):
    n = len(C)
    Mo = [[Fr(0)] * n for _ in range(n)]
    P = [[Fr(1) if i == j else Fr(0) for j in range(n)]
         for i in range(n)]
    for c in coefs:
        if c:
            for i in range(n):
                for j in range(n):
                    Mo[i][j] += c * P[i][j]
        P = [[sum(P[i][k] * C[k][j] for k in range(n))
              for j in range(n)] for i in range(n)]
    return Mo


@lru_cache(maxsize=1)
def blocks():
    """The four Hecke blocks as exact bases of the 26-dim homology.

    Returns dict name -> tuple of basis vectors (each a 26-tuple of
    Fractions), with dimensions [2, 4, 8, 12].
    """
    C2 = _as_list(cuspidal_hecke(2))
    n = len(C2)
    out = {
        "ell": _nullspace(C2),
        "old": _nullspace([[C2[i][j] + (2 if i == j else 0)
                            for j in range(n)] for i in range(n)]),
        "q4": _nullspace(_poly_of_op(G4, C2)),
        "q6": _nullspace(_poly_of_op(H6, C2)),
    }
    return {k: tuple(tuple(v) for v in vs) for k, vs in out.items()}


def charpoly_T2_factorization():
    """The certified factorization data of charpoly(T_2 | H_1)."""
    return {"factors": "x^2 (x+2)^4 g4^2 h6^2", "g4": G4, "h6": H6,
            "dims": dict(BLOCK_DIMS)}


@lru_cache(maxsize=1)
def star_involution():
    """iota* on the 26-dimensional cuspidal homology.  EXACT."""
    m = model()
    E, iota, erep, eid, esign = (m["E"], m["iota"], m["erep"],
                                 m["eid"], m["esign"])
    IS = [[0] * E for _ in range(E)]
    for k in range(E):
        y = iota[erep[k]]
        IS[eid[y]][k] += esign[y]
    return tuple(tuple(r) for r in
                 _restrict26(m, _quotient(m, IS)))


def harmonic_density(name):
    """Presence density rho of a block on the 84 edges.  Pr.

    rho is the diagonal of the orthogonal projector onto the block's
    harmonic representatives: basis-independent, iota*-invariant, and
    with trace exactly the block dimension.
    """
    m = model()
    E, D2, K, free, nq = (m["E"], m["D2"], m["K"], m["free"], m["nq"])
    G56 = [[sum(D2[e][i] * D2[e][j] for e in range(E))
            for j in range(56)] for i in range(56)]

    def harmonic(v26):
        v = [sum(K[a][j] * v26[a] for a in range(len(K)))
             for j in range(nq)]
        v84 = [Fr(0)] * E
        for j, e in enumerate(free):
            v84[e] = v[j]
        rhs = [sum(Fr(D2[e][i]) * v84[e] for e in range(E))
               for i in range(56)]
        Aug = [[Fr(G56[i][j]) for j in range(56)] + [rhs[i]]
               for i in range(56)]
        R, pv = _rref(Aug)
        alpha = [Fr(0)] * 56
        for r_, c_ in enumerate(pv):
            if c_ < 56:
                alpha[c_] = R[r_][56]
        return [v84[e] - sum(Fr(D2[e][i]) * alpha[i]
                             for i in range(56)) for e in range(E)]

    W = [harmonic(list(v)) for v in blocks()[name]]
    d = len(W)
    Gm = [[sum(W[a][e] * W[b][e] for e in range(E)) for b in range(d)]
          for a in range(d)]
    Aug = [row[:] + [Fr(1) if i == j else Fr(0) for j in range(d)]
           for i, row in enumerate(Gm)]
    for c in range(d):
        pr = next(i for i in range(c, d) if Aug[i][c] != 0)
        Aug[c], Aug[pr] = Aug[pr], Aug[c]
        pv = Aug[c][c]
        Aug[c] = [x / pv for x in Aug[c]]
        for i in range(d):
            if i != c and Aug[i][c] != 0:
                f = Aug[i][c]
                Aug[i] = [a - f * b for a, b in zip(Aug[i], Aug[c])]
    Gi = [row[d:] for row in Aug]
    return [sum(W[a][e] * Gi[a][b] * W[b][e]
                for a in range(d) for b in range(d)) for e in range(E)]
