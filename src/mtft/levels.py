"""mtft.levels — level-generic modular symbols for squarefree N (v0.24.0).

The v0.23.0 tree could only see X_0(143): ``mtft.hecke`` hardwired
``LEVEL = 143``, the P^1 validity test named 11 and 13 explicitly, and the
triangle/cusp counts were baked in as the literals 56 and 4.  That made
*level-universality* — the program's largest standing risk, flagged after
Wave 7 — untestable in code.

This module rebuilds the same construction with N threaded through, so any
result proved at 143 can be replayed at other levels and either promoted to
a general theorem or exposed as a property of 143 alone.

Scope (enforced, not assumed).  Squarefree N with no elliptic points, i.e.
nu_2 = nu_3 = 0.  The Manin relations implemented here use free S- and
R-orbits; a level with elliptic points needs the extra torsion relations
(x + Sx = 0 with Sx = x forces 2x = 0), which are NOT implemented.  Such
levels raise :class:`UnsupportedLevelError` rather than returning a wrong
answer silently.  N = 21 is the smallest excluded case (two order-3 points).

Verification.  ``manin_model(143)`` reproduces ``mtft.hecke.model()`` field
by field; that equality is a shipped test, and it is what licenses reading
any level-generic result back onto the certified 143 work.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from functools import lru_cache
from math import gcd

__all__ = [
    "UnsupportedLevelError", "level_data", "is_supported", "check_supported",
    "supported_levels", "manin_model", "hecke_matrix", "cuspidal_hecke",
    "boundary_matrix", "cusp_labels", "genus",
]


class UnsupportedLevelError(ValueError):
    """Raised for levels outside the implemented Manin-relation scope."""


# ── level arithmetic ────────────────────────────────────────────────

def _prime_factors(N):
    fs, n, p = [], N, 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            fs.append((p, e))
        p += 1
    if n > 1:
        fs.append((n, 1))
    return fs


def _divisors(N):
    return sorted(d for d in range(1, N + 1) if N % d == 0)


def level_data(N: int) -> dict:
    """Index, elliptic-point counts, cusp count and genus of X_0(N).

    nu_2 and nu_3 are counted directly as solutions of x^2 + 1 = 0 and
    x^2 + x + 1 = 0 mod N, which sidesteps Kronecker-symbol conventions
    and is exact for every N.
    """
    if N < 1:
        raise UnsupportedLevelError(f"level must be positive, got {N}")
    fs = _prime_factors(N)
    squarefree = all(e == 1 for _, e in fs)
    index = 1
    for p, e in fs:
        index *= (p + 1) * p ** (e - 1)
    nu2 = sum(1 for x in range(N) if (x * x + 1) % N == 0)
    nu3 = sum(1 for x in range(N) if (x * x + x + 1) % N == 0)
    ncusps = len(_divisors(N)) if squarefree else None
    g = None
    if ncusps is not None:
        g12 = 12 + index - 3 * nu2 - 4 * nu3 - 6 * ncusps
        assert g12 % 12 == 0, "genus formula non-integral"
        g = g12 // 12
    return dict(N=N, factors=fs, squarefree=squarefree, index=index,
                nu2=nu2, nu3=nu3, ncusps=ncusps, genus=g)


def genus(N: int) -> int:
    return level_data(N)["genus"]


def is_supported(N: int) -> bool:
    try:
        check_supported(N)
        return True
    except UnsupportedLevelError:
        return False


def check_supported(N: int) -> dict:
    """Return level data, or raise with the specific reason for exclusion."""
    d = level_data(N)
    if not d["squarefree"]:
        raise UnsupportedLevelError(
            f"N = {N} is not squarefree (factors {d['factors']}); the cusp "
            "classification and boundary map here assume squarefree level")
    if d["nu2"] or d["nu3"]:
        raise UnsupportedLevelError(
            f"N = {N} has elliptic points (nu_2 = {d['nu2']}, nu_3 = "
            f"{d['nu3']}); the torsion Manin relations are not implemented")
    if d["genus"] == 0:
        raise UnsupportedLevelError(
            f"N = {N} has genus 0; there is no cuspidal homology to build")
    return d


def supported_levels(limit: int = 200) -> list:
    """Squarefree, elliptic-point-free, positive-genus levels up to *limit*."""
    return [N for N in range(1, limit + 1) if is_supported(N)]


# ── the model ───────────────────────────────────────────────────────

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


def merel(p: int):
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


@lru_cache(maxsize=8)
def manin_model(N: int) -> dict:
    """Exact P^1(Z/N) Manin model for squarefree, elliptic-free level N.

    Mirrors ``mtft.hecke.model()`` with every 143-specific literal replaced
    by its level-generic value: the validity test becomes gcd(c, d, N) = 1,
    the triangle count becomes len(tris), and the cusp count becomes the
    number of divisors of N.  Cusps additionally carry their canonical
    divisor label d = gcd(c, N).
    """
    check_supported(N)
    units = [u for u in range(1, N) if gcd(u, N) == 1]

    def valid(c, d):
        return gcd(gcd(c, d), N) == 1

    cache = {}

    def canon(c, d):
        key = (c % N, d % N)
        hit = cache.get(key)
        if hit is None:
            hit = min(((key[0] * u) % N, (key[1] * u) % N) for u in units)
            cache[key] = hit
        return hit

    P1 = sorted({canon(c, d) for c in range(N) for d in range(N)
                 if valid(c, d)})
    idx = {p: i for i, p in enumerate(P1)}

    def perm(f):
        return [idx[canon(*f(*p))] for p in P1]

    sS = perm(lambda c, d: (d % N, (-c) % N))
    sT = perm(lambda c, d: (c, (c + d) % N))
    sR = perm(lambda c, d: (d % N, (d - c) % N))
    iota = perm(lambda c, d: ((-c) % N, d))

    tris = _orbits(sR)
    ntri = len(tris)
    tri_of = {}
    for ti, t in enumerate(tris):
        for f in t:
            tri_of[f] = ti
    fans = sorted(_orbits(sT), key=len)
    ncusp = len(fans)
    cusp_of = {}
    for k, o in enumerate(fans):
        for f in o:
            cusp_of[f] = k
    # canonical divisor label of each cusp: d = gcd(c, N), constant on a
    # T-orbit because T fixes the first coordinate.
    cusp_label = []
    for o in fans:
        ds = {gcd(P1[f][0], N) for f in o}
        assert len(ds) == 1, "cusp label not constant on T-orbit"
        cusp_label.append(ds.pop())
    assert sorted(cusp_label) == _divisors(N), \
        f"cusp labels {sorted(cusp_label)} != divisors of {N}"

    Eorb = _orbits(sS)
    E = len(Eorb)
    erep, eid, esign = {}, {}, {}
    for k, o in enumerate(Eorb):
        assert len(o) == 2, "S-orbit of length != 2 (elliptic point?)"
        r = min(o)
        erep[k] = r
        for f in o:
            eid[f] = k
            esign[f] = 1 if f == r else -1

    D2 = [[0] * ntri for _ in range(E)]
    for ti, t in enumerate(tris):
        assert len(t) == 3, "R-orbit of length != 3 (elliptic point?)"
        for f in t:
            D2[eid[f]][ti] += esign[f]
    DEL = [[0] * ncusp for _ in range(E)]
    for k in range(E):
        x = erep[k]
        DEL[k][cusp_of[x]] += 1
        DEL[k][cusp_of[sS[x]]] -= 1

    _, cols = _rref([[D2[e][t] for t in range(ntri)] for e in range(E)])
    D2r = [[D2[e][t] for t in cols] for e in range(E)]
    _, pivE = _rref([[D2r[e][j] for e in range(E)] for j in range(len(cols))])
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

    DB = [[Fr(DEL[free[j]][c]) for j in range(nq)] for c in range(ncusp)]
    Mr, pivK = _rref(DB)
    freeK = [j for j in range(nq) if j not in set(pivK)]
    K = []
    for fj in freeK:
        v = [Fr(0)] * nq
        v[fj] = Fr(1)
        for r_, c_ in enumerate(pivK):
            v[c_] = -Mr[r_][fj]
        K.append(v)

    g = level_data(N)["genus"]
    assert len(K) == 2 * g, f"cuspidal rank {len(K)} != 2g = {2 * g}"
    assert nq == 2 * g + ncusp - 1, \
        f"relative rank {nq} != 2g + #cusps - 1 = {2 * g + ncusp - 1}"

    return dict(N=N, P1=P1, idx=idx, canon=canon, sS=sS, sT=sT, sR=sR,
                iota=iota, tris=tris, tri_of=tri_of, fans=fans,
                cusp_of=cusp_of, cusp_label=cusp_label, ncusp=ncusp,
                E=E, erep=erep, eid=eid, esign=esign, D2=D2, DEL=DEL,
                cols=cols, D2r=D2r, free=free, Binv=Binv, nq=nq,
                K=K, freeK=freeK, genus=g)


def cusp_labels(N: int) -> list:
    """Divisor label d | N of each cusp, in model index order."""
    return list(manin_model(N)["cusp_label"])


# ── operators ───────────────────────────────────────────────────────

def _qcoords(m, v):
    E, Binv, ncols = m["E"], m["Binv"], len(m["cols"])
    x = [sum(Binv[i][k] * v[k] for k in range(E) if v[k]) for i in range(E)]
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
                continue
            y = idx[canon(u, v)]
            T[eid[y]][k] += esign[y]
    return T


def _quotient(m, T):
    nq, free, E = m["nq"], m["free"], m["E"]
    cs = [_qcoords(m, [Fr(T[e][fj]) for e in range(E)]) for fj in free]
    return [[cs[j][i] for j in range(nq)] for i in range(nq)]


def _restrict_cuspidal(m, Q):
    K, freeK, nq = m["K"], m["freeK"], m["nq"]
    out = []
    for k in K:
        img = [sum(Q[i][j] * k[j] for j in range(nq)) for i in range(nq)]
        out.append([img[fj] for fj in freeK])
    d = len(K)
    return [[out[j][i] for j in range(d)] for i in range(d)]


@lru_cache(maxsize=128)
def hecke_matrix(N: int, p: int):
    """T_p (p prime to N) or U_p (p | N) on the relative model.  EXACT."""
    m = manin_model(N)
    return _quotient(m, _raw(m, merel(p)))


@lru_cache(maxsize=128)
def cuspidal_hecke(N: int, p: int):
    """T_p on the 2g-dimensional cuspidal homology H_1(X_0(N), Z).  EXACT.

    Returns tuple-of-tuples, matching ``mtft.hecke.cuspidal_hecke``.
    """
    m = manin_model(N)
    return tuple(tuple(r) for r in _restrict_cuspidal(m, hecke_matrix(N, p)))


def boundary_matrix(N: int):
    """Boundary map to the degree-zero cusp lattice, in the relative model.

    Rows are the basis (C_{d_1} - C_{d_0}, ..., C_{d_k} - C_{d_0}) of the
    degree-zero cusp divisors, where d_0 is the first cusp in model order.
    """
    m = manin_model(N)
    DEL, free, ncusp = m["DEL"], m["free"], m["ncusp"]
    D4 = [[DEL[free[j]][c] for j in range(m["nq"])] for c in range(ncusp)]
    return [row[:] for row in D4[1:]]
