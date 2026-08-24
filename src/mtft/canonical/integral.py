"""mtft.canonical.integral — integral models, saturation, and the 2026-08-24 arc.

Everything certified in certificates v1-v9 of the integral-model arc,
recomputed at call time from the frozen q-expansions (house policy).
Deterministic: saturation pivots are unimodular away from p and primes
are processed in a fixed order, so every derived lattice is reproducible.

The three integral models of V(I_2) mod 2 (the Integral Model Gate story):
    packaged s2 model            ->  7 points  (a plane; s2 lattice v_2 = 15)
    adapted coords + sat. ideal  ->  3 points  (cusps collapse; v_2 = 25)
    fully saturated model        ->  4 points  (the curve; = the four cusps)
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import gcd

import numpy as np

from .. import integral_lattice as lat
from . import (COORDINATE_LABELS, GENUS, MONOMIALS, adapted_qexpansions,
               ideal_basis, ideal_basis_adapted)

__all__ = [
    "SATURATION_PRIMES", "LEDGER",
    "adapted_matrix", "saturated_qexpansions", "sector_columns",
    "count_points_modp", "points_modp", "cusp_reductions", "ci_a_codifferent",
    "quadratic_saturation_obstruction", "al_splitting", "al_denominator",
]

SATURATION_PRIMES = (2, 3, 5, 7, 13, 19, 103, 5560463)

#: Certified instance values (v1-v9); gates re-derive every entry.
LEDGER = {
    "full_saturation_steps": {2: 25, 3: 8, 5: 1, 7: 1, 13: 1, 19: 1,
                              103: 1, 5560463: 1},
    "counts_mod2": {"packaged_s2": 7, "adapted_mixed": 3, "saturated": 4},
    "counts_mod3": {"saturated": 4},
    "ci_a_packaged_a": -2439613813,
    "ci_a_codifferent_a": -637,
    "f2_orbit_saturation_steps": {2: 6, 3: 2, 19: 1, 103: 1},
    "Q_full_sector": [1] * 9 + [13],
    "Q_f2_codifferent": [1] * 9 + [637],
    "product_index_full_over_f2": 49,
    "al_four_sector_snf": [1] * 6 + [2] * 6 + [52],
    "al_W11_split_snf": [1] * 7 + [2] * 6,
    "al_W13_in_W11plus": [1] * 6 + [2],
    "al_W13_in_W11minus": [1] * 5 + [26],
    "al_denominators": {11: 1, 13: 13},
}

_SLICES = {"(+,+)": [0], "(+,-)": list(range(1, 7)),
           "(-,+)": list(range(7, 12)), "(-,-)": [12]}


def adapted_matrix():
    """The 141 x 13 adapted q-coefficient matrix as an object array."""
    A = np.array([[int(v) for v in row]
                  for row in np.array(adapted_qexpansions(), dtype=object)],
                 dtype=object)
    return A if A.shape[0] > A.shape[1] else A.T


@lru_cache(maxsize=None)
def _saturated():
    B, steps = lat.saturate(adapted_matrix(), SATURATION_PRIMES)
    return B, steps


def saturated_qexpansions():
    """(fully saturated 141 x 13 basis, {p: steps})."""
    B, steps = _saturated()
    return B.copy(), dict(steps)


def sector_columns(sector, saturated=True):
    """The (individually saturated) lattice of one AL sector."""
    A = adapted_matrix()[:, _SLICES[sector]]
    if not saturated:
        return A
    B, _ = lat.saturate(A, SATURATION_PRIMES)
    return B


def _pairs(k):
    return [(i, j) for i in range(k) for j in range(i, k)]


def _sym2_kernel_modp(B, p):
    Bp = np.array([[int(x) % p for x in row] for row in B], dtype=np.int64)
    cols = []
    for i, j in _pairs(B.shape[1]):
        cols.append(list(np.convolve(Bp[:, i], Bp[:, j])[:141] % p))
    _, ker = lat.kernel_modp(np.array(cols, dtype=np.int64).T, p)
    return ker


def count_points_modp(p, model="saturated"):
    """#V(F_p) of the chosen integral model.  p in {2, 3} (larger p is a
    gprun-sized job; refuse rather than stall)."""
    if p not in (2, 3):
        raise ValueError("p**13 evaluations: size this as a gprun job")
    if model == "saturated":
        B, _ = _saturated()
        ker = _sym2_kernel_modp(B, p)
        Q = np.array(ker, dtype=np.int64).T % p          # 91 x 55
        pair_list = _pairs(13)
    elif model == "packaged_s2":
        Q = np.array(ideal_basis(), dtype=np.int64)
        cols = Q.T.copy()
        for k in range(cols.shape[0]):
            g = 0
            for v in cols[k]:
                g = gcd(g, abs(int(v)))
            cols[k] //= g
        Q = (cols.T % p)
        pair_list = list(MONOMIALS)
    elif model == "adapted_mixed":
        Q = np.array(ideal_basis_adapted(), dtype=np.int64) % p
        pair_list = list(MONOMIALS)
    else:
        raise ValueError(model)
    N = p ** GENUS
    powers = np.array([p ** t for t in range(GENUS)], dtype=np.int64)
    X = ((np.arange(N, dtype=np.int64)[:, None] // powers[None, :]) % p)
    Mon = np.empty((N, len(pair_list)), dtype=np.float64)
    for m, (i, j) in enumerate(pair_list):
        Mon[:, m] = (X[:, i] * X[:, j]) % p
    good = (np.rint(Mon @ Q.astype(np.float64)).astype(np.int64) % p == 0
            ).all(axis=1)
    return (int(good.sum()) - 1) // (p - 1)


def points_modp(p):
    """Projective representatives of the SATURATED model's F_p points."""
    B, _ = _saturated()
    ker = _sym2_kernel_modp(B, p)
    Q = np.array(ker, dtype=np.int64).T % p
    N = p ** GENUS
    powers = np.array([p ** t for t in range(GENUS)], dtype=np.int64)
    X = ((np.arange(N, dtype=np.int64)[:, None] // powers[None, :]) % p)
    Mon = np.empty((N, Q.shape[0]), dtype=np.float64)
    for m, (i, j) in enumerate(_pairs(GENUS)):
        Mon[:, m] = (X[:, i] * X[:, j]) % p
    good = (np.rint(Mon @ Q.astype(np.float64)).astype(np.int64) % p == 0
            ).all(axis=1)
    reps, seen = [], set()
    for row in X[good]:
        t = tuple(int(x) for x in row)
        if not any(t):
            continue
        if all(tuple((s * np.array(t)) % p) not in seen for s in range(1, p)):
            seen.add(t)
            reps.append(t)
    return sorted(reps)


def cusp_reductions(p):
    """{cusp d: F_p point of the SATURATED model}, plus a bijection flag
    against the counted points (certificates v2/v3)."""
    A = adapted_matrix()
    B, _ = _saturated()
    a1 = [int(x) for x in A[1]]
    sig = [(1 if s[1] == "+" else -1, 1 if s[3] == "+" else -1)
           for _, s in COORDINATE_LABELS]
    H = lat.hnf(B)
    out = {}
    for d, eps in ((1, lambda i: 1), (11, lambda i: sig[i][0]),
                   (13, lambda i: sig[i][1]),
                   (143, lambda i: sig[i][0] * sig[i][1])):
        x = [eps(i) * a1[i] for i in range(13)]
        # coordinates of the cusp in the saturated basis: solve over Q
        y = _rational_coords(B, x)
        den = 1
        for t in y:
            den = den * t.denominator // gcd(den, t.denominator)
        yi = [int(t * den) for t in y]
        g = 0
        for t in yi:
            g = gcd(g, abs(t))
        yi = [t // g for t in yi]
        yp = [t % p for t in yi]
        fn = next(i for i, v in enumerate(yp) if v)
        inv = pow(yp[fn], p - 2, p)
        out[d] = tuple((v * inv) % p for v in yp)
    return out


def _rational_coords(B, coeff_vector):
    """Solve X with A = B X exactly (normal equations over Q), then return
    y = (X^T)^{-1} x: the cusp's homogeneous coordinates in the basis B
    itself — the same basis `points_modp` reduces, so the two agree."""
    A = adapted_matrix()
    n = 13
    G = [[sum(int(B[k][i]) * int(B[k][j]) for k in range(B.shape[0]))
          for j in range(n)] for i in range(n)]
    R = [[sum(int(B[k][i]) * int(A[k][j]) for k in range(B.shape[0]))
          for j in range(n)] for i in range(n)]
    X = _solve_square(G, R)                                  # A = B X
    # point transforms contravariantly: y = (X^T)^{-1} x
    M = [[X[j][i] for j in range(n)] for i in range(n)]      # X^T
    aug = [row[:] + [Fraction(coeff_vector[i])] for i, row in enumerate(M)]
    for c in range(n):
        pr = next(i for i in range(c, n) if aug[i][c] != 0)
        aug[c], aug[pr] = aug[pr], aug[c]
        pivinv = 1 / aug[c][c]
        aug[c] = [t * pivinv for t in aug[c]]
        for i in range(n):
            if i != c and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [t - f * s for t, s in zip(aug[i], aug[c])]
    return [aug[i][n] for i in range(n)]



def _solve_square(G, R):
    """Solve G X = R exactly for square integer G (invertible), Fractions."""
    from fractions import Fraction as F
    n = len(G)
    aug = [[F(G[i][j]) for j in range(n)] + [F(R[i][j]) for j in range(n)]
           for i in range(n)]
    for c in range(n):
        pr = next(i for i in range(c, n) if aug[i][c] != 0)
        aug[c], aug[pr] = aug[pr], aug[c]
        inv = 1 / aug[c][c]
        aug[c] = [t * inv for t in aug[c]]
        for i in range(n):
            if i != c and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [t - f * s for t, s in zip(aug[i], aug[c])]
    return [[aug[i][n + j] for j in range(n)] for i in range(n)]

def ci_a_codifferent():
    """The codifferent-normalized CI-A relation (certificate v7).

    Returns (a, quadric_dict) with a = -637 = -13 * C_Eis(f2); the
    quadric is expressed on {f1_primitive} u {saturated f2 orbit}.
    """
    A = adapted_matrix()
    f1 = A[:, 0]
    g = 0
    for v in f1:
        g = gcd(g, abs(int(v)))
    f1p = np.array([int(v) // g for v in f1], dtype=object)
    F2, _ = lat.saturate(A[:, [8, 9, 10, 11]], (2, 3, 19, 103))
    basis = [f1p] + [F2[:, k] for k in range(4)]
    pairs = _pairs(5)
    cols = [[int(x) for x in np.convolve(basis[i], basis[j])[:141]]
            for (i, j) in pairs]
    ker = lat.rational_kernel(np.array(cols, dtype=object).T)
    assert len(ker) == 1, "CI-A uniqueness"
    v = ker[0]
    den = 1
    for e in v:
        den = den * e.denominator // gcd(den, e.denominator)
    vi = [int(e * den) for e in v]
    g = 0
    for e in vi:
        g = gcd(g, abs(e))
    vi = [e // g for e in vi]
    if vi[0] > 0:
        vi = [-e for e in vi]
    return vi[0], dict(zip(pairs, vi))


def quadratic_saturation_obstruction(basis, primes):
    """Q(L) = Sat(im Sym^2 basis)/im, as Smith invariants of the quotient.

    `primes` must cover the primes of the index (instance: (7, 13) for the
    (-,+) sector and the f2 orbit; certificates v8).  Route: column HNF of
    the product lattice, saturation of that small basis, then Smith of the
    10 x 10 coordinate matrix — never Smith of the raw 141-row matrix.
    """
    k = basis.shape[1]
    cols = [[int(x) for x in np.convolve(basis[:, i], basis[:, j])[:141]]
            for (i, j) in _pairs(k)]
    M = np.array(cols, dtype=object).T
    H = lat.hnf(M)
    Hs, _ = lat.saturate(H, primes)
    inv = lat.quotient_invariants(Hs, H)
    r = H.shape[1]
    return [1] * (r - len(inv)) + inv


def al_splitting():
    """The Atkin-Lehner splitting obstruction filtration (certificate v9)."""
    L, _ = _saturated()
    sec = {s: sector_columns(s) for s in _SLICES}
    D4 = np.concatenate([sec["(+,+)"], sec["(+,-)"],
                         sec["(-,+)"], sec["(-,-)"]], axis=1)
    W11p, _ = lat.saturate(adapted_matrix()[:, list(range(0, 7))],
                           SATURATION_PRIMES)
    W11m, _ = lat.saturate(adapted_matrix()[:, list(range(7, 13))],
                           SATURATION_PRIMES)
    return {
        "four_sector": lat.quotient_invariants(L, D4),
        "W11_split": lat.quotient_invariants(
            L, np.concatenate([W11p, W11m], axis=1)),
        "W13_in_W11plus": lat.quotient_invariants(
            W11p, np.concatenate([sec["(+,+)"], sec["(+,-)"]], axis=1)),
        "W13_in_W11minus": lat.quotient_invariants(
            W11m, np.concatenate([sec["(-,+)"], sec["(-,-)"]], axis=1)),
    }


def al_denominator(p):
    """Integral denominator of W_p on the saturated H^0(K) lattice.

    W11 -> 1 (an honest integral involution: parity lemma applies);
    W13 -> 13 (oldspace scaling g(q) +- 13 g(q^13); certificate v9).
    """
    signs = {11: [1] * 7 + [-1] * 6,
             13: [1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1]}[p]
    L, _ = _saturated()
    _, den = lat.operator_matrix(L, adapted_matrix(), signs)
    return den
