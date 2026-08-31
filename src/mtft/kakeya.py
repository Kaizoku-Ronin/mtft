"""mtft.kakeya — Arf parity, direction sets, and finite Kakeya geometry.

Three things, one theorem among them.

1.  THE ARF DIRECTION THEOREM (certified here, exact).
    Let A be the G-invariant theta-characteristic affine space of X_0(143):
    128 characteristics, an affine F_2-space of dimension 7, carrying the Arf
    parity q with 96 even and 32 odd.  Let V = Delta(A) be its direction space
    and R the radical of the bilinear form associated with q.

    The observed 96/32 split *forces* q|_R = 0.  (For r in R the bilinear form
    vanishes, so q(a+r) = q(a) + q(r); if q|_R were a nonzero functional it
    would be surjective and the odd count would be |A|/2 = 64, not 32.)
    Hence q descends to V/R and Theta_odd is a single coset of R, giving

        Delta(Theta_odd)  = R   ->  2^5 - 1  =  31 nonzero directions
        Delta(Theta_even) = V   ->  2^7 - 1  = 127 nonzero directions

    Odd spin states are directionally *confined to the radical*; even states
    are directionally *saturated*.  The premise is not assumed: it is derived
    from the counts, so the theorem stands on the census numbers alone.

    This is Kakeya-adjacent, not Kakeya: it concerns difference directions,
    not the existence of a full line in every direction.  Do not conflate the
    two -- see ``arf_direction_theorem`` return field ``is_kakeya_theorem``.

2.  FINITE KAKEYA over a prime field, for the genuine article.  Besicovitch
    sets in F_p^2 built as unions of tangent lines to a parabola, checked
    against Dvir's bound  |K| >= binom(p+1, 2) = p(p+1)/2.

3.  THE CRT DIRECTION BRIDGE.  Z/143 is a product ring with zero divisors, so
    Dvir does not apply to it directly; the honest statement splits first:

        P^1(Z/143)  ~  P^1(F_11) x P^1(F_13),    12 * 14 = 168,

    and 168 is exactly the index [PSL_2(Z) : Gamma_0(143)] and the Manin
    symbol count.  The direction set of the modular symbol model factors as a
    product of two finite projective lines.

OPEN (not shipped): the Atkin-Lehner permutation of P^1(Z/N).  W_Q has
determinant Q, not 1, so the naive row action mod N is not well defined on
P^1 -- a correct implementation needs the coset-level action.  Deliberately
omitted rather than shipped unverified.  Target for the next wave.
"""

from __future__ import annotations

from math import gcd, comb
from itertools import product

import numpy as np

__all__ = [
    "direction_set", "affine_frame", "radical_of_parity",
    "arf_direction_theorem",
    "dvir_bound", "besicovitch_set", "is_kakeya_set", "kakeya_report",
    "p1_points", "crt_direction_bridge",
]


# ---------------------------------------------------------------- F_2 tools

def _f2(rows) -> np.ndarray:
    return np.array(rows, dtype=np.int64) % 2


def _rank2(M) -> int:
    M = _f2(M).copy()
    if M.size == 0:
        return 0
    m, n = M.shape
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, m) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(m):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    return r


def direction_set(rows) -> set:
    """Nonzero pairwise differences of F_2 vectors.

    Over F_2 subtraction is addition, so this is the set of sums x + y for
    distinct x, y in ``rows``.  Returned as a set of byte-tuples.
    """
    R = _f2(rows)
    out = set()
    n = len(R)
    for i in range(n):
        d = (R[i] ^ R) % 2          # all differences against row i at once
        for j in range(i + 1, n):
            out.add(tuple(int(t) for t in d[j]))
    out.discard(tuple([0] * R.shape[1]))
    return out


def affine_frame(rows):
    """Base point and a basis of the direction space of an affine F_2 set.

    Returns ``(t0, B, dim)`` with ``B`` a list of basis vectors spanning
    Delta(rows) as a linear space.
    """
    R = _f2(rows)
    t0 = R[0].copy()
    diffs = (R ^ t0) % 2
    basis, cur = [], []
    for v in diffs:
        trial = cur + [v]
        if _rank2(trial) > len(cur):
            cur = trial
            basis.append(v.copy())
    return t0, basis, len(basis)


def _coords_in_basis(v, basis):
    """Express v in terms of ``basis`` over F_2 (assumes v is in the span)."""
    B = _f2(basis).T                      # columns = basis vectors
    M = np.hstack([B, _f2(v).reshape(-1, 1)])
    m, n = M.shape
    ncols = n - 1
    piv = {}
    r = 0
    for c in range(ncols):
        p = next((i for i in range(r, m) if M[i, c]), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        for i in range(m):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        piv[c] = r
        r += 1
    for i in range(r, m):
        if M[i, ncols]:
            raise ValueError("vector not in span")
    x = np.zeros(ncols, dtype=np.int64)
    for c, ri in piv.items():
        x[c] = M[ri, ncols]
    return x


def radical_of_parity(rows, parity):
    """Radical of the bilinear form attached to an Arf parity function.

    ``rows`` is an affine F_2 set, ``parity`` the 0/1 value of q on each row.
    Returns ``(radical_vectors, dim_V, dim_R)`` with radical vectors given in
    the ambient coordinates.
    """
    R = _f2(rows)
    par = np.array(parity, dtype=np.int64) % 2
    t0, basis, dimV = affine_frame(R)
    index = {tuple(int(t) for t in row): i for i, row in enumerate(R)}

    def q(vec):
        """Parity at base point + vec."""
        key = tuple(int(t) for t in ((t0 ^ _f2(vec)) % 2))
        return int(par[index[key]])

    q0 = q(np.zeros_like(t0))
    # bilinear form on the basis
    G = np.zeros((dimV, dimV), dtype=np.int64)
    for i in range(dimV):
        for j in range(dimV):
            s = (basis[i] ^ basis[j]) % 2
            G[i, j] = (q(s) + q(basis[i]) + q(basis[j]) + q0) % 2
    # radical = kernel of G
    M = G.copy()
    m, n = M.shape
    piv = {}
    r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if M[i, c]), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        for i in range(m):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        piv[c] = r
        r += 1
    free = [c for c in range(n) if c not in piv]
    rad = []
    for fc in free:
        coef = np.zeros(n, dtype=np.int64)
        coef[fc] = 1
        for c, ri in piv.items():
            coef[c] = M[ri, fc]
        vec = np.zeros_like(t0)
        for k, ck in enumerate(coef):
            if ck:
                vec = (vec ^ basis[k]) % 2
        rad.append(vec)
    return rad, dimV, len(rad)


def arf_direction_theorem(rows=None, parity=None):
    """Certify Delta(Theta_odd) = R and Delta(Theta_even) = V.

    With no arguments, loads the X_0(143) G-invariant theta characteristics
    from :mod:`mtft.thetachar`.  Every assertion in the docstring is checked
    and reported; nothing is taken on faith.
    """
    if rows is None:
        from mtft import thetachar as TC
        S = TC.x0143_periods_frame()
        rows, parity = S.invariant_characteristics()
    R = _f2(rows)
    par = np.array(parity, dtype=np.int64) % 2

    n_even = int((par == 0).sum())
    n_odd = int((par == 1).sum())
    rad, dimV, dimR = radical_of_parity(R, par)

    odd_rows = R[par == 1]
    even_rows = R[par == 0]
    D_odd = direction_set(odd_rows)
    D_even = direction_set(even_rows)

    rad_span = set()
    if rad:
        for coeffs in product((0, 1), repeat=len(rad)):
            v = np.zeros(R.shape[1], dtype=np.int64)
            for c, rv in zip(coeffs, rad):
                if c:
                    v = (v ^ rv) % 2
            rad_span.add(tuple(int(t) for t in v))
    rad_span.discard(tuple([0] * R.shape[1]))

    t0, basis, _ = affine_frame(R)
    V_span = set()
    for coeffs in product((0, 1), repeat=len(basis)):
        v = np.zeros(R.shape[1], dtype=np.int64)
        for c, bv in zip(coeffs, basis):
            if c:
                v = (v ^ bv) % 2
        V_span.add(tuple(int(t) for t in v))
    V_span.discard(tuple([0] * R.shape[1]))

    # the forcing argument: q|_R = 0 is implied by n_odd == |R|
    q_restricted_trivial = (n_odd == 2 ** dimR)
    odd_is_single_coset = (n_odd == 2 ** dimR)

    return {
        "n_characteristics": int(len(R)),
        "n_even": n_even,
        "n_odd": n_odd,
        "dim_V": dimV,
        "dim_R": dimR,
        "q_restricted_to_radical_is_zero": bool(q_restricted_trivial),
        "odd_locus_is_single_R_coset": bool(odd_is_single_coset),
        "n_directions_odd": len(D_odd),
        "n_directions_even": len(D_even),
        "expected_odd": 2 ** dimR - 1,
        "expected_even": 2 ** dimV - 1,
        "delta_odd_equals_radical": D_odd == rad_span,
        "delta_even_equals_V": D_even == V_span,
        "odd_confined": len(D_odd) < len(D_even),
        "is_kakeya_theorem": False,
        "note": ("difference-direction statement, not a Kakeya line-existence "
                 "statement; see module docstring"),
    }


# ------------------------------------------------------- finite Kakeya, F_p

def dvir_bound(p: int) -> int:
    """Dvir's lower bound for a Kakeya set in F_p^2:  binom(p+1, 2)."""
    return comb(p + 1, 2)


def besicovitch_set(p: int) -> set:
    """A near-extremal Kakeya set in F_p^2 (union of tangents to a parabola).

    The lines L_m = {(t, m t - m^2)} for m in F_p, together with one vertical
    line for the direction at infinity.  Two tangents meet at (m+m', m m'), so
    the union is {(s, u) : s^2 - 4u is a square}, of size p(p+1)/2, which meets
    Dvir's bound up to the single extra vertical line.
    """
    K = set()
    for m in range(p):
        for t in range(p):
            K.add((t % p, (m * t - m * m) % p))
    for y in range(p):                      # vertical direction
        K.add((0, y))
    return K


def is_kakeya_set(K: set, p: int) -> bool:
    """True iff K contains a full line in every direction of P^1(F_p)."""
    for m in range(p):                      # slope m
        if not any(all(((t, (m * t + b) % p) in K) for t in range(p))
                   for b in range(p)):
            return False
    return any(all(((x0, y) in K) for y in range(p)) for x0 in range(p))


def kakeya_report(p: int) -> dict:
    K = besicovitch_set(p)
    return {
        "p": p,
        "size": len(K),
        "dvir_bound": dvir_bound(p),
        "ambient": p * p,
        "is_kakeya": is_kakeya_set(K, p),
        "meets_bound": len(K) >= dvir_bound(p),
        "excess_over_bound": len(K) - dvir_bound(p),
        "fraction_of_plane": len(K) / (p * p),
    }


# ------------------------------------------------------- the CRT bridge

def p1_points(N: int):
    """Canonical representatives of P^1(Z/N)."""
    units = [u for u in range(N) if gcd(u, N) == 1]
    seen, out = set(), []
    for c in range(N):
        for d in range(N):
            if gcd(gcd(c, d), N) != 1:
                continue
            if (c, d) in seen:
                continue
            orbit = {((c * u) % N, (d * u) % N) for u in units}
            seen |= orbit
            out.append(min(orbit))
    return sorted(out)


def crt_direction_bridge(p: int, q: int) -> dict:
    """Verify P^1(Z/pq) ~ P^1(F_p) x P^1(F_q) as sets, explicitly."""
    N = p * q
    big = p1_points(N)
    small_p = p1_points(p)
    small_q = p1_points(q)

    def reduce_to(pt, m):
        c, d = pt[0] % m, pt[1] % m
        units = [u for u in range(m) if gcd(u, m) == 1]
        return min(((c * u) % m, (d * u) % m) for u in units)

    image = {(reduce_to(pt, p), reduce_to(pt, q)) for pt in big}
    return {
        "N": N,
        "size_P1_N": len(big),
        "size_P1_p": len(small_p),
        "size_P1_q": len(small_q),
        "product": len(small_p) * len(small_q),
        "bijective": len(image) == len(big) == len(small_p) * len(small_q),
        "index_Gamma0": (p + 1) * (q + 1),
        "equals_index_Gamma0": len(big) == (p + 1) * (q + 1),
    }
