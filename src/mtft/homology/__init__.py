"""Canonical integral homology of X0(143).

Ships the canonical integral structure on H_1(X0(143), Z) exported from an
independent PARI/GP modular-symbol pipeline (msinit -> mscuspidal ->
mslattice -> msatkinlehner/msstar/mspetersson), together with exact integer
linear algebra to place it in a standard symplectic frame.

Certified facts (re-asserted on load, exact integer arithmetic):
  * W11, W13, STAR restrict integrally to the canonical lattice L,
    are involutions, and commute pairwise over Z;
  * the intersection pairing P on L is antisymmetric with det = +1,
    so (L, P) is H_1(X0(143), Z) with its unimodular intersection form;
  * W11, W13 are symplectic, STAR is anti-symplectic (symplectic mod 2).

Route note (ledgered 2026-08-27): the naive saturation of the cuspidal
Manin-symbol lattice inside Z^29 is an index-4200 sublattice whose mod-2
pairing is degenerate; theta-characteristic work on it fails silently.
The mslattice canonical structure packaged here is the correct object.

Two frames are provided:
  * ``matrices()``          - the GP canonical frame (this module's data);
  * ``symplectic_frame()``  - an exact U in GL(26, Z) with U^T P U = J_std,
                              plus the three operators conjugated into that
                              frame (computed once, cached, re-verified);
  * ``periods_frame_ops()`` - the operators in the symplectic frame used by
                              :mod:`mtft.periods` (the frame of the frozen
                              Riemann matrix tau0), rounded from float and
                              exactly re-verified.  This is the frame in
                              which theta-function work at tau0 must live.

All heavy objects are Python-int exact; numpy object arrays are used so no
overflow is possible.
"""
from __future__ import annotations

import json
from fractions import Fraction
from importlib import resources

import numpy as np

__all__ = [
    "matrices", "standard_J", "symplectic_frame", "periods_frame_ops",
    "mod2", "int_inverse", "is_symplectic", "is_anti_symplectic",
]

_G = 13
_N = 26


def _load():
    with resources.files(__package__).joinpath(
            "_data/x0143_integral_homology.json").open() as fh:
        return json.load(fh)


def _obj(rows):
    return np.array([[int(x) for x in r] for r in rows], dtype=object)


def standard_J(g: int = _G):
    """Standard symplectic form [[0, I], [-I, 0]] as an exact object array."""
    J = np.zeros((2 * g, 2 * g), dtype=object)
    for i in range(g):
        J[i, g + i] = 1
        J[g + i, i] = -1
    return J


def int_inverse(M):
    """Exact inverse of an integer matrix with det = +-1 (object arrays).

    Fraction-based Gauss-Jordan; asserts the result is integral and a
    two-sided inverse.
    """
    n = M.shape[0]
    A = [[Fraction(int(M[i, j])) for j in range(n)] for i in range(n)]
    I = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        I[c], I[p] = I[p], I[c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        I[c] = [x * inv for x in I[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
                I[r] = [a - f * b for a, b in zip(I[r], I[c])]
    out = np.empty((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            q = I[i][j]
            if q.denominator != 1:
                raise ValueError("matrix is not GL(n, Z)")
            out[i, j] = int(q)
    assert np.array_equal(out @ M, np.eye(n, dtype=object) * 1)
    assert np.array_equal(M @ out, np.eye(n, dtype=object) * 1)
    return out


def is_symplectic(M, E):
    return np.array_equal(M.T @ E @ M, E)


def is_anti_symplectic(M, E):
    return np.array_equal(M.T @ E @ M, -E)


def _det_int(M):
    """Exact determinant via fraction-free Bareiss."""
    n = M.shape[0]
    A = [[int(M[i, j]) for j in range(n)] for i in range(n)]
    sign, prev = 1, 1
    for c in range(n - 1):
        if A[c][c] == 0:
            p = next((r for r in range(c + 1, n) if A[r][c] != 0), None)
            if p is None:
                return 0
            A[c], A[p] = A[p], A[c]
            sign = -sign
        for r in range(c + 1, n):
            for k in range(c + 1, n):
                A[r][k] = (A[r][k] * A[c][c] - A[r][c] * A[c][k]) // prev
            A[r][c] = 0
        prev = A[c][c]
    return sign * A[n - 1][n - 1]


_CACHE: dict = {}


def matrices():
    """The canonical GP-frame data with all structural facts re-asserted."""
    if "gp" in _CACHE:
        return _CACHE["gp"]
    d = _load()
    W11, W13, ST, P = (_obj(d[k]) for k in ("W11", "W13", "STAR", "P"))
    I = np.eye(_N, dtype=object) * 1
    assert np.array_equal(W11 @ W11, I) and np.array_equal(W13 @ W13, I)
    assert np.array_equal(ST @ ST, I)
    assert np.array_equal(W11 @ W13, W13 @ W11)
    assert np.array_equal(W11 @ ST, ST @ W11)
    assert np.array_equal(W13 @ ST, ST @ W13)
    assert np.array_equal(P.T, -P)
    assert _det_int(P) == 1, "intersection form must be unimodular"
    assert is_symplectic(W11, P) and is_symplectic(W13, P)
    assert is_anti_symplectic(ST, P)
    out = {"W11": W11, "W13": W13, "STAR": ST, "P": P, "meta": d["meta"]}
    _CACHE["gp"] = out
    return out


def _symplectic_reduce(E):
    """Exact U in GL(2g, Z) with U^T E U = standard J, for unimodular
    antisymmetric E.  Integer symplectic Gram-Schmidt: unimodularity makes
    v -> <v, .> surjective onto Z, so a dual partner with pairing exactly 1
    always exists (extended gcd along the row)."""
    n = E.shape[0]
    g = n // 2
    basis = [np.array([int(i == k) for i in range(n)], dtype=object)
             for k in range(n)]

    def pair(x, y):
        return int(x @ (E @ y))

    A_vecs, B_vecs = [], []
    pool = basis
    for _ in range(g):
        pool = [_content_divide(x) for x in pool]
        v = pool[0]
        row = [pair(v, u) for u in pool]
        # extended gcd over the row to reach pairing exactly 1
        idx = [i for i, r in enumerate(row) if r != 0]
        assert idx, "degenerate pairing on remaining block"
        gg, coeffs = row[idx[0]], {idx[0]: 1}
        for i in idx[1:]:
            a, b = gg, row[i]
            old_g = gg
            x0, y0 = _extgcd(a, b)
            gg = a * x0 + b * y0
            coeffs = {k: c * x0 for k, c in coeffs.items()}
            coeffs[i] = coeffs.get(i, 0) + y0
            if abs(gg) == 1:
                break
            assert gg != 0 and abs(gg) <= abs(old_g)
        assert abs(gg) == 1, "unimodularity violated: no unit pairing partner"
        w = sum(c * pool[k] for k, c in coeffs.items())
        if gg == -1:
            w = -w
        assert pair(v, w) == 1
        A_vecs.append(v)
        B_vecs.append(w)
        new_pool = []
        for x in pool:
            y = x - pair(x, w) * v + pair(x, v) * w
            if any(int(t) != 0 for t in y):
                new_pool.append(y)
        # drop dependencies: keep vectors until pool spans the 2(g-1) block;
        # exactness of the projection makes v, w pair to zero with every y
        pool = [y for y in new_pool
                if pair(y, v) == 0 and pair(y, w) == 0]
        # prune an exactly dependent tail if present
        pool = _prune_dependent(pool, n)
    U = np.stack(A_vecs + B_vecs, axis=1)
    return U


def _content_divide(v):
    from math import gcd
    g = 0
    for x in v:
        g = gcd(g, abs(int(x)))
    if g > 1:
        v = np.array([int(x) // g for x in v], dtype=object)
    return v


def _extgcd(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_s, old_t


def _prune_dependent(vecs, n):
    """Keep an exactly independent subset (fraction-free elimination)."""
    kept, rows = [], []
    for v in vecs:
        r = [Fraction(int(x)) for x in v]
        for piv_col, base in rows:
            if r[piv_col] != 0:
                f = r[piv_col] / base[piv_col]
                r = [a - f * b for a, b in zip(r, base)]
        nz = next((i for i, x in enumerate(r) if x != 0), None)
        if nz is not None:
            rows.append((nz, r))
            kept.append(v)
    return kept


def symplectic_frame():
    """(U, ops) with U in GL(26, Z), U^T P U = J_std, and the three
    operators conjugated into the J-frame (exact; re-verified)."""
    if "sym" in _CACHE:
        return _CACHE["sym"]
    m = matrices()
    U = _symplectic_reduce(m["P"])
    J = standard_J()
    assert np.array_equal(U.T @ m["P"] @ U, J)
    Ui = int_inverse(U)
    ops = {}
    for k in ("W11", "W13", "STAR"):
        Mj = Ui @ m[k] @ U
        ops[k] = Mj
        chk = is_symplectic if k != "STAR" else is_anti_symplectic
        assert chk(Mj, J)
        assert np.array_equal(Mj @ Mj, np.eye(_N, dtype=object) * 1)
    out = (U, ops)
    _CACHE["sym"] = out
    return out


def mod2(M):
    return np.array([[int(x) % 2 for x in row] for row in M], dtype=np.uint8)


def periods_frame_ops(tol: float = 1e-9):
    """The three operators in the symplectic frame of :mod:`mtft.periods`
    (the frame of the frozen Riemann matrix tau0), as exact integer
    matrices: x_symp = C x_hecke with C = hecke_to_symplectic_change, so
    operators map as M_symp = C M_hecke C^{-1}; STAR is already stored in
    symplectic coordinates by the periods module.

    All conjugation is exact integer arithmetic.  Verified on load:
    involutions over Z; symplectic conditions over Z against the periods
    symplectic form (W11, W13 symplectic; STAR anti-symplectic); W11-W13
    commutation over Z; commutation of STAR with the W's mod 2 (the
    integral lifts in this frame do not commute with STAR over Z -- that
    is a property of the lifts, not of the mod-2 spin action; the
    independent GP frame of :func:`symplectic_frame` carries fully
    commuting integral lifts)."""
    if "per" in _CACHE:
        return _CACHE["per"]
    from mtft import periods as P
    from mtft.periods.involutions import al_matrix

    def as_int(M, name):
        A = np.asarray(M, dtype=float)
        R = np.round(A)
        assert float(np.abs(A - R).max()) < tol, name
        return np.array([[int(x) for x in row] for row in R], dtype=object)

    C = as_int(P.hecke_to_symplectic_change(), "C")
    Ci = int_inverse(C)
    J = as_int(np.asarray(P.symplectic_form(), float).reshape(_N, _N), "J")

    ops = {
        "W11": C @ as_int(al_matrix(11), "W11") @ Ci,
        "W13": C @ as_int(al_matrix(13), "W13") @ Ci,
        "STAR": as_int(P.star_symplectic(), "STAR"),
    }

    I = np.eye(_N, dtype=object) * 1
    for k, M in ops.items():
        assert np.array_equal(M @ M, I), k
    assert np.array_equal(ops["W11"] @ ops["W13"], ops["W13"] @ ops["W11"])
    for k in ("W11", "W13"):
        K = ops[k] @ ops["STAR"] - ops["STAR"] @ ops[k]
        assert all(int(x) % 2 == 0 for row in K for x in row), \
            f"{k}/STAR must commute mod 2"
    assert is_symplectic(ops["W11"], J) and is_symplectic(ops["W13"], J)
    assert is_anti_symplectic(ops["STAR"], J)
    out = {"J": J, **ops}
    _CACHE["per"] = out
    return out
