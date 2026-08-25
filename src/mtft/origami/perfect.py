"""mtft.origami.perfect — perfect t-embeddings and their branch structure.

Implements Galashin (9.2)-(9.3):

    t_i^lambda = <i+1, i-1> / ( <i-1, i> <i, i+1> )
    Theta(lambda)_i = (-1)^{i-1} t_i^lambda ( lambda_{1,i}, -lambda_{2,i} )^T

Lemma 9.2: Theta is an involution on Gr^diamond(2,n) with lambda _|_ Theta(lambda).
Lemma 9.3: on a valid pair, n = 2k and every t_i^lambda < 0.
Proposition 9.4: a t-immersion is PERFECT iff lambda-tilde = Theta(lambda).

The perfect system for a fixed C in Gr(3,n) is therefore

    Theta(lambda) C^T = 0,     lambda subset C,

which we parametrize by the annihilator: rows of A span n-perp, lambda = A C,
so the unknown is n in P^2.  (k = 3 case; documented as such.)
"""
from __future__ import annotations

import numpy as np
import sympy as sp

__all__ = [
    "bracket", "t_coefficients", "Theta", "winding", "is_valid_pair",
    "lambda_from_annihilator", "perfect_residual", "solve_perfect_branches",
    "cyclic_matrix", "orbit_structure", "equivariant_kasteleyn_factor",
]


# ----------------------------------------------------------------- brackets
def bracket(L, i, j, twist=1):
    """<i j>_L = det(L_i | L_j), 1-based, n-periodic with the twisted sign.

    For type (k, n) the paper's twisted cyclic symmetry is
    lambda_{i+n} = (-1)^{k-1} lambda_i; ``twist`` carries that sign
    (twist = +1 for odd k, e.g. k = 3).
    """
    n = L.shape[1]
    ii, jj = (i - 1) % n, (j - 1) % n
    si = twist ** ((i - 1) // n)
    sj = twist ** ((j - 1) // n)
    if hasattr(L, "det") or isinstance(L, sp.Matrix):
        return sp.expand(si * sj * (L[0, ii] * L[1, jj] - L[1, ii] * L[0, jj]))
    return si * sj * (L[0, ii] * L[1, jj] - L[1, ii] * L[0, jj])


def t_coefficients(L, twist=1):
    """[t_1, ..., t_n] per (9.2)."""
    n = L.shape[1]
    return [bracket(L, i + 1, i - 1, twist)
            / (bracket(L, i - 1, i, twist) * bracket(L, i, i + 1, twist))
            for i in range(1, n + 1)]


def Theta(L, twist=1):
    """The involution (9.3).  Accepts and returns sympy or numpy 2 x n."""
    n = L.shape[1]
    ts = t_coefficients(L, twist)
    rows = [[], []]
    for i in range(1, n + 1):
        s = (-1) ** (i - 1)
        rows[0].append(s * ts[i - 1] * L[0, i - 1])
        rows[1].append(-s * ts[i - 1] * L[1, i - 1])
    if isinstance(L, sp.Matrix):
        return sp.Matrix(rows)
    return np.array(rows, dtype=float)


def winding(L, twist=1):
    """Total turning angle of the columns, in units of pi."""
    n = L.shape[1]
    tot = 0.0
    for i in range(1, n + 1):
        u = (float(L[0, (i - 1) % n]), float(L[1, (i - 1) % n]))
        v0, v1 = (i % n), None
        v = (float(L[0, v0]), float(L[1, v0]))
        if i == n:
            v = (twist * v[0], twist * v[1])
        tot += np.arctan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1])
    return tot / np.pi


def is_valid_pair(L, twist=1):
    """Sign-flip validity: <i i+1> > 0, [i i+1]_Theta > 0, all t_i < 0."""
    n = L.shape[1]
    T = Theta(L, twist)
    b = [bracket(L, i, i + 1, twist) for i in range(1, n + 1)]
    bt = [bracket(T, i, i + 1, twist) for i in range(1, n + 1)]
    ts = t_coefficients(L, twist)
    return (all(float(x) > 0 for x in b)
            and all(float(x) > 0 for x in bt)
            and all(float(x) < 0 for x in ts))


# ------------------------------------------------------- the perfect system
def lambda_from_annihilator(nvec, C):
    """lambda = A C where the rows of A span nvec-perp (k = 3)."""
    n0, n1, n2 = nvec
    if isinstance(C, sp.Matrix):
        A = sp.Matrix([[n2, 0, -n0], [0, n2, -n1]])
    else:
        A = np.array([[n2, 0, -n0], [0, n2, -n1]], dtype=float)
    return A @ C if not isinstance(C, sp.Matrix) else A * C


def perfect_residual(nvec, C):
    """max |Theta(lambda) C^T| — zero exactly on a perfect branch."""
    L = lambda_from_annihilator(nvec, C)
    try:
        T = Theta(L)
    except ZeroDivisionError:
        return 1e12
    if not np.all(np.isfinite(T)):
        return 1e12
    return float(np.abs(T @ C.T).max())


def solve_perfect_branches(C, n_starts=4000, seed=20260825, tol=1e-9):
    """All admissible perfect branches for C in Gr(3,n).

    Returns list of dicts with keys: n (normalized annihilator), lam, theta,
    valid, wind_lam, wind_theta, residual.
    """
    from scipy.optimize import fsolve

    def F(v):
        L = lambda_from_annihilator([v[0], v[1], 1.0], C)
        try:
            T = Theta(L)
        except ZeroDivisionError:
            return [1e6, 1e6]
        if not np.all(np.isfinite(T)):
            return [1e6, 1e6]
        M = T @ C.T
        return [M[0, 0] + M[1, 1], M[0, 1] - M[1, 2]]

    rng = np.random.default_rng(seed)
    found = []
    for _ in range(n_starts):
        s = rng.uniform(-6, 6, 2)
        try:
            x, _info, ier, _msg = fsolve(F, s, full_output=True)
        except Exception:
            continue
        if ier != 1:
            continue
        if perfect_residual([x[0], x[1], 1.0], C) > tol:
            continue
        if not any(np.allclose(x, y, atol=1e-6) for y in found):
            found.append(x)

    out = []
    for x in found:
        nv = np.array([x[0], x[1], 1.0])
        nv = nv / np.abs(nv).max()
        L = lambda_from_annihilator(nv, C)
        out.append({
            "n": nv, "lam": L, "theta": Theta(L),
            "valid": is_valid_pair(L),
            "wind_lam": round(winding(L), 6),
            "wind_theta": round(winding(Theta(L)), 6),
            "residual": perfect_residual(nv, C),
        })
    return out


# ------------------------------------------------------------- cyclic action
def cyclic_matrix(C, shift):
    """Matrix of the column-shift action on the row-coordinate space of C."""
    n = C.shape[1]
    Csh = np.array([[C[r, (i - shift) % n] for i in range(n)]
                    for r in range(C.shape[0])], dtype=float)
    return (Csh @ C.T) @ np.linalg.inv(C @ C.T)


def orbit_structure(branches, R, atol=1e-5):
    """Classify branches under the cyclic action R: returns (fixed, cycles)."""
    def key(v):
        v = np.asarray(v, dtype=float)
        return v / np.abs(v).max()

    keys = [key(b["n"]) for b in branches]

    def find(v):
        for j, k in enumerate(keys):
            if np.allclose(np.abs(key(v)), np.abs(k), atol=atol):
                return j
        return None

    perm = [find(R @ b["n"]) for b in branches]
    fixed = [i for i, p in enumerate(perm) if p == i]
    cycles, seen = [], set(fixed)
    for i in range(len(branches)):
        if i in seen or perm[i] is None:
            continue
        cyc, j = [], i
        while j is not None and j not in seen:
            cyc.append(j)
            seen.add(j)
            j = perm[j]
        if len(cyc) > 1:
            cycles.append(cyc)
    return fixed, cycles, perm


# ------------------------------------- C3-equivariant Kasteleyn factorization
def equivariant_kasteleyn_factor(K, order=3):
    """Block-diagonalize a Kasteleyn matrix under a cyclic symmetry.

    ``K`` must be indexed so that rows/columns run over ``order``-orbits in
    rotation order, orbit-major.  Returns (blocks, dets, product_identity),
    where the determinant factors as prod(dets) — the trivial-character factor
    times a norm from the conjugate characters.
    """
    m = K.shape[0]
    orbits = m // order
    om = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2 if order == 3 else \
        sp.exp(2 * sp.pi * sp.I / order)
    U = sp.zeros(m, m)
    for k in range(order):
        for j in range(order):
            for o in range(orbits):
                U[o * order + j, o * order + k] = om ** (k * j) / sp.sqrt(order)
    Kb = sp.simplify(sp.conjugate(U.T) * K * U)
    blocks, dets = [], []
    for k in range(order):
        idx = [o * order + k for o in range(orbits)]
        blk = Kb[idx, idx]
        blocks.append(sp.simplify(blk))
        dets.append(sp.expand(sp.simplify(blk.det())))
    prod_ok = sp.expand(sp.simplify(K.det() - sp.prod(dets))) == 0
    return blocks, dets, prod_ok
