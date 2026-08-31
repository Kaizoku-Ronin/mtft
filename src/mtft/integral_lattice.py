"""mtft.lattice — exact integer-lattice toolkit (v0.19.0).

Pure Python integers + numpy object arrays; no sympy.  Everything here
was battle-tested in the 2026-08-24 integral-model arc (certificates
v1-v9): p-saturation with unimodular pivots, Smith/Hermite forms,
quotient invariants, HNF membership/order tests, and operator matrices
on lattices.

Conventions.  A "lattice" is a numpy object array of shape (m, n) whose
COLUMNS generate a sublattice of Z^m.  All routines are deterministic:
saturation pivots are normalized to c_j = 1 (unimodular away from p),
and prime lists are processed in the given order.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd

import numpy as np

__all__ = [
    "rank_modp", "kernel_modp", "p_saturate", "saturate",
    "smith_invariants", "hnf", "quotient_invariants",
    "solve_in_lattice", "class_order", "rational_kernel",
    "operator_matrix",
    # v0.24.0 additions
    "snf_transform", "int_kernel", "clear_denominators",
]


class InexactInputError(TypeError):
    """Raised when an exact-integer routine is handed a non-integral value.

    Corrections CC-11..CC-14 (2026-08-24 / 2026-08-29).  The v0.23.0 tree
    coerced inputs with ``int(x)``, which *truncates* Fractions silently:
    ``int(Fraction(1, 2)) == 0``.  A rational vector fed to ``saturate``,
    ``solve_in_lattice``, ``rational_kernel`` or ``class_order`` was
    therefore zeroed without warning, and the guard inside ``saturate``
    was dead because ``_as_obj`` had already destroyed the evidence.
    Every coercion path now raises instead of truncating.
    """


def _exact(x):
    """Coerce to Fraction without ever truncating.  Floats must be integral."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, (int, np.integer)):
        return Fraction(int(x))
    if isinstance(x, float):
        if x != int(x):
            raise InexactInputError(
                f"non-integral float {x!r}; exact routines take int or Fraction")
        return Fraction(int(x))
    num, den = getattr(x, "numerator", None), getattr(x, "denominator", None)
    if num is not None and den is not None:
        return Fraction(int(num), int(den))
    raise InexactInputError(f"cannot exactly coerce {type(x).__name__}: {x!r}")


def _as_int(x):
    """Exact integer coercion; raises on anything with a denominator > 1."""
    f = _exact(x)
    if f.denominator != 1:
        raise InexactInputError(
            f"expected an integer, got {f}; clear denominators first "
            "(see integral_lattice.clear_denominators)")
    return int(f)


def _as_obj(B):
    """Integer object array.  Raises (never truncates) on rational input."""
    return np.array([[_as_int(x) for x in row] for row in B], dtype=object)


def _as_frac(B):
    return [[_exact(x) for x in row] for row in B]


def clear_denominators(vectors, by_column=False):
    """Scale Fraction vectors to primitive integer vectors.

    Companion to CC-11: ``rational_kernel`` returns Fraction vectors and
    ``saturate`` requires integers, so this is the sanctioned bridge.
    Returns ``(integer_array, multipliers)`` where multiplier[i] is the
    factor applied to vector i.
    """
    rows = [[_exact(x) for x in v] for v in vectors]
    mults, out = [], []
    for v in rows:
        den = 1
        for x in v:
            den = den * x.denominator // gcd(den, x.denominator)
        w = [int(x * den) for x in v]
        g = 0
        for t in w:
            g = gcd(g, abs(t))
        if g > 1:
            w = [t // g for t in w]
            den = Fraction(den, g)
        mults.append(den)
        out.append(w)
    arr = np.array(out, dtype=object)
    return (arr.T if by_column else arr), mults


def kernel_modp(B, p):
    """(rank, kernel basis) of the columns of B over F_p (p prime)."""
    Bp = [[int(x) % p for x in row] for row in B]
    m, n = len(Bp), len(Bp[0])
    piv, r = [], 0
    for c in range(n):
        pr = next((i for i in range(r, m) if Bp[i][c] % p), None)
        if pr is None:
            continue
        Bp[r], Bp[pr] = Bp[pr], Bp[r]
        inv = pow(Bp[r][c], p - 2, p)
        Bp[r] = [(x * inv) % p for x in Bp[r]]
        for i in range(m):
            if i != r and Bp[i][c]:
                f = Bp[i][c]
                Bp[i] = [(x - f * y) % p for x, y in zip(Bp[i], Bp[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(n) if c not in piv]
    ker = []
    for fc in free:
        v = [0] * n
        v[fc] = 1
        for ri, pc in enumerate(piv):
            v[pc] = (-Bp[ri][fc]) % p
        ker.append(v)
    return r, ker


def rank_modp(B, p):
    return kernel_modp(B, p)[0]


def p_saturate(B, p):
    """(B_saturated_at_p, steps).  Column ops unimodular away from p."""
    B = _as_obj(B)
    steps = 0
    while True:
        r, ker = kernel_modp(B, p)
        if r == B.shape[1]:
            return B, steps
        c = ker[0]
        j = next(j for j, x in enumerate(c) if x % p)
        inv = pow(c[j], p - 2, p)
        c = [(x * inv) % p for x in c]
        c[j] = 1
        v = B @ np.array(c, dtype=object)
        assert all(int(x) % p == 0 for x in v), "saturation invariant broken"
        B[:, j] = np.array([int(x) // p for x in v], dtype=object)
        steps += 1


def saturate(B, primes):
    """Saturate at each prime in order.  Returns (B_sat, {p: steps}).

    CC-11: the denominator check runs on the *caller's* matrix, before any
    coercion.  In v0.23.0 it ran after ``_as_obj`` had already truncated,
    so it could never fire.
    """
    for row in B:
        for x in row:
            if getattr(x, "denominator", 1) != 1:
                raise InexactInputError(
                    "saturate expects an integer matrix; clear denominators "
                    "first (rational_kernel returns Fraction vectors; use "
                    "integral_lattice.clear_denominators)")
    B = _as_obj(B)
    log = {}
    for p in primes:
        B, s = p_saturate(B, p)
        if s:
            log[p] = s
    return B, log


def smith_invariants(B):
    """Nonzero Smith invariant factors of the integer matrix B.

    HNF ping-pong: alternate column-HNF on the matrix and its transpose
    until diagonal, then repair the divisibility chain by pairwise
    gcd/lcm.  Avoids the entry explosion of naive two-sided elimination.
    """
    M = _as_obj(B)
    if M.shape[0] == 0 or M.shape[1] == 0:
        return []
    for _ in range(60):
        M = hnf(M)
        if M.shape[1] == 0:
            return []
        diag = all(i == j or M[i][j] == 0
                   for i in range(M.shape[0]) for j in range(M.shape[1]))
        if diag:
            break
        M = hnf(M.T)
    else:
        raise RuntimeError("smith ping-pong did not converge")
    d = [abs(int(M[i][i])) for i in range(min(M.shape))]
    d = [x for x in d if x]
    for i in range(len(d)):
        for j in range(i + 1, len(d)):
            g = gcd(d[i], d[j])
            d[i], d[j] = g, d[i] // g * d[j]
    return sorted(d)


def hnf(B):
    """Column-style Hermite normal form basis of the column lattice.

    Returns an (m, r) object array, columns echelon with positive pivots.
    """
    A = _as_obj(B)
    m, n = A.shape
    cols = [list(A[:, j]) for j in range(n)]
    basis = []
    row = 0
    while cols and row < m:
        live = [c for c in cols if any(c[row:])]
        if not live:
            break
        while True:
            nz = [c for c in cols if c[row] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda c: abs(c[row]))
            a = nz[0]
            for c in nz[1:]:
                q = c[row] // a[row]
                for i in range(m):
                    c[i] -= q * a[i]
        nz = [c for c in cols if c[row] != 0]
        if nz:
            piv = nz[0]
            if piv[row] < 0:
                for i in range(m):
                    piv[i] = -piv[i]
            basis.append(piv)
            cols = [c for c in cols if c is not piv]
        row += 1
    return np.array(basis, dtype=object).T if basis else np.zeros((m, 0), dtype=object)


def solve_in_lattice(H, v):
    """Solve H x = v for x over Q given an HNF basis H.  None if v is
    outside the Q-span; otherwise the Fraction vector x (integral iff
    v lies in the lattice)."""
    m, r = H.shape
    v = [_exact(t) for t in v]          # CC-12: was Fraction(int(t)) — truncating
    x = [Fraction(0)] * r
    row = 0
    for j in range(r):
        while row < m and H[row][j] == 0:
            if v[row] != 0:
                return None
            row += 1
        if row == m:
            return None
        x[j] = v[row] / int(H[row][j])
        for i in range(m):
            v[i] -= x[j] * int(H[i][j])
        row += 1
    if any(t != 0 for t in v):
        return None
    return x


def class_order(H, v, divisors):
    """Smallest d in sorted(divisors) with d*v in the lattice, else None."""
    for d in sorted(divisors):
        # CC-14: was d * int(t) — truncating on Fraction input
        x = solve_in_lattice(H, [d * _exact(t) for t in v])
        if x is not None and all(t.denominator == 1 for t in x):
            return d
    return None


def rational_kernel(M):
    """Basis of the rational kernel of the columns of M (list of Fraction
    vectors), by exact Gauss elimination.  For the 141 x 15 CI-A systems."""
    A = _as_frac(M)                    # CC-13: was Fraction(int(x)) — truncating
    m, n = len(A), len(A[0])
    piv, r = [], 0
    for c in range(n):
        pr = next((i for i in range(r, m) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(n) if c not in piv]
    ker = []
    for fc in free:
        v = [Fraction(0)] * n
        v[fc] = Fraction(1)
        for ri, pc in enumerate(piv):
            v[pc] = -A[ri][fc]
        ker.append(v)
    return ker


def snf_transform(A):
    """Smith normal form *with* transforms: returns (U, S, V), U A V = S.

    U and V are unimodular; S is diagonal with the divisibility chain
    s_1 | s_2 | ... .  Unlike :func:`smith_invariants` this keeps the
    change-of-basis data, which is what quotient-group *generators*
    (as opposed to orders) require.  Promoted from the Wave-8 toolkit.
    """
    M = [[_as_int(x) for x in row] for row in A]
    m, n = len(M), (len(M[0]) if M else 0)
    U = [[int(i == j) for j in range(m)] for i in range(m)]
    V = [[int(i == j) for j in range(n)] for i in range(n)]

    def swap_rows(X, i, j):
        X[i], X[j] = X[j], X[i]

    def swap_cols(X, i, j):
        for r in X:
            r[i], r[j] = r[j], r[i]

    def addrow(X, d, s, c):
        X[d] = [a + c * b for a, b in zip(X[d], X[s])]

    def addcol(X, d, s, c):
        for r in X:
            r[d] = r[d] + c * r[s]

    t = 0
    while t < min(m, n):
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if M[i][j] != 0 and (piv is None
                                     or abs(M[i][j]) < abs(M[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None:
            break
        i, j = piv
        if i != t:
            swap_rows(M, t, i), swap_rows(U, t, i)
        if j != t:
            swap_cols(M, t, j), swap_cols(V, t, j)
        prog = True
        while prog:
            prog = False
            for i in range(t + 1, m):
                if M[i][t]:
                    q = M[i][t] // M[t][t]
                    addrow(M, i, t, -q), addrow(U, i, t, -q)
                    if M[i][t]:
                        swap_rows(M, t, i), swap_rows(U, t, i)
                        prog = True
            for j in range(t + 1, n):
                if M[t][j]:
                    q = M[t][j] // M[t][t]
                    addcol(M, j, t, -q), addcol(V, j, t, -q)
                    if M[t][j]:
                        swap_cols(M, t, j), swap_cols(V, t, j)
                        prog = True
        fixed = False
        while not fixed:
            fixed = True
            for i in range(t + 1, m):
                if any(M[i][j] % M[t][t] for j in range(t + 1, n)):
                    addrow(M, t, i, 1), addrow(U, t, i, 1)
                    fixed = False
                    prog = True
                    while prog:
                        prog = False
                        for ii in range(t + 1, m):
                            if M[ii][t]:
                                q = M[ii][t] // M[t][t]
                                addrow(M, ii, t, -q), addrow(U, ii, t, -q)
                                if M[ii][t]:
                                    swap_rows(M, t, ii), swap_rows(U, t, ii)
                                    prog = True
                        for jj in range(t + 1, n):
                            if M[t][jj]:
                                q = M[t][jj] // M[t][t]
                                addcol(M, jj, t, -q), addcol(V, jj, t, -q)
                                if M[t][jj]:
                                    swap_cols(M, t, jj), swap_cols(V, t, jj)
                                    prog = True
                    break
        if M[t][t] < 0:
            M[t] = [-x for x in M[t]]
            U[t] = [-x for x in U[t]]
        t += 1
    return (np.array(U, dtype=object), np.array(M, dtype=object),
            np.array(V, dtype=object))


def int_kernel(A):
    """Basis (as columns) of the integer kernel {x in Z^n : A x = 0}.

    Saturated by construction: obtained from the trailing columns of the
    SNF right transform, so the returned lattice is primitive.
    """
    _, S, V = snf_transform(A)
    m, n = S.shape
    r = sum(1 for i in range(min(m, n)) if S[i][i] != 0)
    return V[:, r:].copy()


def quotient_invariants(ambient, sub):
    """Smith invariants of ambient/sub for two bases of the same Q-span.

    Solves each sub column in the ambient basis (must be integral) and
    returns the invariant factors of the coordinate matrix.
    """
    H = hnf(ambient)
    X = []
    for j in range(sub.shape[1]):
        x = solve_in_lattice(H, list(sub[:, j]))
        assert x is not None, "sub not in ambient span"
        assert all(t.denominator == 1 for t in x), "sub not integral in ambient"
        X.append([int(t) for t in x])
    return smith_invariants(np.array(X, dtype=object).T)


def operator_matrix(L, reference, diag_signs):
    """Matrix (and content denominator) of the operator that is
    diag(diag_signs) in the reference basis, expressed on the lattice L.

    Both L and reference are (m, k) bases of the same Q-span.  Exact
    normal-equations route (Gram matrices over Z, Fraction solves):
    C = (R^T R)^{-1} R^T L, then Y = (L^T L)^{-1} L^T (R diag C).
    Returns (Y as list of Fraction columns, lcm of denominators).
    """
    m, k = L.shape
    R = reference

    def gram_solve(Bm, T):
        G = [[sum(int(Bm[r][i]) * int(Bm[r][j]) for r in range(m))
              for j in range(k)] for i in range(k)]
        RHS = [[T[i][j] for j in range(len(T[0]))] for i in range(len(T))]
        n = k
        aug = [[Fraction(G[i][j]) for j in range(n)] + list(RHS[i])
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
        return [[aug[i][n + j] for j in range(len(RHS[0]))]
                for i in range(n)]

    RtL = [[sum(Fraction(int(R[r][i])) * int(L[r][j]) for r in range(m))
            for j in range(k)] for i in range(k)]
    C = gram_solve(R, RtL)                       # L = R C
    W = [[sum(Fraction(int(R[r][i])) * diag_signs[i] * C[i][j]
              for i in range(k)) for j in range(k)] for r in range(m)]
    LtW = [[sum(Fraction(int(L[r][i])) * W[r][j] for r in range(m))
            for j in range(k)] for i in range(k)]
    Y = gram_solve(L, LtW)                       # W = L Y
    den = 1
    for row in Y:
        for t in row:
            den = den * t.denominator // gcd(den, t.denominator)
    return Y, den
