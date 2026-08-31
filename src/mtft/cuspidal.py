"""mtft.cuspidal — cuspidal subgroups and Eisenstein torsion (v0.24.0).

Packages the Wave-8 arc (2026-08-29) as reusable, level-generic code.

The wave computed, for X_0(143), the cuspidal subgroup C(J_0(143)), its
2-torsion image in J[2], the mod-2 Eisenstein plane E_Eis, and established
that C[2] = E_Eis = J[m_2] — from which the lambda_2 residual Galois module
is SPLIT and its extension class vanishes.  Every step there is repeated
here with the level threaded through, so the same question can be asked at
other squarefree levels.  :func:`cross_level_control` does exactly that; it
is the answer to the level-universality risk flagged after Wave 7.

Two independent routes to C are implemented and agreement between them is a
gate, not a comment (protocol E2):

  route "projector"  Manin-Drinfeld splitting by a CRT-corrected spectral
                     projector pi = h(T_2) u(T_2), h the cuspidal charpoly
                     of T_2, u the truncated inverse of h at the Eisenstein
                     eigenvalue; C is read from the Abel-Jacobi coordinates
                     of the cusp divisors.
  route "cokernel"   no projector at all: C = Z^(c-1) / D(E cap Z^n), where
                     E is the saturated integral Eisenstein subspace and D
                     the boundary map.

CC-15 (2026-08-30).  The 2026-08-29 probe labelled the four cusps of
X_0(143) by index as (1, 11, 13, 143).  The true divisor labels in model
order are (143, 13, 11, 1): index 0 is the cusp infinity, which has width 1
and c = 0, hence gcd(c, N) = N.  The old names were reversed by d -> N/d.
Because the Atkin-Lehner divisor-toggle commutes with d -> N/d, every gate
of Wave 8 is unaffected and only the *names* of the generators change; the
corrected names are produced by :func:`cuspidal_group`, which reads its
labels from :func:`mtft.levels.cusp_labels`.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from functools import lru_cache
from math import gcd

import numpy as np

from . import levels as LV
from .integral_lattice import snf_transform, int_kernel, clear_denominators

__all__ = [
    "eisenstein_subspace", "cuspidal_group", "eisenstein_kernel_mod2",
    "two_torsion_image", "cross_level_control", "charpoly",
    "GOOD_PRIME_POOL",
]

GOOD_PRIME_POOL = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)


# ── small exact helpers ─────────────────────────────────────────────

def _obj(M):
    return np.array([[Fr(x) for x in row] for row in M], dtype=object)


def _good_primes(N, count=6):
    return [p for p in GOOD_PRIME_POOL if N % p][:count]


def _rational_solve(A, Y):
    """Solve A x = Y exactly for full-column-rank A.  Fractions in/out."""
    nr, nc = len(A), len(A[0])
    k = len(Y[0])
    M = [[Fr(A[i][j]) for j in range(nc)] + [Fr(Y[i][j]) for j in range(k)]
         for i in range(nr)]
    r = 0
    piv = []
    for c in range(nc):
        pr = next((i for i in range(r, nr) if M[i][c] != 0), None)
        if pr is None:
            raise ValueError("matrix is not full column rank")
        M[r], M[pr] = M[pr], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    for i in range(r, nr):
        if any(M[i][nc + j] != 0 for j in range(k)):
            raise ValueError("inconsistent system: target outside the span")
    return [[M[i][nc + j] for j in range(k)] for i in range(nc)]


def charpoly(M):
    """Characteristic polynomial (ascending coefficients) by Faddeev-LeVerrier.

    Exact over Q; the returned coefficients of an integer matrix are integers.
    """
    A = _obj(M)
    n = A.shape[0]
    I = np.eye(n, dtype=object) * Fr(1)
    Mk = np.zeros((n, n), dtype=object) * Fr(0)
    Mk = I.copy()
    coeffs = [Fr(1)]                       # leading, degree n
    Ak = A @ Mk
    for k in range(1, n + 1):
        c = -sum(Ak[i][i] for i in range(n)) / k
        coeffs.append(c)
        Mk = Ak + c * I
        Ak = A @ Mk
    asc = list(reversed(coeffs))
    for c in asc:
        assert c.denominator == 1, "non-integral charpoly coefficient"
    return [int(c) for c in asc]


def _polyval_mat(coeffs, X):
    n = X.shape[0]
    out = np.zeros((n, n), dtype=object)
    for c in reversed(coeffs):
        out = out @ X + np.eye(n, dtype=object) * Fr(c)
    return out


def _polyval(coeffs, x):
    v = Fr(0)
    for c in reversed(coeffs):
        v = v * x + c
    return v


def _shift_coeffs(poly, c, upto):
    """Coefficients of poly(x + c) up to degree *upto*."""
    from math import comb
    out = [Fr(0)] * (upto + 1)
    for k, a in enumerate(poly):
        for j in range(0, min(k, upto) + 1):
            out[j] += Fr(a) * comb(k, j) * Fr(c) ** (k - j)
    return out


def _m2(M):
    return np.array([[int(x) % 2 for x in row] for row in M], dtype=np.int64)


def _f2_kernel(A):
    A = A.copy() % 2
    m, n = A.shape
    piv, r = {}, 0
    for c in range(n):
        pr = next((i for i in range(r, m) if A[i, c]), None)
        if pr is None:
            continue
        A[[r, pr]] = A[[pr, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        piv[c] = r
        r += 1
    out = []
    for fc in [c for c in range(n) if c not in piv]:
        v = np.zeros(n, dtype=np.int64)
        v[fc] = 1
        for pc, ri in piv.items():
            v[pc] = A[ri, fc]
        out.append(v % 2)
    return out


def _f2_rank(A):
    A = np.asarray(A, dtype=np.int64).copy() % 2
    if A.size == 0:
        return 0
    m, n = A.shape
    r = 0
    for c in range(n):
        pr = next((i for i in range(r, m) if A[i, c]), None)
        if pr is None:
            continue
        A[[r, pr]] = A[[pr, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    return r


# ── Eisenstein subspace and cuspidal group ──────────────────────────

@lru_cache(maxsize=16)
def eisenstein_subspace(N: int, nprimes: int = 4):
    """Saturated integral Eisenstein subspace E of the relative model.

    E is the joint kernel of (T_p - (p+1)) over several good primes; T_p
    acts on the boundary as multiplication by p+1, so E is a lift of the
    degree-zero cusp lattice.  Returned as an (nq x (c-1)) integer array
    whose columns are a *primitive* basis.
    """
    m = LV.manin_model(N)
    nq = m["nq"]
    rows = []
    for p in _good_primes(N, nprimes):
        T = LV.hecke_matrix(N, p)
        rows.extend([[Fr(T[i][j]) - (Fr(p + 1) if i == j else 0)
                      for j in range(nq)] for i in range(nq)])
    ker = _f2 = None
    # exact rational kernel of the stacked system
    A = [[Fr(x) for x in row] for row in rows]
    mrows, n = len(A), nq
    piv, r = [], 0
    for c in range(n):
        pr = next((i for i in range(r, mrows) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(mrows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    basis = []
    for fc in [c for c in range(n) if c not in set(piv)]:
        v = [Fr(0)] * n
        v[fc] = Fr(1)
        for ri, pc in enumerate(piv):
            v[pc] = -A[ri][fc]
        basis.append(v)
    expect = m["ncusp"] - 1
    if len(basis) != expect:
        raise AssertionError(
            f"Eisenstein subspace has dimension {len(basis)}, expected "
            f"{expect} = #cusps - 1 at level {N}")
    ints, _ = clear_denominators(basis)
    # saturate: primitive basis from the SNF left transform
    U, S, _ = snf_transform(ints.T.tolist())
    Uinv = np.array(_int_inverse(U), dtype=object)
    return Uinv[:, :expect].copy()


def _int_inverse(U):
    """Exact inverse of a unimodular integer matrix."""
    n = U.shape[0]
    A = [[Fr(int(U[i][j])) for j in range(n)] +
         [Fr(1 if i == j else 0) for j in range(n)] for i in range(n)]
    for c in range(n):
        pr = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[pr] = A[pr], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    out = [[A[i][n + j] for j in range(n)] for i in range(n)]
    for row in out:
        for x in row:
            assert x.denominator == 1, "matrix was not unimodular"
    return [[int(x) for x in row] for row in out]


@lru_cache(maxsize=16)
def cuspidal_group(N: int) -> dict:
    """The cuspidal subgroup C(J_0(N)), by two independent routes.

    Returns structure, order, named generators (with corrected CC-15 cusp
    labels), the Abel-Jacobi order of each individual cusp difference, and
    the gate record.  Agreement of the two routes is a gate.
    """
    m = LV.manin_model(N)
    nq, ncusp = m["nq"], m["ncusp"]
    labels = m["cusp_label"]
    D = np.array(LV.boundary_matrix(N), dtype=object)
    gates = {}

    # ---- route "cokernel" -------------------------------------------
    E = eisenstein_subspace(N)
    DE = D @ E
    _, S_ck, _ = snf_transform(DE.tolist())
    inv_ck = [int(S_ck[i][i]) for i in range(min(S_ck.shape))
              if S_ck[i][i] != 0]
    order_ck = 1
    for s in inv_ck:
        order_ck *= s
    struct_ck = sorted(s for s in inv_ck if s > 1)

    # ---- route "projector" ------------------------------------------
    p0 = _good_primes(N, 1)[0]
    lam = p0 + 1
    Trel = np.array([[Fr(x) for x in row] for row in LV.hecke_matrix(N, p0)],
                    dtype=object)
    Tcusp = LV.cuspidal_hecke(N, p0)
    h = charpoly(Tcusp)
    h_at = _polyval(h, Fr(lam))
    gates["charpoly_nonvanishing"] = (h_at != 0)
    if h_at == 0:
        raise AssertionError(
            f"cuspidal charpoly vanishes at the Eisenstein eigenvalue {lam} "
            f"for p = {p0} at level {N}; try another prime")
    k = ncusp - 1
    sh = _shift_coeffs(h, lam, k - 1)
    # truncated inverse of h at x = lam to order k
    u = [Fr(1) / sh[0]]
    for j in range(1, k):
        acc = Fr(0)
        for i in range(1, j + 1):
            acc += sh[i] * u[j - i]
        u.append(-acc / sh[0])
    I = np.eye(nq, dtype=object) * Fr(1)
    Sh = Trel - lam * I
    uT = np.zeros((nq, nq), dtype=object)
    P = I.copy()
    for c in u:
        uT = uT + c * P
        P = P @ Sh
    pi = _polyval_mat(h, Trel) @ uT

    Kc = [[Fr(m["K"][j][i]) for j in range(len(m["K"]))] for i in range(nq)]
    gates["pi_idempotent"] = bool(np.array_equal(pi @ pi, pi))
    gates["pi_kills_cuspidal"] = all(
        x == 0 for x in (pi @ np.array(Kc, dtype=object)).flat)
    gates["D_pi_equals_D"] = bool(np.array_equal(D @ pi, D))

    # integral section of D
    U, S, V = snf_transform(D.tolist())
    sec_inv = [int(S[i][i]) for i in range(ncusp - 1)]
    gates["boundary_surjective"] = all(s == 1 for s in sec_inv)
    Vn = np.array(V, dtype=object)
    sigma = Vn[:, :ncusp - 1] @ np.array(
        [[Fr(1, sec_inv[i]) if i == j else Fr(0) for j in range(ncusp - 1)]
         for i in range(ncusp - 1)], dtype=object) @ np.array(U, dtype=object)
    gates["sigma_is_section"] = bool(
        np.array_equal(D @ sigma, np.eye(ncusp - 1, dtype=object)))

    Cq = (I - pi) @ sigma
    X = _rational_solve(Kc, [[Cq[i][j] for j in range(ncusp - 1)]
                             for i in range(nq)])
    d = 1
    for row in X:
        for x in row:
            d = d * x.denominator // gcd(d, x.denominator)
    Y = [[int(x * d) for x in row] for row in X]
    U2, S2, V2 = snf_transform(Y)
    sv = [int(S2[i][i]) if i < min(S2.shape) else 0 for i in range(ncusp - 1)]
    orders = [d // gcd(d, s) if s else 1 for s in sv]
    order_pr = 1
    for o in orders:
        order_pr *= o
    struct_pr = sorted(o for o in orders if o > 1)

    gates["two_route_agreement"] = (struct_pr == struct_ck
                                    and order_pr == order_ck)

    # named generators, CC-15 corrected labels
    base = labels[0]
    names = [f"C{labels[i]}" for i in range(1, ncusp)]
    V2n = np.array(V2, dtype=object)
    gens = {}
    for i, o in enumerate(orders):
        if o > 1:
            terms = " + ".join(
                f"{int(V2n[j][i])}*[{names[j]} - C{base}]"
                for j in range(ncusp - 1) if V2n[j][i] != 0)
            gens[o] = terms
    per_cusp = {}
    for j in range(ncusp - 1):
        e = [Fr(1) if i == j else Fr(0) for i in range(ncusp - 1)]
        pt = [sum(X[i][t] * e[t] for t in range(ncusp - 1)) for i in range(len(Kc[0]))]
        den = 1
        for z in pt:
            den = den * z.denominator // gcd(den, z.denominator)
        per_cusp[f"[{names[j]} - C{base}]"] = int(den)

    return dict(
        N=N, structure=struct_pr, order=order_pr,
        invariants_projector=struct_pr, invariants_cokernel=struct_ck,
        order_cokernel=order_ck, generators=gens, per_cusp_orders=per_cusp,
        cusp_labels=list(labels), gates=gates,
        coords=[[str(x) for x in row] for row in X], denominator=d,
        smith_projector=sv, smith_cokernel=inv_ck,
    )


# ── mod-2 Eisenstein plane and the two-torsion image ────────────────

@lru_cache(maxsize=16)
def eisenstein_kernel_mod2(N: int, nprimes: int = 8):
    """Joint kernel of {T_p - (1+p) mod 2} on H_1(X_0(N), F_2).

    For N = 143 this is the 2-dimensional plane E_Eis of Waves 6-8, and
    its dimension is exactly the mod-2 Eisenstein multiplicity.
    """
    g = LV.level_data(N)["genus"]
    dim = 2 * g
    rows = []
    for p in _good_primes(N, nprimes):
        T = LV.cuspidal_hecke(N, p)
        rows.append(_m2([[T[i][j] - ((1 + p) if i == j else 0)
                          for j in range(dim)] for i in range(dim)]))
    return [v.tolist() for v in _f2_kernel(np.vstack(rows))]


@lru_cache(maxsize=16)
def two_torsion_image(N: int) -> dict:
    """Image of C[2] in J[2] = H_1(X_0(N), F_2), with its relation to E_Eis.

    This is the Wave-8 comparison, level-generic: ``equals_eisenstein`` is
    the statement that drove the SPLIT verdict at 143.
    """
    cg = cuspidal_group(N)
    m = LV.manin_model(N)
    ncusp = m["ncusp"]
    X = [[Fr(s) for s in row] for row in cg["coords"]]
    d = cg["denominator"]
    Y = [[int(x * d) for x in row] for row in X]
    _, S2, V2 = snf_transform(Y)
    sv = [int(S2[i][i]) if i < min(S2.shape) else 0 for i in range(ncusp - 1)]
    orders = [d // gcd(d, s) if s else 1 for s in sv]
    V2n = np.array(V2, dtype=object)
    vecs = []
    for i, o in enumerate(orders):
        if o % 2 == 0:
            n = (o // 2) * V2n[:, i]
            pt = [sum(X[r][t] * Fr(int(n[t])) for t in range(ncusp - 1))
                  for r in range(len(X))]
            two = [2 * z for z in pt]
            assert all(z.denominator == 1 for z in two), "not 2-torsion"
            vecs.append([int(z) % 2 for z in two])
    E = eisenstein_kernel_mod2(N)
    dim_c2 = _f2_rank(vecs) if vecs else 0
    inside = (_f2_rank(list(E) + list(vecs)) == _f2_rank(E)) if vecs else True
    return dict(N=N, dim_C2_image=dim_c2, dim_E_Eis=len(E),
                contained_in_eisenstein=bool(inside),
                equals_eisenstein=bool(inside and dim_c2 == len(E)),
                vectors=vecs, eisenstein_basis=E)


def lambda2_congruence_scan(Lmax: int = 199, verbose: bool = False) -> dict:
    """Wave-8 open conjecture: a_ell = 1 + ell (mod lambda_2^2) for X_0(143).

    lambda_2 = (2, a_2 - 1) is the prime of O_K = Z[a_2] above 2 attached to
    the newform f_3 (the degree-6 block q6).  Wave 8 verified the congruence
    one level *above* the defining Eisenstein congruence for every good odd
    prime ell <= 199.  No proof is known: ordinary Sturm bounds control
    congruences mod lambda_2, not mod lambda_2^2, so this is a genuine open
    statement and every new prime is a real test.

    Returns the scan table and, crucially, ``counterexamples`` — the whole
    point of shipping this is that a future prime may falsify it.  The
    2-maximality of Z[a_2] (odd discriminant) is checked as a gate, since
    O/lambda_2^2 = Z/4 with a_2 -> 3 depends on it.
    """
    from . import hecke as _HK
    H6 = list(_HK.H6)
    disc = _discriminant(H6)
    if disc % 2 == 0:
        raise AssertionError(
            f"disc(h6) = {disc} is even; Z[a_2] need not be 2-maximal and "
            "the map O/lambda_2^2 = Z/4 is not justified")

    T2 = np.array([[Fr(x) for x in row] for row in _HK.cuspidal_hecke(2)],
                  dtype=object)
    dim = T2.shape[0]
    h6T = _polyval_mat(H6, T2)
    ker = []
    A = [[h6T[i][j] for j in range(dim)] for i in range(dim)]
    piv, r = [], 0
    for c in range(dim):
        pr = next((i for i in range(r, dim) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(dim):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    for fc in [c for c in range(dim) if c not in set(piv)]:
        v = [Fr(0)] * dim
        v[fc] = Fr(1)
        for ri, pc in enumerate(piv):
            v[pc] = -A[ri][fc]
        ker.append(v)
    ints, _ = clear_denominators(ker)
    U, _, _ = snf_transform(ints.T.tolist())
    B = np.array(_int_inverse(U), dtype=object)[:, :len(ker)]
    Bl = [[Fr(B[i][j]) for j in range(B.shape[1])] for i in range(dim)]

    def restrict(M):
        Y = [[sum(Fr(M[i][k]) * Bl[k][j] for k in range(dim))
              for j in range(B.shape[1])] for i in range(dim)]
        X = _rational_solve(Bl, Y)
        return np.array([[x for x in row] for row in X], dtype=object)

    R2 = restrict(_HK.cuspidal_hecke(2))
    nb = R2.shape[0]
    deg = len(H6) - 1
    pows = [np.eye(nb, dtype=object) * Fr(1)]
    for _ in range(deg - 1):
        pows.append(pows[-1] @ R2)
    cyc = None
    for j in range(nb):
        e = np.array([Fr(1) if i == j else Fr(0) for i in range(nb)],
                     dtype=object)
        cols = [[(P @ e)[i] for P in pows] for i in range(nb)]
        try:
            _rational_solve(cols, [[Fr(0)] for _ in range(nb)])
        except ValueError:
            continue
        cyc = (e, cols)
        break
    if cyc is None:
        raise AssertionError("no cyclic vector found for the q6 block")
    e, Acyc = cyc

    rows, counter = [], []
    for ell in _primes_upto(Lmax):
        if ell == 2 or 143 % ell == 0:
            continue
        Rl = restrict(_HK.cuspidal_hecke(ell))
        b = [[(Rl @ e)[i]] for i in range(nb)]
        cs = [row[0] for row in _rational_solve(Acyc, b)]
        M = np.zeros((nb, nb), dtype=object)
        for i, ci in enumerate(cs):
            M = M + pows[i] * ci
        assert np.array_equal(M, Rl), "polynomial identity for a_ell failed"
        val = Fr(0)
        for i, ci in enumerate(cs):
            val += ci * Fr(3) ** i
        if val.denominator % 2 == 0:
            raise AssertionError(
                f"even denominator {val.denominator} at ell = {ell}; "
                "the map to Z/4 is undefined")
        v4 = (val.numerator % 4) * pow(val.denominator % 4, -1, 4) % 4
        ok = ((ell + 1 - v4) % 4 == 0)
        rows.append(dict(ell=ell, a_ell_mod_lambda2sq=int(v4), holds=bool(ok)))
        if not ok:
            counter.append(ell)
        if verbose:
            print(f"  ell={ell:4d}  a_ell = {v4} (mod lambda_2^2)  "
                  f"{'ok' if ok else 'COUNTEREXAMPLE'}")
    return dict(Lmax=Lmax, disc_h6=disc, primes_tested=len(rows),
                counterexamples=counter, holds=(not counter), table=rows,
                status=("CONFIRMED to Lmax (open in general)" if not counter
                        else "FALSIFIED"))


def _primes_upto(n):
    s = [True] * (n + 1)
    out = []
    for p in range(2, n + 1):
        if s[p]:
            out.append(p)
            for q in range(p * p, n + 1, p):
                s[q] = False
    return out


def _discriminant(coeffs):
    """Discriminant of a monic integer polynomial via the resultant with f'."""
    n = len(coeffs) - 1
    f = [Fr(c) for c in coeffs]
    df = [Fr(i) * f[i] for i in range(1, len(f))]
    res = _resultant(f, df)
    sign = (-1) ** (n * (n - 1) // 2)
    val = sign * res / f[-1]
    assert val.denominator == 1
    return int(val)


def _resultant(a, b):
    a, b = list(a), list(b)
    res = Fr(1)
    while True:
        da, db = len(a) - 1, len(b) - 1
        if db < 0:
            return Fr(0)
        if db == 0:
            return res * b[0] ** da
        q_deg = da - db
        if q_deg < 0:
            a, b = b, a
            if (len(a) - 1) % 2 and (len(b) - 1) % 2:
                res = -res
            continue
        lead = b[-1]
        rem = list(a)
        for i in range(q_deg, -1, -1):
            c = rem[i + db] / lead
            if c:
                for j in range(db + 1):
                    rem[i + j] -= c * b[j]
        while rem and rem[-1] == 0:
            rem.pop()
        res = res * lead ** (len(a) - 1 - (len(rem) - 1 if rem else 0))
        if (len(a) - 1) % 2 and db % 2:
            res = -res
        a, b = b, rem
        if not rem:
            return Fr(0)


def cross_level_control(levels=(15, 33, 35, 105, 143)) -> list:
    """Replay the Wave-8 cuspidal/Eisenstein comparison across levels.

    This is the level-universality control: it decides whether the Wave-8
    identity C[2] = E_Eis is a general theorem about Eisenstein ideals or a
    property of 143.  Returns one row per level; no interpretation is
    attached here.

    NON-VACUITY (CC-16).  dim H_1 (x) F_2 = 2g, so at genus 1 the Eisenstein
    plane E_Eis (always 2-dimensional here) is the *whole* mod-2 homology and
    ``C2_equals_E_Eis`` is true for trivial reasons -- it tests nothing.  Each
    row therefore carries ``informative``, and callers counting evidence must
    count informative rows only.  N = 15 is the vacuous case; it is retained
    as a boundary check, not as support.

    N = 105 = 3*5*7 is the first control with three prime factors and eight
    cusps, so the cusp set is no longer a (Z/2)^2-torsor.  It carries the same
    genus as 143, which isolates cusp structure from genus.
    """
    out = []
    for N in levels:
        ld = LV.level_data(N)
        cg = cuspidal_group(N)
        tt = two_torsion_image(N)
        dim_h1_mod2 = 2 * ld["genus"]
        informative = tt["dim_E_Eis"] < dim_h1_mod2
        out.append(dict(
            N=N, genus=ld["genus"], ncusps=ld["ncusps"],
            n_prime_factors=len(ld["factors"]),
            structure=cg["structure"], order=cg["order"],
            routes_agree=cg["gates"]["two_route_agreement"],
            dim_H1_mod2=dim_h1_mod2,
            dim_E_Eis=tt["dim_E_Eis"], dim_C2=tt["dim_C2_image"],
            C2_equals_E_Eis=tt["equals_eisenstein"],
            informative=bool(informative),
            vacuity_reason=None if informative else
            "E_Eis is all of H_1 (x) F_2 at genus 1; identity is trivial",
        ))
    return out
