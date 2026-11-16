r"""mtft.lines.periods — the 26 periods of a class by Hecke block, by two routes, and the lattice of 143a1 (LN-01, 2026-11-16).

H₁ of the punctured curve is ℤ²⁹ (``mtft.lines.darts``).  Three of its dimensions wind round the
cusps (``mtft.lines.rademacher.windings``).  The other 26 are H₁ of the closed curve X₀(143), on
which the Hecke algebra acts with four blocks, the kernels of the irreducible factors of the
characteristic polynomial x²(x + 2)⁴g₄²h₆² of T₂:

    ell   2   143a1                      old   4   the level-11 oldspace (the ghost of X₀(11))
    q4    8   the quartic newform f2     q6   12   the sextic newform f3

The component of a class in a block is its *period* there: the integrals of that block's forms
along the class.  Projected along the other blocks, H₁(X₀(143), ℤ) is a lattice in each block, and
in a ℤ-basis of that lattice a class has 2 + 4 + 8 + 12 integer coordinates.

Route 1 (EXACT).  ``coordinates(cls)``: the lab's integer coordinates of the 29 generators (frozen,
``mtft.lines.table.coord``) applied to the class.  They are the numbers of the lab's table.

Route 2 (EXACT).  ``components(γ)``: the modular symbol {∞, γ∞} written in Manin symbols by the
continued fraction of γ∞ and resolved in the package's own model, ``mtft.hecke`` (its own
relations, its own blocks, no file).  The bases differ, so the routes are compared through what no
basis changes: ``compare_bases`` solves for the matrix between them, block by block, exactly, and
checks that the lab's coordinates of the generators span all of ℤ^dim.

The lattice of 143a1 (CERTIFIED, 1e-8).  ``lattice_143a1`` integrates the newform of the curve
y² + y = x³ − x² − x − 2 along twelve small elements, with coefficients from point counts, and
solves for Ω₁, Ω₂ with ∫_γ 2πi f dz = n₁Ω₁ + n₂Ω₂, (n₁, n₂) the class's coordinates on the block
``ell``.  The second route is the package's own frozen period matrix (``mtft.periods``, from
PARI/GP): ``stored_periods_in_lattice``.

Provenance: the Modular Plane lab, record MP12 (12 November 2026); its Route B ran on ``mtft.hecke``.
"""
from __future__ import annotations

import math
from fractions import Fraction as Fr
from functools import lru_cache

import numpy as np

from . import darts as _d
from . import table as _t

BLOCKS = _t.BLOCKS
BLOCK_DIMS = _t.BLOCK_DIMS
CURVE_J = Fr(64 ** 3, -1859)
"""j of 143a1: c4³/Δ for y² + y = x³ − x² − x − 2 (c4 = 64, Δ = −1859 = −11·13²)."""


# ------------------------------------------------------------------ route 1: the lab's coordinates
def coordinates(cls) -> dict:
    """The integer coordinates of a class of ℤ²⁹ on the four blocks, in the lab's bases."""
    C = _t.coord()
    if len(cls) != 29:
        raise ValueError("a class of level 143 has 29 coordinates")
    return {name: tuple(sum(int(cls[e]) * C[name][e][i] for e in range(29)) for i in range(BLOCK_DIMS[name])) for name in BLOCKS}


def coordinates_of(gamma) -> dict:
    """The block coordinates of an element of Γ₀(143)."""
    return coordinates(_d.graph(143).cls(_d.as_matrix(gamma)))


def flat(coords: dict) -> tuple:
    """The 26 coordinates in the order ell, old, q4, q6."""
    return tuple(v for name in BLOCKS for v in coords[name])


def thirty(gamma) -> dict:
    """The thirty integers of an element of Γ₀(143): Ψ, the three windings and the 26 periods by block."""
    from . import rademacher as _r
    gamma = _d.as_matrix(gamma)
    cls = _d.graph(143).cls(gamma)
    co = coordinates(cls)
    return {"psi": _r.Psi(gamma), "windings": _r.windings(gamma), "periods": co, "cls": cls,
            "vector": (_r.Psi(gamma),) + _r.windings(gamma) + flat(co)}


# ------------------------------------------------------------------ route 2: the package's own Manin symbols
def manin_symbol(g) -> list:
    """{∞, g∞} as a vector on the 84 edges of ``mtft.hecke``'s model: the Manin symbols of the convergents of g∞."""
    import mtft.hecke as H
    m = H.model()
    N, E = m["N"], m["E"]
    idx, canon, eid, esign = m["idx"], m["canon"], m["eid"], m["esign"]
    (a, _), (c, _) = g
    if c < 0:
        a, c = -a, -c
    v = [0] * E

    def put(u, w):
        y = idx[canon(u % N, w % N)]
        v[eid[y]] += esign[y]

    if c:
        pm, qm, p_, q_ = 1, 0, a // c, 1                                # convergents p_{-1}/q_{-1}, p_0/q_0
        num, den = a - (a // c) * c, c
        put(-1, 0)                                                      # j = 0: the symbol (-q_0 : q_{-1}) = (-1 : 0)
        j = 0
        while num:
            j += 1
            aj = den // num
            num, den = den - aj * num, num
            pm, qm, p_, q_ = p_, q_, aj * p_ + pm, aj * q_ + qm
            put((-1) ** (j - 1) * q_, qm)
        if (p_, q_) != (a, c):
            raise ArithmeticError("the continued fraction did not end at the entry it started from")
    return v


def _inverse(M):
    n = len(M)
    A = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        pr = next((i for i in range(c, n) if A[i][c] != 0), None)
        if pr is None:
            raise ArithmeticError("the matrix is singular")
        A[c], A[pr] = A[pr], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


@lru_cache(maxsize=1)
def _winv():
    import mtft.hecke as H
    blocks = H.blocks()
    cols = [vec for name in BLOCKS for vec in blocks[name]]             # 26 block vectors, each of 26 cuspidal coordinates
    if [len(blocks[name]) for name in BLOCKS] != [BLOCK_DIMS[name] for name in BLOCKS]:
        raise ArithmeticError("the Hecke blocks do not have dimensions 2, 4, 8, 12")
    W = [[cols[j][i] for j in range(26)] for i in range(26)]
    return _inverse(W)


def components(g) -> dict:
    """The components of the class of an element of Γ₀(143) on the four blocks, in ``mtft.hecke``'s bases (Fractions)."""
    import mtft.hecke as H
    g = _d.as_matrix(g)
    if g[1][0] % 143:
        raise ValueError("the matrix is not in the group of level 143")
    m = H.model()
    q = H._qcoords(m, [Fr(x) for x in manin_symbol(g)])
    K, freeK, nq = m["K"], m["freeK"], m["nq"]
    back = [sum(K[j][i] * q[freeK[j]] for j in range(len(K))) for i in range(nq)]
    if any(back[i] != q[i] for i in range(nq)):
        raise ArithmeticError("the class is not cuspidal in the package's model")
    w26 = [q[fj] for fj in freeK]
    comp = [sum(row[j] * w26[j] for j in range(26)) for row in _winv()]
    out, off = {}, 0
    for name in BLOCKS:
        out[name] = tuple(comp[off:off + BLOCK_DIMS[name]])
        off += BLOCK_DIMS[name]
    return out


@lru_cache(maxsize=1)
def generator_components() -> dict:
    """Block name → 29 rows: the components of the 29 generators by route 2."""
    comps = [components(g) for g in _d.graph(143).generators]
    return {name: tuple(c[name] for c in comps) for name in BLOCKS}


def _rref(rows):
    """Reduced row echelon form over the rationals: (rows, pivot columns)."""
    A = [[Fr(x) for x in r] for r in rows]
    piv, r = [], 0
    ncol = len(A[0]) if A else 0
    for c in range(ncol):
        pr = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == len(A):
            break
    return A, piv


def rank(rows) -> int:
    """The rank of a list of rational rows (exact)."""
    return len(_rref(rows)[1]) if rows else 0


def relations(rows) -> list:
    """The reduced basis of the linear relations among the rows (exact): vectors v with Σ v_i row_i = 0."""
    if not rows:
        return []
    n, m = len(rows), len(rows[0])
    red, piv = _rref([[rows[i][j] for i in range(n)] for j in range(m)])       # the columns are the rows given
    ps = set(piv)
    ns = []
    for fc in [c for c in range(n) if c not in ps]:
        v = [Fr(0)] * n
        v[fc] = Fr(1)
        for r_, c_ in enumerate(piv):
            v[c_] = -red[r_][fc]
        ns.append(v)
    if not ns:
        return []
    out, _ = _rref(ns)
    return [tuple(r) for r in out if any(r)]


@lru_cache(maxsize=1)
def compare_bases() -> dict:
    """Route 1 against route 2 on the 29 generators, block by block.

    For each block: the matrix U with (lab's coordinates) U = (package's components), solved on dim
    independent generators and checked on all 29; whether it is invertible; and the Smith invariants
    of the lab's 29 × dim integer matrix (all 1 exactly when its rows span ℤ^dim).
    """
    from mtft.integral_lattice import smith_invariants
    lab, own = _t.coord(), generator_components()
    out = {}
    for name in BLOCKS:
        dm = BLOCK_DIMS[name]
        rows, pick = [], []
        for e in range(29):                                             # dim independent rows of the lab's matrix
            if rank(rows + [lab[name][e]]) > len(rows):
                rows.append(lab[name][e])
                pick.append(e)
            if len(rows) == dm:
                break
        if len(rows) != dm:
            out[name] = {"rank": len(rows), "consistent": False, "invertible": False, "smith": []}
            continue
        Linv = _inverse(rows)
        U = [[sum(Linv[i][k] * own[name][pick[k]][j] for k in range(dm)) for j in range(dm)] for i in range(dm)]
        consistent = all(tuple(sum(lab[name][e][i] * U[i][j] for i in range(dm)) for j in range(dm)) == tuple(own[name][e])
                         for e in range(29))
        smith = [int(v) for v in smith_invariants([list(r) for r in lab[name]])]
        out[name] = {"rank": dm, "consistent": consistent, "invertible": rank(U) == dm, "smith": smith,
                     "primitive": smith == [1] * dm}
    return out


# ------------------------------------------------------------------ the lattice of 143a1, numerically
def ap_143a1(p: int) -> int:
    """a_p of 143a1 by counting the points of y² + y = x³ − x² − x − 2 modulo the prime p."""
    if p == 2:
        return 2 + 1 - (1 + sum(1 for x in range(2) for y in range(2) if (y * y + y - (x ** 3 - x * x - x - 2)) % 2 == 0))
    s = 0
    for x in range(p):
        v = (4 * x ** 3 - 4 * x * x - 4 * x - 7) % p
        if v:
            s += 1 if pow(v, (p - 1) // 2, p) == 1 else -1
    return -s


@lru_cache(maxsize=4)
def an_143a1(nmax: int = 6000) -> tuple:
    """a_1, …, a_nmax of the newform of 143a1 (a_p from point counts, multiplicativity, bad primes 11 and 13)."""
    sieve = bytearray([1]) * (nmax + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(nmax ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    an = [0, 1] + [0] * (nmax - 1)
    for p in range(2, nmax + 1):
        if sieve[p]:
            ap = ap_143a1(p)
            pk, prev, cur = p, 1, ap
            while pk <= nmax:
                an[pk] = cur
                prev, cur = cur, (ap * cur - (0 if p in (11, 13) else p) * prev)
                pk *= p
    for n in range(2, nmax + 1):
        if an[n] == 0 and not sieve[n]:
            m, p = n, next(p for p in range(2, n) if n % p == 0)
            pk = 1
            while m % p == 0:
                m //= p
                pk *= p
            if m > 1:
                an[n] = an[pk] * an[m]
    return tuple(an[1:])


def period_143a1(g, nmax: int = 6000) -> complex:
    """∫ 2πi f(z) dz of the newform of 143a1 from z₀ = (−d + i)/c to g z₀ = (a + i)/c, for g in Γ₀(143) with c ≠ 0."""
    (a, _), (c, d) = g
    if c == 0:
        return 0j
    anf = np.array(an_143a1(nmax), dtype=float) / np.arange(1, nmax + 1)
    narr = np.arange(1, nmax + 1)
    z0, z1 = complex(-d, 1) / c, complex(a, 1) / c
    return complex(np.sum(anf * (np.exp(2j * np.pi * narr * z1) - np.exp(2j * np.pi * narr * z0))))


def small_elements() -> list:
    """The twelve small elements of Γ₀(143) the lab integrated along."""
    out = [((1, 0), (143, 1)), ((12, 1), (143, 12)), ((84, -47), (143, -80)), ((1, 0), (286, 1)), ((72, 1), (143, 2)), ((191, 2), (286, 3))]
    for d in (4, 5, 6, 7, 9, 10):
        a = pow(d, -1, 143)
        out.append(((a, (a * d - 1) // 143), (143, d)))
    return out


def _klein_j(tau: complex) -> complex:
    import mpmath as mp
    with mp.workdps(30):
        return complex(1728 * mp.kleinj(mp.mpc(tau)))


@lru_cache(maxsize=1)
def lattice_143a1() -> dict:
    """Ω₁, Ω₂ of 143a1 from twelve small elements, with the worst residual and the lattice's j against the curve's."""
    small = small_elements()
    per = [period_143a1(g) for g in small]
    ell = [coordinates_of(g)["ell"] for g in small]
    n = len(small)
    i_, j_ = max(((i, j) for i in range(n) for j in range(i + 1, n)),
                 key=lambda ij: abs(ell[ij[0]][0] * ell[ij[1]][1] - ell[ij[0]][1] * ell[ij[1]][0]))
    det = ell[i_][0] * ell[j_][1] - ell[i_][1] * ell[j_][0]
    if det == 0:
        raise ArithmeticError("the small elements do not span the block")
    om1 = (per[i_] * ell[j_][1] - per[j_] * ell[i_][1]) / det
    om2 = (-per[i_] * ell[j_][0] + per[j_] * ell[i_][0]) / det
    err = max(abs(per[k] - (ell[k][0] * om1 + ell[k][1] * om2)) for k in range(n))
    tau = om2 / om1 if (om2 / om1).imag > 0 else om1 / om2
    jl = _klein_j(tau)
    return {"Omega1": om1, "Omega2": om2, "residual": err, "j_lattice": jl, "j_curve": float(CURVE_J),
            "j_error": abs(jl - float(CURVE_J)) / abs(float(CURVE_J)), "pair": (i_, j_), "det": det,
            "coordinates": ell, "periods": per}


def period_of(cls) -> complex:
    """∫ 2πi f dz of the newform of 143a1 along a class of ℤ²⁹: n₁Ω₁ + n₂Ω₂."""
    lat = lattice_143a1()
    n1, n2 = coordinates(cls)["ell"]
    return n1 * lat["Omega1"] + n2 * lat["Omega2"]


def stored_periods_in_lattice() -> dict:
    """The second route to the lattice: ``mtft.periods``' frozen periods of 143a1 (PARI/GP) written in Ω₁, Ω₂.

    The stored entries are ∫ f dz over a basis of H₁(X₀(143), ℤ); times 2πi they must be integer
    combinations of Ω₁ and Ω₂, and together generate the whole lattice (index 1).
    """
    import mtft.periods.core as PC
    rec = PC.period_record()
    row = rec["rows"].index("f1_143a")
    O = PC.omega_cusp()
    lat = lattice_143a1()
    om1, om2 = lat["Omega1"], lat["Omega2"]
    M = np.array([[om1.real, om2.real], [om1.imag, om2.imag]])
    vals = np.array([complex(O[row, j]) * 2j * math.pi for j in range(O.cols)])
    sol = np.linalg.solve(M, np.vstack([vals.real, vals.imag]))
    near = np.rint(sol)
    ints = [(int(near[0, j]), int(near[1, j])) for j in range(O.cols)]
    index = 0
    for i in range(len(ints)):
        for j in range(i + 1, len(ints)):
            index = math.gcd(index, ints[i][0] * ints[j][1] - ints[i][1] * ints[j][0])
    return {"worst": float(np.abs(sol - near).max()), "integers": ints, "index": index, "count": len(ints)}


def tau_143a1() -> complex:
    """The lattice's modulus in the upper half-plane."""
    lat = lattice_143a1()
    t = lat["Omega2"] / lat["Omega1"]
    return t if t.imag > 0 else 1 / t
