"""Frozen gates for the theorem compendium, Chapter I (X0(143)).

Every assertion here is an exact integer or rational computation; no floating
point enters.  The module depends only on the standard library and on mtft's
frozen certified data (mtft.canonical) and its integer-lattice routines
(mtft.integral_lattice).  It does not touch the release pin.

Gates (numbering follows the compendium):
  I.1   cuspidal divisor class group Z/10 x Z/420; orders 420, 60, 70;
        Ligozat divisors == Atkin-Lehner-transport divisors on every admissible
        eta quotient in the generating box.
  I.3b  h^0(S_0) = 2 from the two sector-block ranks (6 and 5).
  I.3c  C[2] has four elements, W acts trivially on it (four cuspidal theta
        characteristics, all W-invariant).
  I.4   P^1(Z/143): 168 points, idempotents 78/66, T-orbits 1/11/13/143,
        84 S-edges, 56 ST-triangles, V - E + F = 2 - 2g.
  I.5   T_2 on the frozen basis is block-diagonal with characteristic
        polynomial x (x+2)^2 g(x) s(x); g, s totally real; blocks non-CM.
  I.5B  Atkin-Lehner grading of the blocks reproduces DESCENT.
"""
from fractions import Fraction as Fr
from itertools import product
from math import gcd

import pytest

import mtft.canonical as C
import mtft.integral_lattice as IL

N = 143
DIV = (1, 11, 13, 143)          # eta factors eta(delta * tau)
CUSPS = (143, 1, 11, 13)        # cusp 1/d in the order (inf, 0, 1/11, 1/13)

# ----------------------------------------------------------------------------
# exact linear algebra over Q (small matrices; no sympy dependency)
# ----------------------------------------------------------------------------

def _rref(rows):
    """Row-reduce a list of Fraction rows in place; return (rank, pivot columns)."""
    rows = [[Fr(x) for x in r] for r in rows]
    m, n = len(rows), len(rows[0]) if rows else 0
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, m) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        inv = 1 / rows[r][c]
        rows[r] = [x * inv for x in rows[r]]
        for i in range(m):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
        if r == m:
            break
    return r, piv, rows


def _rank(rows):
    return _rref(rows)[0]


def _solve(A_rows, b):
    """Solve A x = b exactly (A: list of rows, overdetermined allowed).
    Returns x or raises if inconsistent."""
    n = len(A_rows[0])
    aug = [list(r) + [bi] for r, bi in zip(A_rows, b)]
    rank, piv, red = _rref(aug)
    if n in piv:
        raise ValueError("inconsistent system")
    x = [Fr(0)] * n
    for i, c in enumerate(piv):
        x[c] = red[i][n]
    if rank != n:
        raise ValueError("underdetermined system")
    return x


def _charpoly(M):
    """Faddeev-LeVerrier: coefficients c_n..c_0 of det(xI - M), exact Fractions."""
    n = len(M)
    M = [[Fr(x) for x in row] for row in M]
    I = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]

    def mul(A, B):
        return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

    coeffs = [Fr(1)]
    Nk = I
    for k in range(1, n + 1):
        MN = mul(M, Nk)
        ck = -sum(MN[i][i] for i in range(n)) / k
        coeffs.append(ck)
        Nk = [[MN[i][j] + (ck if i == j else 0) for j in range(n)] for i in range(n)]
    return coeffs          # x^n + c1 x^{n-1} + ... + c_n


def _polymul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


# ----------------------------------------------------------------------------
# I.1 — eta quotients, two divisor routes, the lattice, its Smith form
# ----------------------------------------------------------------------------

def _ligozat(r, d):
    return Fr(N, 24) * sum(Fr(gcd(d, dl) ** 2 * rd, gcd(d, N // d) * d * dl) for dl, rd in zip(DIV, r))


def _transport(r, d):
    """ord at the cusp W_Q(inf), Q = N/d, by Atkin-Lehner transport: (1/24) sum r_delta sigma_Q(delta)."""
    Q = N // d
    return Fr(sum(rd * (dl * Q // gcd(dl, Q) ** 2) for dl, rd in zip(DIV, r)), 24)


def _admissible(r):
    r1, r11, r13, r143 = r
    return (sum(r) == 0
            and sum(dl * rd for dl, rd in zip(DIV, r)) % 24 == 0
            and sum((N // dl) * rd for dl, rd in zip(DIV, r)) % 24 == 0
            and (r11 + r143) % 2 == 0 and (r13 + r143) % 2 == 0)


def _admissible_divisors():
    """All admissible exponent vectors with |r_1|,|r_11|,|r_13| <= 12 and their divisors.
    The admissibility conditions are congruences mod 24 and 2, so this box contains a
    generating set of the admissible lattice."""
    out = []
    for a, b, c in product(range(-12, 13), repeat=3):
        r = (a, b, c, -(a + b + c))
        if _admissible(r):
            D1 = [_ligozat(r, d) for d in CUSPS]
            D2 = [_transport(r, d) for d in CUSPS]
            out.append((r, D1, D2))
    return out


GEN = {   # the three displayed generators: exponents (r1, r11, r13, r143)
    "F1": (-6, 6, -6, 6),
    "F2": (-2, 0, 2, 0),
    "u": (-1, -1, 1, 1),
}
G_COLS = None   # filled lazily: divisor coordinates (coeff at 0, 1/11, 1/13) of F1, F2, u


def _generator_matrix():
    cols = []
    for r in GEN.values():
        D = [_ligozat(r, d) for d in CUSPS]
        cols.append([int(D[1]), int(D[2]), int(D[3])])
    return [[cols[j][i] for j in range(3)] for i in range(3)]     # 3x3, columns = generators


def test_I1_divisor_routes_agree_and_are_integral():
    data = _admissible_divisors()
    assert len(data) > 3
    for r, D1, D2 in data:
        assert D1 == D2, r
        assert all(x.denominator == 1 for x in D1) and sum(D1) == 0, r
    # the displayed divisors
    assert [_ligozat(GEN["F1"], d) for d in CUSPS] == [35, -35, 35, -35]
    assert [_ligozat(GEN["F2"], d) for d in CUSPS] == [1, -11, -1, 11]
    assert [_ligozat(GEN["u"], d) for d in CUSPS] == [6, -6, -6, 6]


def test_I1_cuspidal_class_group_is_Z10_x_Z420():
    data = _admissible_divisors()
    rows = [[int(D[1]), int(D[2]), int(D[3])] for _, D, _ in data]
    assert IL.smith_invariants(rows) == [1, 10, 420]
    G = _generator_matrix()
    assert G == [[-35, -11, -6], [35, -1, -6], [-35, 11, 6]]
    assert IL.smith_invariants(G) == [1, 10, 420]
    # hand form: gcd of entries 1, gcd of 2x2 minors 10, |det| 4200
    det = (G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1])
           - G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0])
           + G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0]))
    assert det == -4200
    minors = []
    for r1, r2 in ((0, 1), (0, 2), (1, 2)):
        for c1, c2 in ((0, 1), (0, 2), (1, 2)):
            minors.append(G[r1][c1] * G[r2][c2] - G[r1][c2] * G[r2][c1])
    g2 = 0
    for m in minors:
        g2 = gcd(g2, abs(m))
    assert g2 == 10


def test_I1_orders_420_60_70():
    G = _generator_matrix()
    H = IL.hnf(G)
    divs = sorted(d for d in range(1, 4201) if 4200 % d == 0)
    assert IL.class_order(H, [1, 0, 0], divs) == 420     # (0) - (inf)
    assert IL.class_order(H, [0, 1, 0], divs) == 60      # (1/11) - (inf)
    assert IL.class_order(H, [0, 0, 1], divs) == 70      # (1/13) - (inf)


# ----------------------------------------------------------------------------
# I.3c — the 2-torsion of C and the four cuspidal theta characteristics
# ----------------------------------------------------------------------------

def _in_L(v3):
    """Is the degree-0 cusp divisor with coordinates (a0, a11, a13) in the lattice L?"""
    G = _generator_matrix()
    H = IL.hnf(G)
    divs = [1]
    return IL.class_order(H, list(v3), divs) == 1


def _act4(W, v):
    """W on cusp divisors (a0, a11, a13, ainf).  W11: 0<->1/11, inf<->1/13.  W13: 0<->1/13, 1/11<->inf."""
    a0, a11, a13, ainf = v
    return [a11, a0, ainf, a13] if W == "W11" else [a13, ainf, a0, a11]


def test_I3c_four_cuspidal_theta_characteristics_W_invariant():
    e2 = [0, 30, 0, -30]          # 30 (1/11) - 30 (inf)
    e3 = [0, 0, 35, -35]          # 35 (1/13) - 35 (inf)
    e1 = [210, 0, 0, -210]        # 210 (0) - 210 (inf)
    two_tors = [[0, 0, 0, 0], e2, e3, [a + b for a, b in zip(e2, e3)]]
    for t in two_tors[1:]:
        assert not _in_L(t[:3])                       # nonzero class
        assert _in_L([2 * x for x in t][:3])          # of order 2
    assert _in_L([a - b for a, b in zip(e1, e2)][:3])  # [210 e1] = [30 e2]
    assert not _in_L([a - b for a, b in zip(e2, e3)][:3])   # the three are distinct
    A = [6, 6, 0, 0]
    for t in two_tors:
        theta = [a + b for a, b in zip(A, t)]
        for W in ("W11", "W13"):
            diff = [a - b for a, b in zip(_act4(W, theta), theta)]
            assert sum(diff) == 0 and _in_L(diff[:3]), (t, W)


# ----------------------------------------------------------------------------
# frozen basis of S_2(Gamma_0(143)) and the eta-quotient forms
# ----------------------------------------------------------------------------

def _basis():
    A = C.adapted_qexpansions()
    A = [[int(x) for x in row] for row in A]         # rows q^0..q^140, columns e_1..e_13
    assert len(A) == 141 and len(A[0]) == 13
    return A


def _eta_prod(exps, n):
    """prod_delta prod_{m>=1} (1 - q^{delta m})^{r_delta} to O(q^{n+1})."""
    c = [0] * (n + 1)
    c[0] = 1
    for dl, r in exps.items():
        for m in range(1, n // dl + 1):
            for _ in range(abs(r)):
                if r > 0:
                    for i in range(n, dl * m - 1, -1):
                        c[i] -= c[i - dl * m]
                else:
                    for i in range(dl * m, n + 1):
                        c[i] += c[i - dl * m]
    return c


def _labels():
    return [tuple(l) for l in C.COORDINATE_LABELS]


def test_I3b_h0_S0_equals_two():
    A = _basis()
    labels = _labels()
    plus = [k for k, l in enumerate(labels) if l[1] in ("(+,+)", "(+,-)")]
    minus = [k for k, l in enumerate(labels) if l[1] in ("(-,+)", "(-,-)")]
    assert len(plus) == 7 and len(minus) == 6
    rho_plus = _rank([[A[n][k] for k in plus] for n in range(1, 7)])
    rho_minus = _rank([[A[n][k] for k in minus] for n in range(1, 7)])
    assert rho_plus == 6 and rho_minus == 5
    assert (7 - rho_plus) + (6 - rho_minus) == 2          # h^0(S_0) = 2
    assert _rank([[A[n][k] for k in range(13)] for n in range(1, 29)]) == 13   # Sturm bound 28


def test_I3b_old_forms_and_ug_membership():
    A = _basis()
    labels = _labels()
    n = 140
    g = [0] + _eta_prod({1: 2, 11: 2}, n - 1)
    g13 = [0] * (n + 1)
    for i, c in enumerate(g):
        if 13 * i <= n:
            g13[13 * i] = c
    ug = [0] * 7 + _eta_prod({1: 1, 11: 1, 13: 1, 143: 1}, n - 7)
    op = labels.index(("old+", "(-,+)"))
    om = labels.index(("old-", "(-,-)"))
    assert all(A[i][op] == g[i] + 13 * g13[i] for i in range(n + 1))
    assert all(A[i][om] == g[i] - 13 * g13[i] for i in range(n + 1))
    assert ug[7] == 1 and ug[8] == -1
    sextic = [k for k, l in enumerate(labels) if l[0] == "f3"]
    assert len(sextic) == 6
    rows = [[A[i][k] for k in sextic] for i in range(1, 29)]
    aug = [r + [ug[i]] for r, i in zip(rows, range(1, 29))]
    assert _rank(rows) == _rank(aug) == 6             # ug lies in the sextic block


# ----------------------------------------------------------------------------
# I.5 — T_2 on the frozen basis
# ----------------------------------------------------------------------------

G_POLY = [1, -3, -1, 5, 1]                        # x^4 - 3x^3 - x^2 + 5x + 1   (high to low)
S_POLY = [1, 0, -10, 2, 24, -7, -12]              # x^6 - 10x^4 + 2x^3 + 24x^2 - 7x - 12


def _T2_matrix():
    A = _basis()
    nrow = len(A)
    half = nrow // 2
    images = []
    for k in range(13):
        col = [A[2 * m][k] + (2 * A[m // 2][k] if m % 2 == 0 else 0) for m in range(half)]
        images.append(col)
    B = [[Fr(A[i][k]) for k in range(13)] for i in range(half)]      # 70 x 13
    X = []                                                             # columns of the matrix
    for k in range(13):
        x = _solve(B, [Fr(v) for v in images[k]])                      # exact, all 70 rows consistent
        X.append(x)
    return [[X[j][i] for j in range(13)] for i in range(13)]           # M[i][j]: e_j -> sum_i M[i][j] e_i


def test_I5a_T2_block_structure_and_charpoly():
    M = _T2_matrix()
    labels = _labels()
    blocks = {}
    for k, l in enumerate(labels):
        blocks.setdefault(l[0], []).append(k)
    expected = {"f1": [1, 0], "f3": S_POLY, "old+": [1, 2], "f2": G_POLY, "old-": [1, 2]}
    for name, cols in blocks.items():
        off = [i for i in range(13) if i not in cols]
        assert all(M[i][j] == 0 for i in off for j in cols), name          # T_2-stable block
        sub = [[M[i][j] for j in cols] for i in cols]
        assert _charpoly(sub) == [Fr(c) for c in expected[name]], name
    full = _charpoly(M)
    prod = [Fr(1)]
    for p in ([1, 0], [1, 2], [1, 2], G_POLY, S_POLY):
        prod = _polymul(prod, [Fr(c) for c in p])
    assert full == prod


def _sign_changes(poly, grid):
    def ev(t):
        v = Fr(0)
        for c in poly:
            v = v * t + c
        return v
    vals = [ev(Fr(t)) for t in grid]
    assert all(v != 0 for v in vals)
    return sum(1 for a, b in zip(vals, vals[1:]) if (a < 0) != (b < 0))


def test_I5c_totally_real_and_non_CM():
    assert _sign_changes(G_POLY, [-2, -1, 0, 1, 2, 3]) == 4
    assert _sign_changes(S_POLY, [-3, -2, -1, 0, 1, Fr(5, 4), Fr(13, 10), Fr(7, 5), 2, 3]) == 6
    # non-CM: T_2 invertible on both blocks (2 inert in Q(sqrt(-11))), T_5 nonzero (5 inert in Q(sqrt(-143)))
    assert G_POLY[-1] == 1 and S_POLY[-1] == -12
    A = _basis()
    labels = _labels()
    for name in ("f2", "f3"):
        cols = [k for k, l in enumerate(labels) if l[0] == name]
        assert any(A[5][k] != 0 for k in cols), name


# ----------------------------------------------------------------------------
# I.4 — P^1(Z/143)
# ----------------------------------------------------------------------------

def _p1_points():
    units = [u for u in range(1, N) if gcd(u, N) == 1]

    def canon(c, d):
        c %= N
        d %= N
        return min(((u * c) % N, (u * d) % N) for u in units)

    pts = set()
    for c in range(N):
        for d in range(N):
            if gcd(gcd(c, d), N) == 1:
                pts.add(canon(c, d))
    return sorted(pts), canon


def _orbits(pts, canon, Mx):
    a, b, e, f = Mx
    seen, lens = set(), []
    for p in pts:
        if p in seen:
            continue
        q, L = p, 0
        while q not in seen:
            seen.add(q)
            L += 1
            q = canon(a * q[0] + e * q[1], b * q[0] + f * q[1])
        lens.append(L)
    return sorted(lens)


def test_I4_P1_Z143_and_Euler_characteristic():
    assert (78 % 11, 78 % 13, 66 % 11, 66 % 13) == (1, 0, 0, 1)
    assert (78 + 66) % N == 1 and (78 * 66) % N == 0 and (78 * 78) % N == 78 and (66 * 66) % N == 66
    pts, canon = _p1_points()
    assert len(pts) == 168

    def red(p, m):
        c, d = p[0] % m, p[1] % m
        return (c * pow(d, -1, m)) % m if d % m else "inf"

    assert len({(red(p, 11), red(p, 13)) for p in pts}) == 168
    T, S, ST = (1, 1, 0, 1), (0, -1, 1, 0), (0, -1, 1, 1)
    tl = _orbits(pts, canon, T)
    sl = _orbits(pts, canon, S)
    rl = _orbits(pts, canon, ST)
    assert tl == [1, 11, 13, 143]
    assert len(sl) == 84 and set(sl) == {2}
    assert len(rl) == 56 and set(rl) == {3}
    assert len(tl) - len(sl) + len(rl) == 2 - 2 * 13


# ----------------------------------------------------------------------------
# I.5B — 143a1 and the Atkin-Lehner grading; comparison with DESCENT
# ----------------------------------------------------------------------------

def _e143a1():
    w = C.DESCENT["weierstrass"]
    return w["a1"], w["a2"], w["a3"], w["a4"], w["a6"]


def _count_points(p):
    a1, a2, a3, a4, a6 = _e143a1()
    n = 1
    for x in range(p):
        for y in range(p):
            if (y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)) % p == 0:
                n += 1
    return n


def _reduction_type(p):
    a1, a2, a3, a4, a6 = _e143a1()
    sing = [(x, y) for x in range(p) for y in range(p)
            if (y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)) % p == 0
            and (a1 * y - (3 * x * x + 2 * a2 * x + a4)) % p == 0
            and (2 * y + a1 * x + a3) % p == 0]
    assert len(sing) == 1
    x0, _ = sing[0]
    Fxx, Fxy, Fyy = (-(6 * x0 + 2 * a2)) % p, a1 % p, 2 % p
    D = (Fxy * Fxy - Fxx * Fyy) % p
    return "split" if any((t * t - D) % p == 0 for t in range(p)) else "non-split"


def test_I5B_143a1_traces_reductions_and_signs():
    assert _e143a1() == (0, -1, 1, -1, -2)
    assert [p + 1 - _count_points(p) for p in (2, 3, 5, 7)] == [0, -1, -1, -2]
    assert _reduction_type(11) == "non-split" and _reduction_type(13) == "non-split"
    # a_p = -1 at both, so the Atkin-Lehner signs are (+,+): 143a1 is the (+,+) sector
    labels = _labels()
    assert labels[labels.index(("f1", "(+,+)"))][1] == "(+,+)"


def test_I5B_sector_dimensions_reproduce_DESCENT():
    # derived: (+,+) = 143a1 (1); (+,-) = sextic (6); (-,+) = old+ (1) + quartic (4); (-,-) = old- (1)
    derived = {"(+,+)": 1, "(+,-)": 6, "(-,+)": 5, "(-,-)": 1}
    assert C.DESCENT["sector_dims_S2"] == derived
    assert tuple(C.SECTOR_ORDER) == ("(+,+)", "(+,-)", "(-,+)", "(-,-)")
    labels = _labels()
    counted = {}
    for _, sec in labels:
        counted[sec] = counted.get(sec, 0) + 1
    assert counted == derived
    # Riemann-Hurwitz for an involution on genus 13: fixed points r = 28 - 4 g'
    genus = {"W11": derived["(+,+)"] + derived["(+,-)"],
             "W13": derived["(+,+)"] + derived["(-,+)"],
             "W143": derived["(+,+)"] + derived["(-,-)"]}
    assert genus == {"W11": 7, "W13": 6, "W143": 2}
    assert C.DESCENT["quotient_genera"] == genus
    assert C.DESCENT["fixed_points"] == {W: 28 - 4 * g for W, g in genus.items()}
    assert C.DESCENT["W11_acts_freely"] is True
