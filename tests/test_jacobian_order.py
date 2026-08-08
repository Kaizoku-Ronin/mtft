"""tests/test_jacobian_order.py — gates for mtft.crypto.jacobian_order
(v0.13.0 ancestry+crypto wave).

K3 audit tier (disclosed completion BU-2): the developer shipped the
module + 93-prime CSV without a gate; this file encodes the session's
independent verification as permanent regression tests.

  A_1 orders  — all 93 primes vs direct Legendre-symbol point counting
                on 143a1 (y^2 + y = x^3 - x^2 - x - 2)            [EXACT]
  A_2/A_3     — charpoly route on the repo's own integer Hecke
                matrices: |A_i(F_q)| = F_i,q(q+1) at q in {2,3,5,7} [EXACT]
  field polys — discriminants 1957 / 194616205 match the repo's T2
                charpoly factors (same quartic/sextic fields)     [EXACT]
  integrity   — Hasse-Weil bounds, embedding-degree minimality,
                candidates()/total_jacobian_order consistency
"""
import math

import pytest

from mtft.crypto import jacobian_order as jo


@pytest.fixture(scope="module")
def table():
    t = jo.JacobianOrderTable()
    assert len(t._records) == 279, "CSV should hold 93 primes x 3 orbits"
    return t


def _count_ec_points(q):
    """#E(F_q) for E: y^2 + y = x^3 - x^2 - x - 2, exact."""
    if q == 2:
        c = 1
        for x in (0, 1):
            rhs = (x**3 - x*x - x - 2) % 2
            for y in (0, 1):
                c += ((y*y + y - rhs) % 2 == 0)
        return c
    c = 1  # the point at infinity
    for x in range(q):
        rhs = (x**3 - x*x - x - 2) % q
        d = (4 * rhs + 1) % q              # (2y+1)^2 = 4(x^3-x^2-x-2)+1
        if d == 0:
            c += 1
        else:
            c += 1 + (1 if pow(d, (q - 1) // 2, q) == 1 else -1)
    return c


def test_table_integrity(table):
    assert len(table.primes()) == 93
    assert {r.dim for r in table.by_orbit(1)} == {1}
    assert {r.dim for r in table.by_orbit(2)} == {4}
    assert {r.dim for r in table.by_orbit(3)} == {6}
    for r in table._records:
        lo, hi = jo.hasse_weil_bounds(r.q, r.dim)
        assert lo <= r.order <= hi, f"Hasse-Weil violated at q={r.q} orbit {r.orbit}"
        assert 0.0 < r.security_ratio() <= 2.0


def test_a1_orders_exact(table):
    for q in table.primes():
        rec = table.get(q, 1)
        assert rec.order == _count_ec_points(q), f"A_1 mismatch at q={q}"


def test_embedding_degrees(table):
    for r in table._records:
        if r.embedding_degree <= 0:
            continue
        l, q, e = r.largest_prime, r.q, r.embedding_degree
        assert pow(q, e, l) == 1, f"q^e != 1 mod l at q={q} orbit {r.orbit}"
        # minimality: e minimal iff q^(e/p) != 1 for every prime p | e
        m, p = e, 2
        primes = set()
        while p * p <= m:
            while m % p == 0:
                primes.add(p); m //= p
            p += 1
        if m > 1:
            primes.add(m)
        for p in primes:
            assert pow(q, e // p, l) != 1, f"embedding degree not minimal at q={q}"


def test_candidates_and_totals(table):
    cands = table.candidates(orbit=3, min_ratio=0.70, min_bits=30)
    assert all(c.orbit == 3 and c.security_ratio() >= 0.70
               and c.largest_prime_bits >= 30 for c in cands)
    q = table.primes()[10]
    prod = math.prod(table.get(q, i).order for i in (1, 2, 3))
    assert jo.total_jacobian_order(q, table) == prod
    assert jo.total_jacobian_order(999983, table) is None
    assert jo.LEVEL == 143 and jo.ORBIT_DIMS == (1, 4, 6)


def _poly_from_leading(coeffs, x):
    sp = __import__("sympy")
    n = len(coeffs) - 1
    return sum(c * x ** (n - i) for i, c in enumerate(coeffs))


def test_field_polys_and_al():
    sp = pytest.importorskip("sympy")
    x = sp.Symbol("x")
    d2 = int(sp.discriminant(_poly_from_leading(jo.ORBIT_FIELD_POLYS[2], x), x))
    d3 = int(sp.discriminant(_poly_from_leading(jo.ORBIT_FIELD_POLYS[3], x), x))
    assert d2 == 1957            # == disc of the repo T2 charpoly factor F2
    assert d3 == 194616205       # == disc of F3  (same fields, EXACT)
    assert jo.ORBIT_FIELD_POLYS[1] is None
    assert jo.ATKIN_LEHNER == {1: (1, 1), 2: (-1, 1), 3: (1, -1)}


def test_a23_charpoly_spot():
    """|A_i(F_q)| = F_i,q(q+1) from the repo's integer Hecke matrices."""
    sp = pytest.importorskip("sympy")
    import os
    import sys
    import numpy as np
    studies = os.path.join(os.path.dirname(__file__), "..", "studies")
    sys.path.insert(0, os.path.abspath(studies))
    from x0143_particle_box import P1, tessellation, ModularSymbols
    p1 = P1(143); tessellation(143); ms = ModularSymbols(p1)
    x = sp.Symbol("x")
    table = jo.JacobianOrderTable()
    for q in (2, 3, 5, 7):
        M = ms.hecke_on_quotient(q)
        Mint = sp.Matrix(np.rint(np.array(M, dtype=float)).astype(int).tolist())
        for f, _m in sp.factor_list(sp.Poly(Mint.charpoly(x).as_expr(), x))[1]:
            deg = sp.Poly(f, x).degree()
            if deg == 4:
                assert int(f.subs(x, q + 1)) == table.get(q, 2).order
            elif deg == 6:
                assert int(f.subs(x, q + 1)) == table.get(q, 3).order
