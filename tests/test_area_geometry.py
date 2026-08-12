"""Tests for the arithmetic area geometry promotion and CC-04/CC-03."""
from itertools import combinations

import mpmath as mp
import pytest

from mtft.curvature import (brioschi, cold_amplitude, finite_atom_curvature,
                            rigidity_class, CURVATURE)
from mtft.moments import weight_susceptibility


def w(n):
    return sum(mp.log(d) / d for d in range(2, n + 1) if n % d == 0)


def X(n):
    return (mp.log(n), w(n)) if n > 1 else (mp.mpf(0), mp.mpf(0))


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def geom(beta, atoms):
    b = mp.mpf(beta)
    P = [X(n) for n in atoms]
    q = [mp.e ** (-b * x[0]) for x in P]
    Z = sum(q)
    p = [x / Z for x in q]
    m = len(P)
    mu = (sum(p[i] * P[i][0] for i in range(m)),
          sum(p[i] * P[i][1] for i in range(m)))
    C = [(x[0] - mu[0], x[1] - mu[1]) for x in P]
    E = sum(p[i] * C[i][0] ** 2 for i in range(m))
    F = sum(p[i] * C[i][0] * C[i][1] for i in range(m))
    G = sum(p[i] * C[i][1] ** 2 for i in range(m))
    tri = sum(p[i] * p[j] * p[k]
              * cross((P[j][0] - P[i][0], P[j][1] - P[i][1]),
                      (P[k][0] - P[i][0], P[k][1] - P[i][1])) ** 2
              for i, j, k in combinations(range(m), 3))
    return E * G - F * F, tri


# ── the two identities ──────────────────────────────────────────────

@pytest.mark.parametrize("atoms", [(1, 2, 3, 4, 5, 6), (1, 2, 3, 5, 7, 10)])
def test_fisher_volume_is_squared_arithmetic_area(atoms):
    mp.mp.dps = 80
    for beta in ("2.7", "5.3"):
        d, t = geom(beta, atoms)
        assert abs(d - t) < mp.mpf("1e-60")


# ── the rigidity theorem ────────────────────────────────────────────

@pytest.mark.parametrize("atoms", [(1, 2, 3), (1, 2, 3, 4),
                                   (1, 2, 4, 8), (1, 2, 4, 16)])
def test_rigidity_class_is_one_quarter(atoms):
    """One off-line point => K = 1/4 identically. {1,2,4,8} and
    {1,2,4,16} are out-of-sample: they did not motivate the theorem."""
    mp.mp.dps = 60
    assert rigidity_class([X(n) for n in atoms]) == 1 or len(atoms) == 3
    for beta in (7, 19, 41):
        assert abs(finite_atom_curvature(beta, atoms)
                   - mp.mpf(1) / 4) < mp.mpf("1e-60")


def test_rigidity_hypothesis_is_sharp():
    mp.mp.dps = 60
    with pytest.raises(ZeroDivisionError):
        finite_atom_curvature(6, (1, 2, 4))          # zero off-line
    for atoms in ((1, 2, 4, 8, 16), (1, 2, 3, 4, 5)):  # two off-line
        assert rigidity_class([X(n) for n in atoms]) == 2
        assert abs(finite_atom_curvature(17, atoms)
                   - mp.mpf(1) / 4) > 1


def test_dyadic_collinearity_is_exact():
    mp.mp.dps = 60
    assert abs(X(4)[0] - 2 * X(2)[0]) < mp.mpf("1e-50")
    assert abs(X(4)[1] - 2 * X(2)[1]) < mp.mpf("1e-50")


# ── the prime boundary, with the n = 1 case ─────────────────────────

def test_prime_boundary_including_n_equals_one():
    mp.mp.dps = 50
    for n in range(1, 200):
        lo = mp.log(n) / n if n > 1 else mp.mpf(0)
        assert w(n) >= lo - mp.mpf("1e-40")
        eq = abs(w(n) - lo) < mp.mpf("1e-40")
        isp = n > 1 and all(n % d for d in range(2, int(n ** .5) + 1))
        assert eq == (isp or n == 1)


# ── CC-04: the corrected cold amplitude ─────────────────────────────

def test_cold_amplitude_closed_form():
    mp.mp.dps = 60
    c = cold_amplitude()
    assert mp.nstr(c, 32) == CURVATURE["cold_amplitude"]
    k = finite_atom_curvature(200, (1, 2, 3, 5))
    assert abs(-k / (mp.mpf(6) / 5) ** 200 - c) < mp.mpf("1e-16")


def test_cc02_retraction_is_preserved_and_diagnosed():
    """The retracted value must still be reproducible from the six-atom
    model, which is what proves it was atom-6 contamination."""
    mp.mp.dps = 60
    old = mp.mpf(CURVATURE["cold_amplitude_RETRACTED_CC04"])
    k6 = finite_atom_curvature(200, (1, 2, 3, 4, 5, 6))
    assert abs(-k6 / (mp.mpf(6) / 5) ** 200 - old) < mp.mpf("1e-22")
    assert abs(old - cold_amplitude()) < 3 * (mp.mpf(5) / 6) ** 200


# ── CC-03: precision may only be raised ─────────────────────────────

def test_cc03_low_dps_override_cannot_corrupt():
    """This exact call returned ~-2e+19 before the adaptive fix."""
    mp.mp.dps = 60
    assert abs(finite_atom_curvature(80, (1, 2, 4, 8), dps=90)
               - mp.mpf(1) / 4) < mp.mpf("1e-60")


# ── the arithmetic area zeta ────────────────────────────────────────

def test_area_zeta_wall_is_fifth_order():
    mp.mp.dps = 40
    s = 1 + mp.mpf("1e-4")
    E = mp.diff(lambda x: mp.log(mp.zeta(x)), s, 2)
    G = weight_susceptibility(s)
    F = mp.diff(mp.zeta, s + 1, 2)
    A = mp.zeta(s) ** 3 * (E * G - F ** 2)
    assert abs(mp.mpf("1e-4") ** 5 * A / weight_susceptibility(1) - 1) \
        < mp.mpf("2e-4")
