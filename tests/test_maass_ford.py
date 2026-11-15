"""mtft.maass.ford: the Ford domain of Gamma_0*(143), its quadrature, norms, overlaps and participation ratios."""
import math

import numpy as np
import pytest

from mtft.maass import ford as F
from mtft.maass import spectrum as SP


def test_the_domain_at_143():
    D = F.domain(143)
    s = D.summary()
    assert (s["tops"], s["vertices"], s["distinct heights"]) == (17, 24, 16)                  # the record MP20
    assert abs(s["lowest vertex"][0] - 0.121503) < 1e-6 and abs(s["lowest vertex"][1] - 0.009851) < 1e-6
    assert abs(s["highest top"] - 1 / math.sqrt(143)) < 1e-15                                  # the Fricke circle
    assert abs(F.area(143) - 14 * math.pi) < 1e-12
    assert abs(D.area_by_quadrature(16) / (14 * math.pi) - 1) < 1e-9                          # the record: 2e-11
    assert abs(D.area_by_quadrature(4) / (14 * math.pi) - 1) > 1e-9                           # a coarse rule misses the gate
    assert D.intervals(0.5) == [(0.0, 0.5)] and D.intervals(0.009) == []


def test_other_levels():
    for N, low in ((1, math.sqrt(3) / 2), (11, 0.08963757247), (13, 0.06661733875)):
        D = F.domain(N)
        assert abs(D.lowest_vertex[1] - low) < 1e-10
        assert abs(D.area_by_quadrature(16) / D.area - 1) < 1e-12
    assert abs(F.area(1) - math.pi / 3) < 1e-15 and abs(F.area(11) - 2 * math.pi) < 1e-15
    with pytest.raises(ValueError):
        F.area(12)
    with pytest.raises(ValueError):
        F.circles(12)


def test_norms_and_orthogonality():
    a, b = SP.form("u.e.mm.1"), SP.form("u.e.mm.2")               # two forms of one sector: nothing but the
    G = F.gram([a, b])                                            # eigenvalue makes them orthogonal
    assert G[0, 0] > 0 and G[1, 1] > 0
    assert abs(G[0, 1]) / math.sqrt(G[0, 0] * G[1, 1]) < 1e-8     # the record: 2e-10
    assert abs(G[0, 1] - G[1, 0]) == 0


def test_participation_ratios_of_mp20():
    """MP20's seven ratios, I = area * int u^4 / (int u^2)^2, with the 400 frozen coefficients."""
    want = {"u.o.mm.1": 1.8905, "u.e.mm.1": 2.3288, "u.o.pp.1": 2.2821, "t.e.pp.1": 2.1920, "t.o.mm.1": 1.9035,
            "t.e.pm.1": 1.8884, "t.o.pm.1": 2.4722}
    fs = [SP.form(k) for k in want]
    I = F.participation_ratios(fs)
    for k, v in zip(want, I):
        assert abs(v - want[k]) < 1e-3, (k, v)
    assert np.all(I > 1.0)                                        # 1 is a constant function: a cusp form is above it


def test_products_and_their_symmetry():
    """An integral over F* means something only for the trivial character; other products are zero on the curve."""
    o, e1, e2, tw = SP.form("u.o.mm.1"), SP.form("u.e.mm.1"), SP.form("u.e.pp.1"), SP.form("t.o.mm.1")
    assert F.character([o, e1]) == (1, 1, 1, 0) and F.character([e1, e2]) == (0, -1, -1, 0) and F.character([tw]) == (1, -1, -1, 1)
    assert F.invariant([e1, e1]) and F.invariant([e1, e1, e2]) and F.invariant([o, o, e2]) and F.invariant([tw, tw, e2])
    assert not F.invariant([e1, e2, e2]) and not F.invariant([o, e1, e2]) and not F.invariant([tw, e1, e1])
    assert F.triple(e1, e2, e2) == 0.0 and F.triple(o, e1, e2) == 0.0           # forbidden: not integrated
    a, b = F.triple(e1, e1, e2), F.triple(o, o, e2)                             # allowed: integrated
    assert abs(a) > 1e-3 and abs(b) > 1e-3 and abs(a) < 1 and abs(b) < 1
    assert abs(F.triple(e2, e1, e1) - a) < 1e-12                                # the order of the forms does not matter
    G = F.gram([o, e1, e2])
    assert G[0, 1] == 0.0 and G[1, 2] == 0.0 and np.all(np.diag(G) > 0)         # different sectors: zero by symmetry
    # why the rule is needed: the half-domain rule applied to a product of non-trivial character gives a number
    # that is not zero and not an invariant of the curve
    raw = F.integrate([e1, e2], lambda U: U[0] * U[1] * U[1], degree=3)[0]
    assert abs(raw) > 1e-6
