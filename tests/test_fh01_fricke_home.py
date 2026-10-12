"""FH-01: the top at the Fricke point."""
from mtft.matter import fricke_home as F


def test_arithmetic_and_parity():
    assert F.CURVE["w143"] == -F.CURVE["root_number"] == 1 and F.CURVE["rank"] == 1 and not F.is_two_torsion()
    assert 0.98 < F.shortest_dual() < 0.99
    assert 1.6 < F.separation(230) < 1.63 and 1.78 < F.separation(280) < 1.8
    assert 1e-5 < F.torus_suppression(280) < F.torus_suppression(230) < 5e-5


def test_concordance_and_vub():
    assert F.concordance() == (230.0, 280.0)
    sym, data = F.texture_zero_ratio(); assert 0.044 < sym < 0.046 and 0.089 < data < 0.091 and 1.95 < data / sym < 2.05
