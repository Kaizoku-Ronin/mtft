"""GC-01: gauge couplings fix the volume and the Pati–Salam scale, not the shape."""
from mtft.matter import gauge_couplings as G


def test_ratio_and_scales():
    assert abs(G.ratio(91.1876) - 3.991) < 0.005 and abs(G.ratio(1.1236e6) - 3.678) < 0.005
    assert 4.5e13 < G.matching_scale(1.0) < 5.6e13 and 8e16 < G.alpha23_scale() < 1.2e17


def test_quantised_designs():
    assert [round(x, 3) for x in G.design_ratios()] == [0.111, 0.429, 1.0, 2.333, 9.0] and [round(x, 3) for x in G.single_point_ratios()] == [0.2, 1.0, 5.0]
    assert 7e11 < G.matching_scale(7 / 3) < 1.1e12 and G.matching_scale(9.0) is None and G.matching_scale(5.0) is None


def test_b_tau_tension_and_volume():
    assert 0.69 < G.b_tau_at(G.matching_scale(7 / 3)) < 0.73 and abs(G.b_tau_at(1.1236e6) - 1) < 0.01
    assert G.volume(18.0, 1.0, (1 / 3, 1 / 3)) == 18.0 - 1 / 9
