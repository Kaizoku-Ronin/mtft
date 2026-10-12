"""BL-01: a B–L Wilson line on 143a1."""
import math
from mtft.matter import bl_wilson_line as B


def test_cancelling_characters_and_surviving_ones():
    assert abs(B.tau_143a1().imag - 1.02327) < 1e-4
    assert B.ratio(0.7, "1/2", (1, 0)) < 1e-6 and B.ratio(0.5, "tau/2", (0, 1)) < 1e-6 and B.ratio(0.5, "(1+tau)/2", (1, 0)) < 1e-6
    assert 0.92 < B.ratio(0.02, "1/2", (0, 1)) < 0.93 and 0.78 < B.ratio(0.02, "tau/2", (1, 0)) < 0.80


def test_lepton_point_misses_vcb():
    a = B.lepton_point("1/2", (0, 1)); b = B.lepton_point("tau/2", (1, 0))
    assert 0.70 < a["s_X"] < 0.73 and 0.13 < a["vcb_twisted"] < 0.135 and 0.49 < b["s_X"] < 0.51 and 0.18 < b["vcb_twisted"] < 0.19
    assert a["vcb_twisted"] / B.WINDOW[1] > 2.4 and b["vcb_twisted"] / B.WINDOW[0] > 6.1
    lo = B.vcb_point("1/2", (0, 1), B.WINDOW[1])["muon_factor"]; hi = B.vcb_point("tau/2", (1, 0), B.WINDOW[0])["muon_factor"]
    assert 6.5 < lo < 6.8 and 17.5 < hi < 18.5


def test_saddles_and_image_factor():
    pairs = [(s, c) for s in ((1, 0), (0, 1), (1, 1)) for c in ((1, 0), (0, 1), (1, 1)) if s != c]
    for p in (4, 6, 8):
        assert all(min(B.hessian_at(s, c, p)) < 0 for s, c in pairs)
    assert all(B.potential((math.pi * c[0], math.pi * c[1]), s) > B.potential((0.0, 0.0), s) for s, c in pairs)
    assert 3.0 < B.image_factor(0.85, "1/2") < 3.2 and 3.4 < B.image_factor(0.85, "tau/2") < 3.6
