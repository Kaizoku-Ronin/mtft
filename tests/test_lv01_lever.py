"""LV-01: the lever at the jet point."""
from mtft.matter import lever as L


def test_lever_formulas_against_exact_diagonalisation():
    for u, v in ((0.2, 0.0), (0.3, 0.1), (0.09, -0.135)):
        ex = L.down_only(u, v); fo = L.formulas(u, v)
        assert abs(ex[1] - fo[1]) < 0.005 and abs(ex[0] / fo[0] - 1) < 0.07 and abs(ex[3] / fo[3] - 1) < 0.07
    assert L.down_only(0.2, 0.0)[2] < 1e-9 and 0.98 < L.down_only(0.2, 0.0)[3] < 1.03      # pure lever: c = 1, massless first family
    assert L.down_only(0.0, 0.2)[1] < 1e-3                                                   # pure local arm: no down-sector V_ub


def test_what_the_data_ask():
    u, v, r = L.opposite_arms(); assert 0.089 < u < 0.091 and -0.136 < v < -0.134 and 0.041 < r < 0.044
    assert 0.19 < L.up_sector_c() < 0.21
    assert all(abs(f["gamma"] - L.GAMMA_MEASURED) > 40 for f in L.REAL_ARM_FITS)
    assert abs(L.TWIST_FITS[1]["beta"] + 90) < 2 and 0.25 < L.TWIST_FITS[0]["lever_share"] < 0.27
