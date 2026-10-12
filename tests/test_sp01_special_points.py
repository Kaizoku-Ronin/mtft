"""SP-01: canonical-metric distances on X_0(143) and the curve size fixed by V_cb."""
from mtft.matter import special_points as S


def test_canonical_metric_distances():
    A, Y = S.density_setup()
    assert abs(S.arakelov_density(2j, A, Y) / abs(__import__("cmath").exp(-4 * __import__("math").pi)) ** 2 / (4 * __import__("math").pi ** 2) - 0.02434) < 2e-5
    assert abs(S.path_length(0.18363, 0.03346, A, Y) - 0.1781) < 2e-3 and abs(S.path_length(0.0, 1 / 143 ** 0.5, A, Y) - 0.107) < 2e-3


def test_curve_size_and_cusp_mixing():
    AX = S.curve_area_for_Vcb(0.1781, 0.5629, 1 / 3); assert 17 < AX < 19.5 and 0.03 < S.cusp_mixing(AX, 1 / 3) < 0.06
    assert 30 < S.curve_area_for_Vcb(0.1886, 0.4943, 0.55) < 38 and 18 < S.curve_area_for_Vcb(S.CUSP_TO_CUSP, 0.0, 1 / 3) < 21
