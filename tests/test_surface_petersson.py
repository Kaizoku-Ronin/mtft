import pytest

from mtft.surface import petersson as PT


@pytest.mark.slow
def test_bilinear_relations_vs_quadrature():
    r = PT.petersson_gate(Y0=2.0, h=0.3, nx=5)
    assert all(abs(x - 1) < 0.03 for x in r["diag_ratio"]) and r["offdiag_max_dev"] < 0.05 and r["prediction_cross_sector_zero"] and r["computed_cross_sector_max_correlation"] < 0.05  # v0.29.1 API (KK05); disclosed auditor fix
