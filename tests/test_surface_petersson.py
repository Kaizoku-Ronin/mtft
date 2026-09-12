import pytest

from mtft.surface import petersson as PT


@pytest.mark.slow
def test_bilinear_relations_vs_quadrature():
    r = PT.petersson_gate(Y0=2.0, h=0.3, nx=5)
    assert all(abs(x - 1) < 0.03 for x in r["diag_ratio"]) and r["offdiag_max_dev"] < 0.05 and r["cross_sector_vanishing"]
