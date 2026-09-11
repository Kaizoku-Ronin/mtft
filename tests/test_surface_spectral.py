"""SPEC-01: hyperbolic spectrum gates (N=11 fast; N=143 marked slow)."""
import pytest

from mtft.surface import spectral as SP


def test_n11_maass_candidates_and_pseudo_modes():
    r = SP.maass_candidates(11, h=0.16, nx=10, k=14)
    assert all(r["gates"].values())
    assert any(abs(x - 4.40) < 0.05 for x in r["stable"])       # first Y-stable mode of Gamma_0(11)
    assert r["pseudo_modes_Y0a"]                                  # Eisenstein pseudo-modes exist and move


@pytest.mark.slow
def test_n143_lambda1_above_kim_sarnak():
    r = SP.maass_candidates(143, h=0.22, nx=7, k=14)
    assert r["gates"]["lambda0_zero"] and r["gates"]["kim_sarnak_975_4096"]
    assert 0.3 < r["lambda_1"] < 0.5
