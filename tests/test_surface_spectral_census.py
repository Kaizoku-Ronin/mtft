"""SPEC-02: Hecke/AL census on Gamma_0(11). Where each classifier is confident they must agree;
modes that are neither Y-stable nor Eisenstein-like at this mesh are an ambiguity band, reported not asserted."""
from mtft.surface import spectral as SP


def test_n11_census_two_classifiers_agree_low_spectrum():
    r = SP.maass_candidates(11, h=0.2, nx=8, k=12)
    rows = [x for x in SP.arithmetic_census(11, h=0.2, nx=8, k=12, samples_per_face=60, stable_from=r["stable"]) if x["lambda"] <= 7.0]
    stable = [x for x in rows if x["Y_stable"]]
    eis = [x for x in rows if x["eisenstein_like"]]
    assert len(stable) >= 2 and len(eis) >= 2
    for x in stable:
        assert not x["eisenstein_like"] and abs(abs(x["W11"]) - 1) < 0.1   # cusp form: not Eisenstein, clean AL parity
    for x in eis:
        assert not x["Y_stable"]                                          # pseudo-mode: moves with Y
    assert all(abs(x["a2"]) <= 2.2 for x in rows)
