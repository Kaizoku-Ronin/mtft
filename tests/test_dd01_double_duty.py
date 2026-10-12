"""DD-01: the parity selection rule for the CKM matrix."""
from mtft.matter import double_duty as D


def test_parity_rule():
    assert D.parities() == {"top": 1, "x_jet": 1, "w_jet": -1}
    assert D.decoupled_family() == ["w_jet"] and D.allowed("top", "x_jet") and not D.allowed("top", "w_jet")


def test_numbers():
    assert abs(D.required_coefficient() - 0.401) < 0.002
    assert abs(D.separation(1 / 3) - 0.953) < 0.002 and abs(D.separation(0.55) - 1.094) < 0.002
    assert 17.5 < D.curve_area(1 / 3) < 19.0 and 27 < D.curve_area(0.55) < 29
