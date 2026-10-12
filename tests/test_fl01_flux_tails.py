"""FL-01: which tail sets which angle, and the flux bound."""
from mtft.matter import flux_tails as F


def test_required_tails():
    q, d, u = F.required_tails(+1); assert 39.5 < q < 40.5 and 5.30 < d < 5.33 and u == 1.0
    assert 85 < F.required_tails(-1)[0] < 86.5
    a = F.required_angles(+1); assert abs(a["theta_L_d"] - 0.0505) < 2e-4 and abs(a["theta_R"] - 0.381) < 2e-3


def test_curve_flux_falls_short():
    assert 13.5 < F.best_flux_ratio() < 14.2 < F.required_tails(+1)[0]
    assert F.flux_ratio(13, 11, 0.214) < 1
