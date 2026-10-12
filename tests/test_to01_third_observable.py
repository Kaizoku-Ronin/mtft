"""TO-01: the pre-registered third observable."""
from mtft.matter import third_observable as T


def test_prediction_fails_window():
    lo, hi = T.vcb_range(); assert 0.078 < lo < 0.079 and 0.198 < hi < 0.199 and not T.passes()
    assert T.direct_entry_bound(18.6) < 3e-4 < T.sector_ratios()[0]


def test_required_asymmetry():
    thL, thR, ratio = T.lr_asymmetry(); assert abs(thL - 0.05) < 1e-12 and 7.5 < ratio < 7.9
