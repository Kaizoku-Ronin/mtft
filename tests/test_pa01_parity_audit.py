"""PA-01: the parity audit."""
from mtft.matter import parity_audit as P


def test_no_down_sector_odd_source():
    assert P.down_sector_odd_sources() == []
    assert sum(1 for v in P.SOURCES.values() if v["odd"]) == 4 and P.odd_moduli() == {"even": 37, "odd": 13}


def test_loophole_and_verdict():
    lo = P.loophole(); assert 0.17 < lo["v_needed"] < 0.19
    assert all(20 < lo["cases"][s]["ratio_to_data"] < 70 for s in (-1, 1))
    v = P.verdict(); assert 0.90 < v["c_floor"] < 0.92 and 2.2 < v["factor_over_data"] < 2.4
    assert 0.0042 < P.floor_vub() < 0.0044
