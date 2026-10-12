"""DS-01: the down sector is parity-symmetric."""
from mtft.matter import down_parity as D


def test_touch_pattern_and_exact_relation():
    t = D.touched(); assert t["Q"] == t["d^c"] == t["H_d"] == 0 and t["u^c"] == t["e^c"] == 1 and t["nu^c"] == 2 and t["H_u"] == t["L"] == -1
    for b, bp in ((0.04, 0.3), (0.0164, 1.17), (0.02, 1.0)):
        r, t2 = D.angles_2x2(b, bp); assert abs(r / t2 - 1) < 1e-9


def test_floor_and_solution():
    assert 0.46 < D.c_floor() < 0.48 and 0.0042 < D.vub_floor() < 0.0044 and D.c_floor() > 0.45
    s = D.solve(); assert 1.16 < s["tan_thR"] < 1.18 and 0.049 < s["sin_thL_u"] < 0.052 and 0.06 < s["z"] < 0.07 and 47 < s["phase_deg"] < 50
    assert 70 < D.tail_ratio_needed() < 72 and 0.90 < D.c_with_flux_cap() < 0.92
