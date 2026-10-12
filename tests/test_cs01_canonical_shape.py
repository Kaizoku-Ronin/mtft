"""CS-01: the pre-registered canonical-shape test."""
from mtft.matter import canonical_shape as C


def test_candidates_and_rule():
    c = C.candidates()
    assert c["pullback_phi"] == (4, 4) and abs(c["hyperbolic_over_neron"][1] - 44.30) < 0.01 and abs(c["hyperbolic_over_neron"][0] - 22.15) < 0.01
    assert not any(C.hit(v) for v in c.values())


def test_invariants_and_chance():
    assert C.report()["invariants_in_window"] == ["|chi|", "2g", "hyp_area_over_2pi"]
    assert 0.09 < C.chance_per_candidate() < 0.11 and 0.3 < C.report()["expected_accidental"] < 0.4
