"""OM-01 and YB-01: the forced-odd-mode count and the bridged top."""
from mtft.matter import odd_modes as O
from mtft.matter import bridged_top as Y


def test_lefschetz_count_matches_pole_orders():
    for N in range(1, 13):
        ev, od = O.split_by_poles(N); assert ev + od == N and ev - od == O.lefschetz(N)
    assert O.first_forced_odd_degree() == 3 and O.split_by_poles(3) == (2, 1)


def test_no_odd_mode_on_the_b_stack():
    assert sorted(O.b_stack_blocks()) == ["H_d", "d^c", "e^c"] and O.odd_modes_on_b_stack() == []
    assert O.family_cost(3) == 72


def test_bridged_top_fails():
    assert 1.63 < Y.jet_factor() < 1.64 and 0.24 < Y.tail_value(1.79) < Y.tail_value(1.62) < 0.28
    o = Y.outcome(); assert not o["passes"] and 0.28 < max(o["flat"]) < 0.33 and o["shortfall"] > 1.55
    assert 0.16 < min(o["curved"]) < max(o["curved"]) < 0.28
