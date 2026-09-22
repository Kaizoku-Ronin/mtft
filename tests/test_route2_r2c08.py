import sympy as sp
from mtft.research import product_surface as PS

def test_selection_rule_and_scan():
    assert not PS.yukawa_selection_rule((3, -1), (3, -1))["yukawa_allowed"] and PS.yukawa_selection_rule((3, -1), (-1, -3))["yukawa_allowed"]
    assert PS.kunneth_parity(3, -1) == 1 and PS.kunneth_parity(-1, -3) == 0 and PS.kunneth_parity(3, 1) == 0
    cc = PS.curve_curve_no_go(); assert cc["pairs"] == 8 and cc["no_go"]
    sols = PS.mixed_origin_solutions(); assert len(sols) == 8 and all(set(r["origins"]) == {"curve", "torus"} for r in sols)
    loci = {r["slope_free_locus"] / PS.A_E for r in sols}; assert loci == {2, sp.Rational(1, 2)}
    higgs = {r["H"] for r in sols}; assert higgs == {(-4, 2), (-2, 4), (2, -4), (4, -2)}
    # independent bundle-degree count of the selection rule: (0,q1)+(0,q2)+(0,1) must be a (0,2)-form -> q1 + q2 = 1
    assert all((q1 + q2 == 1) == (q1 != q2) for q1 in (0, 1) for q2 in (0, 1))
    r = PS.m4_record(); assert "NONE" in r["yukawa"] and PS.rank_structure()["rank_bound"] == "rank(v)"
