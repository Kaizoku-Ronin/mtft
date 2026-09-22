import sympy as sp
from mtft.research import yukawa_triangle as YT, pipeline as PL

def test_telescoping_theorem_and_m1_triangles():
    t = YT.triangle_closure(3, 3); assert t["d_H"] == -6 and t["higgs_tachyon_m2"] == -sp.Rational(1, 4) and not t["flux_neutral_possible"]
    assert YT.triangle_closure(3, -3)["flux_neutral_possible"]
    m = YT.m1_triangles(); assert m["closed"] and all(sum(v) == 0 for v in m["degrees"].values()) and m["higgs_degree"] == -6
    assert m["higgs_hypercharges"] == {"H_u": sp.Rational(1, 2), "H_d": -sp.Rational(1, 2), "H_d(lep)": -sp.Rational(1, 2), "H_u(nu)": sp.Rational(1, 2)}
    # independent: the hypercharges of the M1 stacks (1/6, 0, -1/2, 1/2, -1/2) give Y(L,a-bar) = 0 - (-1/2) = 1/2 and Y(L,b-bar) = -1/2
    y = {"c": sp.Rational(1, 6), "L": 0, "a": -sp.Rational(1, 2), "b": sp.Rational(1, 2), "d": -sp.Rational(1, 2)}; assert y["L"] - y["a"] == sp.Rational(1, 2) and y["L"] - y["b"] == -sp.Rational(1, 2)

def test_e7_fails_yukawa_gate_and_pipeline_records_no_go():
    assert YT.e7_triangle_test()["yukawa_gate"] == "FAIL" and YT.e7_triangle_test((0, 1, -1, 2, -2))["yukawa_gate"] == "PASS"
    g = PL.m1_gate_report()["gates"]["light_higgs"]; assert g["pass"] is False and g["witness"]["higgs_degree"] == -6
    assert len(YT.light_higgs_no_go()["escapes"]) == 2
