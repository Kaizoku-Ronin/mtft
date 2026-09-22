import sympy as sp
from fractions import Fraction
from mtft.research import gravitational_anomaly as GA, vector_higgs as VH, pipeline as PL

def test_kappa_is_one():
    r = VH.reduction_constants(); assert r["sigma_identity"] and r["kappa"] == 1

def test_tensor_integrality_and_survivors():
    assert GA.tensor_integrality(24)["net_self_dual_tensors_needed"] == 6 and GA.tensor_integrality(16)["net_self_dual_tensors_needed"] == 4
    assert not GA.tensor_integrality(22)["integral"] and not GA.tensor_integrality(18)["integral"]
    s = GA.m1_tensor_survivors(); assert sorted(s["survivors"]) == [(-1, -1, Fraction(4)), (1, 1, Fraction(6))] and not s["gravitino_possible"]
    # independent: the p2 coefficient of the four assignments from R2C-01 signed dimensions
    from mtft.research import chirality as CH
    for cd, ab, p in CH.m1_scalar_higgs_no_go()["ledger_preserving_p2_table"]:
        assert GA.tensor_integrality(p["n_grav"])["p2_fermions"] == p["p2_coefficient"]

def test_colour_cubic_obstruction():
    c = GA.m1_colour_cubic(); fL = GA.F["L"]; assert c["coefficients"]["L"] == -1 and c["f_L_irreducible"]
    assert sp.expand(c["dI8_dD3"] - (-fL - GA.F["a"] / 2 - GA.F["b"] / 2 + sp.Rational(5, 2) * GA.F["c"] - GA.F["d"] / 2)) == 0
    h = GA.hom_type_obstruction(); assert h["unchanged_by_extra_cd"] and len(h["escape_routes"]) == 2
    ff = GA.family_preserving_flux(0); assert ff["family_preserved"] and ff["m"] == {"c": 3, "L": 0, "a": 6, "b": 6, "d": 3} and ff["higgs_degrees"] == (-6, -6)   # CC-29
    w = GA.colour_cubic_without_U1L(); assert not w["escape_b_works"] and set(w["irreducible_terms"]) == {"a", "b"}
    g = PL.m1_gate_report()["gates"]["six_d_anomaly_lift"]; assert g["pass"] is False and g["witness"]["p2_tensor_integrality"]
