import sympy as sp
from fractions import Fraction
from mtft.research import gravitational_anomaly as GA, vector_higgs as VH, pipeline as PL

def test_kappa_is_one():
    r = VH.reduction_constants(); assert r["sigma_identity"] and r["kappa"] == 1

def test_tensor_coefficient_cc33():
    # the constant is re-derived here, independently of the module: -L8/8 with L8 = (7 p2 - p1^2)/45, and A-roof_8 = (7 p1^2 - 4 p2)/5760
    assert GA.P2_WEYL == Fraction(-4, 5760) == Fraction(-1, 1440) and GA.P2_SELF_DUAL_TENSOR == Fraction(-7, 360) and GA.P1SQ_SELF_DUAL_TENSOR == Fraction(1, 360)
    assert GA.P2_SELF_DUAL_TENSOR == 28 * GA.P2_WEYL and GA.P2_GRAVITINO == 245 * GA.P2_WEYL            # Weyl units 1 : 28 : 245
    assert 245 + 28 == 273 and 28 + 1 == 29 and (28 + 2) * (-GA.P2_WEYL) == Fraction(1, 48)              # H - V + 29 T = 273; (2,0) multiplet
    w = GA.tensor_coefficient_witnesses(); assert w["route_i_genera"] and w["route_ii_273"] and w["route_iii_2_0_multiplet"]

def test_tensor_integrality_and_survivors():
    # CC-33: cancellation needs n_grav = 0 mod 28; none of M1's ledger-preserving assignments (24, 22, 18, 16) survives
    for n in (24, 22, 18, 16): assert not GA.tensor_integrality(n)["integral"] and GA.tensor_integrality(n)["net_self_dual_tensors_needed"] == Fraction(-n, 28)
    assert GA.tensor_integrality(28)["integral"] and GA.tensor_integrality(28)["anti_self_dual_tensors_needed"] == 1 and GA.tensor_integrality(56)["integral"]
    s = GA.m1_tensor_survivors(); assert s["survivors"] == [] and len(s["table"]) == 4 and not s["gravitino_possible"]
    assert GA.tensor_integrality(24, gravitino=True)["net_self_dual_tensors_needed"] == Fraction(-269, 28)
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
    g = PL.m1_gate_report()["gates"]["six_d_anomaly_lift"]; assert g["pass"] is False and g["witness"]["p2_tensor_integrality"]["survivors"] == [] and len(g["witness"]["p2_tensor_integrality"]["table"]) == 4
