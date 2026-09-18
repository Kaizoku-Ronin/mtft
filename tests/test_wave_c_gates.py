"""Wave C: model records and the composed gate battery; every gate returns a witness.  Expected values come from independent
routes (charge sums, Smith forms, F2 linear algebra on the recorded Sq^2 matrices, the direct Einstein solve)."""
import sympy as sp
from mtft.research import parents as PR, discrete_anomalies as DA, tensor_gs as TG, compactification as CP, bordism as BO, pipeline as PL


def test_records_are_immutable_and_distinct():
    assert set(PR.RECORDS) == {"M1", "ONE_PARENT_INDEX3", "U8_ADJOINT", "C3X"} and PR.M1["multiplicity_origin"] != PR.C3X["multiplicity_origin"]
    try: PR.M1["ranks"] = (1,); raised = False
    except TypeError: raised = True
    assert raised


def test_z3_fermion_ledger_matches_axg03_and_fails_hsieh():
    z = DA.spin_zn_fermion_test(DA.m1_ledger(), 3); assert (z["S1"], z["S3"]) == (24, 24) and z["cubic_residue_mod_6n"] == 12 and z["S3"] % 9 == 6 and not z["passes"]
    assert DA.m1_z3_generator_certificate()["certifies_Z3"]
    # cross-check S3 against the anomaly polynomial evaluated on the generator charges (independent route)
    from mtft.surface import smflux as SF
    from fractions import Fraction
    S = SF.M1_STACKS; q = dict(zip(S["N"], (0, 0, 1, 0, 1))); assert Fraction(SF.abelian_anomaly_polynomial(S["N"], S["m"], q)) == 24


def test_flux_gates_and_transgression():
    K = sp.Matrix([[0, -2, 1, 1, 0], [3, -4, 0, 0, 1]]); m = (0, -3, 3, 3, 0)
    g = TG.native_scalar_flux_gate(K, m); assert g["K_m"] == [12, 12] and g["k_delta"] == [3, -2, -1, -1, 1] and g["k_delta_m"] == 0
    t = TG.flux_transgression(K); assert t["K_hat"] == sp.Matrix([[3, -2, -1, -1, 1], [0, -6, 3, 3, 0]]) and t["invariant_factors"] == [1, 3] and t["discrete_remnant"] == [3]


def test_c3x_factorization_lattice_bordism():
    P = TG.c3x_anomaly_polynomial(); assert TG.factorization_check(P["I8"], P["X4"], P["Y4"]); C, W, eta, h, x = P["symbols"]
    assert sp.expand(P["I6_flux_reduction"] + 12 * x * P["X4"]) == 0 and P["BF_level"] == -12
    lat = TG.integral_lattice_gate(TG.C3X_LATTICE["Omega"], TG.C3X_LATTICE["vectors"]); assert lat["unimodular"] and lat["even"] and lat["signature"] == (1, 1) and lat["all_norms_integral"]
    b = BO.c3x_spin_bordism_certificate(); assert b["dims"] == (5, 9, 16) and (b["rank_H4_H6"], b["rank_H6_H8"]) == (2, 7) and b["composite_zero"] and b["kernel_equals_image"] and b["E3_6_1"] == 0 and b["omega7_spin_BG"] == 0


def test_background_by_direct_einstein_solve():
    M, U, rF, R = sp.symbols("M U rho_F R", positive=True); lam, k = sp.symbols("lambda4 k")
    T4 = -(U + rF); T2 = rF - U; sol = sp.solve([M * (-lam - k) - T4, M * (-2 * lam) - T2], [lam, k])          # AXG-04's direct Einstein solve
    bg = CP.product_background(M, U, rF); assert sp.simplify(bg["lambda4"] - sol[lam]) == 0 and sp.simplify(bg["k"] - sol[k]) == 0
    c = CP.c3x_ads_control(M, R, rF); assert sp.simplify(c["k"] + 1 / R ** 2) == 0 and sp.simplify(c["lambda4"] + 1 / R ** 2 + 2 * rF / M) == 0
    assert sp.simplify(3 * R ** 2 - c["ell4_squared"]).is_positive and c["radion_mass2"].subs({M: 1, R: 1, rF: 1}) == 8
    assert CP.einstein_frame_potential(1, 1, 1, dim=6)["stationary_point_exists"] is False and CP.einstein_frame_potential(1, 1, 1, dim=7)["stationary_point_exists"] is False


def test_gate_reports_return_witnesses_not_verdicts():
    r = PL.m1_gate_report(); assert all("witness" in g and "status" in g for g in r["gates"].values()) and "viable" not in r
    assert r["gates"]["three_internal_modes"]["pass"] and r["gates"]["4d_local_anomaly"]["pass"] and not r["gates"]["tensor_transgression_Z3"]["pass"]
    assert set(r["route_2_todo"]) >= {"tensor_transgression_Z3", "six_d_anomaly_lift", "classical_radius"}
    c = PL.c3x_gate_report(); assert c["gates"]["local_factorization"]["pass"] and c["gates"]["spin_bordism"]["pass"] and c["gates"]["global_GS"]["pass"] is None and c["gates"]["families"]["pass"] is None
    assert len(PL.route_2_requirements()) >= 8
