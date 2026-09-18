"""Composed gate battery (Wave C): field inventory -> anomaly polynomial -> charge lattice -> discrete anomaly -> flux/scalar gate
-> mode counts -> Yukawa chirality -> compactification.  Every gate returns a witness and its assumptions; there is no
"viable" boolean.  Gate 4, route 2: an MTFT-native parent (one parent, three internal modes) must pass every gate C3X
passes, without input copies; `m1_gate_report` shows which gates M1 fails today."""
import sympy as sp
from . import anomalies as AN, charge_lattices as CL, discrete_anomalies as DA, tensor_gs as TG, mode_operators as MO, compactification as CP, parents as PR

def m1_gate_report():
    N = dict(zip(PR.M1["stacks"], PR.M1["ranks"])); m = dict(zip(PR.M1["stacks"], PR.M1["degrees"])); model = AN.stack_model(N, m, order=PR.M1["stacks"])
    A = AN.anomaly_polynomials(model); K = AN.primitive_shift_matrix(model); q = A["q"]
    Y = sp.Matrix([sp.Rational(v) for v in PR.M1["hypercharge"]]); BL = sp.Matrix([sp.Rational(v) for v in PR.M1["B_minus_L"]]); phase = sp.ones(5, 1)
    gates = {}
    gates["4d_local_anomaly"] = {"pass": (K * sp.Matrix.hstack(phase, Y, BL)) == sp.zeros(2, 3) and AN.polarised_trace(A["cubic"], q, list(Y), list(Y), list(Y)) == 0, "witness": "K kernel = span{phase, Y, B-L}; cubic vanishes on it", "status": "EXACT"}
    gates["charge_lattice"] = {"pass": CL.smith_remnant(K)["invariant_factors"] == [1, 1], "witness": CL.smith_remnant(K), "status": "EXACT"}
    gates["majorana_dressing"] = {"pass": False, "witness": CL.integer_dressing(K, [0, 0, 2, 0, -2]), "meaning": "nu^c nu^c not dressable: bare Majorana mass forbidden", "status": "EXACT"}
    gates["native_scalar_flux"] = {"pass": TG.native_scalar_flux_gate(K, PR.M1["degrees"])["k_delta_m"] == 0, "witness": TG.native_scalar_flux_gate(K, PR.M1["degrees"]), "status": "EXACT (necessary only)"}
    tr = TG.flux_transgression(K); z3 = DA.spin_zn_fermion_test(DA.m1_ledger(), 3)
    gates["tensor_transgression_Z3"] = {"pass": z3["passes"], "witness": {"invariant_factors": tr["invariant_factors"], **z3}, "status": "EXACT (fermion-only)"}
    i13 = sp.I / sp.sqrt(13); fam = MO.s0_twist_cohomology([i13, -i13, i13], 0)
    gates["three_internal_modes"] = {"pass": fam["h0"] == 3 and fam["pure"], "witness": fam, "status": "EXACT"}
    gates["six_d_yukawa_chirality"] = {"pass": None, "witness": "M1 assigns no 6D chiralities; the section-product tensors have no specified 6D interaction (AXG-03 §7)", "status": "NOT_TESTED"}
    gates["six_d_anomaly_lift"] = {"pass": False, "witness": "AXG-02: irreducible p2 and mixed Abelian-color-cubic terms; AXG-03: tensor signature (2,1)", "status": "recorded from frozen studies"}
    V = CP.einstein_frame_potential(1, 1, 1, dim=6); gates["classical_radius"] = {"pass": V["stationary_point_exists"], "witness": str(V["V"]), "status": "EXACT (restricted potential)"}
    return {"model": "M1", "gates": gates, "route_2_todo": [k for k, g in gates.items() if g["pass"] is not True]}

def c3x_gate_report():
    P = TG.c3x_anomaly_polynomial(); lat = TG.integral_lattice_gate(TG.C3X_LATTICE["Omega"], TG.C3X_LATTICE["vectors"]); from .bordism import c3x_spin_bordism_certificate as B
    b = B(); M, R, rF = sp.symbols("M R rho_F", positive=True); bg = CP.c3x_ads_control(M, R, rF)
    gates = {"local_factorization": {"pass": TG.factorization_check(P["I8"], P["X4"], P["Y4"]), "witness": str(P["X4"]) + " * " + str(P["Y4"]), "status": "EXACT"},
             "integral_lattice": {"pass": lat["unimodular"] and lat["even"] and lat["all_norms_integral"], "witness": lat, "status": "EXACT"},
             "spin_bordism": {"pass": b["omega7_spin_BG"] == 0, "witness": b, "status": "EXACT (ordinary category)"},
             "global_GS": {"pass": None, "witness": "not constructed", "status": "OPEN"},
             "one_mode_per_copy": {"pass": MO.s0_twist_cohomology([sp.I / sp.sqrt(13), -sp.I / sp.sqrt(13)], 1)["pure"], "witness": "O(P+Q-R): (1,0)", "status": "EXACT"},
             "families": {"pass": None, "witness": "three INPUT copies", "status": "input, not derived"},
             "background": {"pass": None, "witness": {"lambda4": str(bg["lambda4"]), "ell4_squared": str(bg["ell4_squared"])}, "status": "AdS classical control; no scale separation"}}
    return {"model": "C3X", "gates": gates}

def route_2_requirements():
    """What an MTFT-native parent must supply and pass (Gate 4, route 2)."""
    return ["explicit global gauge group and 6D chiralities for every field", "one internal bundle with h^0 = 3, h^1 = 0 (purity) for each family sector — no input copies",
            "4D local anomaly cancellation on the retained U(1)s (exact)", "6D anomaly polynomial with an integral factorization I8 = X4 Y4 and an even unimodular tensor lattice",
            "Spin x Z_n and bordism gates for the discrete remnants", "opposite 6D chiralities for every Yukawa pair", "a specified action whose reduction yields the family Yukawas from overlaps",
            "a background solving the equations with the flux, with its stability tests declared", "at least one scale fixed or calibrated"]
