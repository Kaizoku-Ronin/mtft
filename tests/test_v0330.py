"""v0.33.0 gates: the theorem compendium's exact inputs applied to the toolset (2026-09-23).

Every assertion is exact (integers, Fractions, sympy) except the theta-function certificates, which are numerical on two routes.
"""
import itertools
from fractions import Fraction as Fr

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

from mtft import hecke as HK
from mtft.canonical import gates as CG
from mtft.research import gravitational_anomaly as GA, product_surface as PS, theta_torus as TT, unified_parent as UP, vector_higgs as VH


# ------------------------------------------------------------------ CC-33: the tensor coefficient
def test_cc33_tensor_units_and_rule():
    w = GA.tensor_coefficient_witnesses()
    assert w["weyl_units"] == {"weyl": 1, "self_dual_tensor": 28, "gravitino": 245}
    assert w["route_i_genera"] and w["route_ii_273"] and w["route_iii_2_0_multiplet"]
    assert [n for n in range(1, 60) if GA.tensor_integrality(n)["integral"]] == [28, 56]
    assert not any(UP.gaugino_p2_integrality(g)["integral"] for g in ("E6", "E7", "E8"))


# ------------------------------------------------------------------ Lemma E.1: point counts and gonality
def test_point_counts_two_routes_and_gonality():
    pc = HK.point_counts_from_T2()
    assert pc["dim_S2"] == 13 and pc["tr_T2"] == -1 and pc["tr_T2_squared"] == 39 and pc["F2"] == 4 and pc["F4"] == 18
    # independent of the module's arithmetic: the counts from the certified characteristic polynomial x^2 (x+2)^4 g^2 s^2
    x = sp.symbols("x"); g = sp.Poly(HK.G4[::-1], x); s6 = sp.Poly(HK.H6[::-1], x)
    roots = [0] + [-2] * 2 + [complex(r) for r in sp.Poly(g, x).nroots()] + [complex(r) for r in sp.Poly(s6, x).nroots()]
    assert len(roots) == 13 and abs(sum(roots) + 1) < 1e-9 and abs(sum(r * r for r in roots) - 39) < 1e-9
    ss = HK.supersingular_lower_bound(); assert ss["psi"] == 168 and ss["mass"] == 7 and ss["with_cusps_at_least"] == 18
    gb = HK.gonality_lower_bound(); assert gb["gonality_at_least"] == 4 and gb["excluded_degrees"] == [2, 3] and gb["two_routes_agree"]


def test_w13_quotient_petri_and_h0_2sumP():
    inv, quads = CG.w13_quotient_quadrics(); assert inv == [0, 7, 8, 9, 10, 11] and len(quads) == 6
    for q in quads:                                                                    # each quadric lies in I_2(X): residual of the products of the adapted forms
        E = __import__("mtft.canonical", fromlist=["adapted_qexpansions"]).adapted_qexpansions()
        for n in range(0, 60):
            assert sum(c * sum(E[a][CG.MONOMIALS[k][0]] * E[n - a][CG.MONOMIALS[k][1]] for a in range(n + 1)) for k, c in q.items()) == 0
    g = CG.gate_petri_w13_quotient()
    assert g["ok"] and (g["dim_I2_Y"], g["rank_S1_I2"], g["rank_S2_I2"]) == (6, 31, 91) and g["gonality_Y_at_least"] == 4 and g["h0_O_2sumP"] == 1


# ------------------------------------------------------------------ V.2–V.3: Kodaira form and the three geometries
def test_kodaira_forms_and_constant_curvature_checks():
    for d in (-6, -3, 3, 6):
        k = VH.kodaira_forms(d); assert k["tachyonic"]["lowest"] == -sp.Rational(abs(d), 24) and k["tachyonic"]["multiplicity"] == abs(d) + 12 and k["metric_independent"]
        c = VH.constant_curvature_checks(d); assert c["all_tachyonic_levels_match"]
        assert c["sphere"]["tachyonic"][1] == abs(d) - 1 and c["torus"]["tachyonic"][1] == abs(d) and not c["sphere"]["level_at_plus_B"] and c["genus13"]["level_at_plus_B"]
    b = VH.block_spectrum(-6); assert VH.kodaira_forms(-6)["tachyonic"]["multiplicity"] == 18


# ------------------------------------------------------------------ VI.3: own enumeration of the 64 ordered pairs
def test_surface_enumeration_from_scratch():
    types = [(a, b) for a in (3, -3, 1, -1) for b in (3, -3, 1, -1) if abs(a * b) == 3]
    sols, n_bideg, n_slope = [], 0, 0
    for f1, f2 in itertools.product(types, types):
        flags = (int(f1[0] < 0) + int(f2[0] < 0), int(f1[1] < 0) + int(f2[1] < 0))
        if flags not in ((1, 0), (0, 1)): n_bideg += 1; continue
        H = (-(f1[0] + f2[0]), -(f1[1] + f2[1]))
        if not H[0] * H[1] < 0: n_slope += 1; continue
        sols.append((f1, f2, H, sp.Rational(-H[0], H[1])))
    assert (n_bideg, n_slope, len(sols)) == (48, 12, 4)
    assert {(q, u, H) for q, u, H, _ in sols} == {(r["Q"], r["u"], r["H"]) for r in PS.surface_solutions()}
    assert all(abs(q[0]) != abs(u[0]) for q, u, _, _ in sols)                     # never curve x curve, never torus x torus
    m = PS.higgs_slope_mass((-4, 2), (1, 0)); assert m["equals_2pi_slope_over_areas"] and m["massless_locus"] == [2 * PS.A_E]
    m = PS.higgs_slope_mass((2, -4), (0, 1)); assert m["equals_2pi_slope_over_areas"] and m["massless_locus"] == [PS.A_E / 2]


# ------------------------------------------------------------------ VI.4: the rank bound, qualified
def test_up_mass_rank_examples():
    r = TT.up_mass_rank_examples(n_random=30)
    assert r["S2_ranks_random_VEVs"] == [2] and r["S1_rank_kunneth_rank_two_VEV"] == 3 and r["S1_ranks_product_VEVs"] == [2]
    rb = TT.up_mass_rank_bound(); assert "unconditional in S2" in rb["consequence"] and rb["rank_bound_general"] == "min(3, 2 * kunneth_rank(v))"


# ------------------------------------------------------------------ VI.5: the closed-form torus factor (no quadrature) against the pinned values
def test_theta_closed_form_pinned():
    tau = TT.tau_143a1_qseries(); assert abs(mp.re(tau) - 0.5) == 0 and abs(mp.im(tau) - mp.mpf("1.0232745926964612055995663")) < 1e-24
    assert abs(TT.theta_norm_closed_form(3, tau) ** 2 - mp.sqrt(mp.im(tau) / 6)) < 1e-25
    s1 = TT.torus_factor_closed_form(1, 2, tau); s2 = TT.torus_factor_closed_form(1, 3, tau)
    assert s1["rank"] == 2 and s2["rank"] == 3
    assert np.allclose(s1["singular_values"], [1.10302625437, 0.859234662985], atol=1e-10)
    assert np.allclose(s2["singular_values"], [1.12362454225, 0.978985971006, 0.843100702144], atol=1e-10)
    assert all(abs(s1["B"][1, b] - s1["B"][2, b]) < 1e-25 for b in range(2))            # rows j = 1, 2 equal (z -> -z)


@pytest.mark.slow
def test_theta_closed_form_matches_quadrature():
    tau = TT.tau_143a1(); q1 = TT.torus_factor(1, 2, tau, N=16); c1 = TT.torus_factor_closed_form(1, 2, tau)
    assert max(abs(q1["B"][i, j] - c1["B"][i, j]) for i in range(3) for j in range(2)) < 1e-12
