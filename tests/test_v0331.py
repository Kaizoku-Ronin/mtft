"""v0.33.1 gates (2026-09-27): the exact block census of split surface backgrounds (compendium V.1–V.3 applied leg by leg)."""
import sympy as sp
from mtft.research import product_surface as PS


def test_s2_census_at_the_wall_reproduces_the_coupled_relaxation_numbers():
    c = PS.block_census(sp.Rational(1, 2), "S2"); s = c["summary"]
    assert s["negative_exact_product_modes"] == 279 and s["harmonic_negative"] == 90          # the study's ">= 279, only 90 Dolbeault-closed"
    assert s["harmonic_massless"] == 2 * (4 + 56) + 14 * 14 and s["critical_r"] == [sp.Rational(1, 3), sp.Rational(1, 2), 3]
    by = {b["bidegree"]: b for b in c["blocks"]}
    assert by[(-3, 1)]["X_leg"]["ground"] == -5 and by[(-3, 1)]["X_leg"]["multiplicity"] == 15 and by[(-3, 1)]["X_leg"]["closed"]   # 6 x 15 = 90
    assert by[(2, -4)]["E_leg"]["ground"] == 0 and by[(2, -4)]["E_leg"]["multiplicity"] == 4 and by[(2, -4)]["E_leg"]["closed"]     # the 4 selected
    assert by[(-2, 4)]["X_leg"]["ground"] == 0 and by[(-2, 4)]["X_leg"]["multiplicity"] == 56 and by[(-2, 4)]["X_leg"]["closed"]    # the 56 opposite-root
    assert by[(-1, -3)]["E_leg"]["ground"] == -1 and not by[(-1, -3)]["E_leg"]["closed"]                                             # non-closed tachyons: 3 x 3 = 9
    assert not s["polystable"] and c["notes"] == []                                                                                    # S2 is fully exact


def test_slope_level_and_ground_levels_are_consistent():
    # a block's ground level equals 2 pi mu/(A_X A_E) on the leg whose form sits on the negative-degree factor (mixed signs); never polystable
    for model in ("S1", "S2"):
        for r in (sp.Rational(1, 3), sp.Rational(1, 2), 1, 2, 3):
            c = PS.block_census(r, model)
            assert not c["summary"]["polystable"]
            for b in c["blocks"]:
                a, bb = b["bidegree"]
                if a < 0 < bb: assert b["X_leg"]["ground"] == b["slope_level"]
                if bb < 0 < a: assert b["E_leg"]["ground"] == b["slope_level"]


def test_s1_wall_has_nine_harmonic_tachyons():
    s = PS.block_census(2, "S1")["summary"]
    assert s["harmonic_negative"] == 9 and s["negative_exact_product_modes"] == 99 and s["critical_r"] == [sp.Rational(1, 3), 2, 3]
    assert PS.block_census(2, "S1")["notes"] == []                                                # the degree-4 caveat is closed (h0_cm_divisor)
    assert PS.curve_h0_flux(6) == (1, None) and PS.curve_h0_flux(4) == (1, None) and PS.curve_h0_flux(-2) == (0, None)


def test_h0_of_cm_divisors_by_w13_eigenspaces():
    assert PS.h0_cm_divisor((1, 1, 1)) == {"h0_K_minus_D": 10, "h0_O_D": 1, "degree": 3}          # Lemma E.1
    assert PS.h0_cm_divisor((2, 2, 2)) == {"h0_K_minus_D": 7, "h0_O_D": 1, "degree": 6}           # Lemma E.2
    assert PS.h0_cm_divisor((2, 1, 1))["h0_O_D"] == 1 and PS.h0_cm_divisor((2, 2, 1))["h0_O_D"] == 1   # no g^1_4, g^1_5 through 2P1+P2+P3, 2P1+2P2+P3
    assert PS.h0_cm_divisor((1, 1, 1, 1)) is None and PS.h0_cm_divisor((3, 1, 1)) is None            # four points / order-3 conditions: not decided here


def test_s1_extension_graph_and_gluing_dimensions():
    g = {e["class_bundle"]: e for e in PS.extension_graph("S1", 2)["edges"]}
    assert g[(1, -3)]["ext1"] == 3 and g[(1, -3)]["level_E_leg"] == sp.Rational(-5, 2)            # Ext^1(L_c, L_1): the 9 tachyons (x3 colours)
    assert g[(-4, 2)]["ext1"] == 32 and g[(-4, 2)]["level_X_leg"] == 0 and g[(4, -2)]["ext1"] == 2  # Higgs classes and their conjugates, massless
    assert g[(-3, -1)]["ext1"] == 0 and g[(-3, -1)]["ext2"] == 15 and g[(3, 1)]["ext1"] == 10       # no colour–doublet recombination at all
    assert [PS.s1_colour_recombination_dims(k)["dim_Ext1_E_k_L_L"] for k in (1, 2, 3)] == [17, 2, 2]


def test_cp_structure_of_the_flux_choice():
    c = PS.cm_point_cp_structure(); assert c["CP_pairs"] == (("P1", "P2"), ("P3", "P4")) and c["W11_pairs"] == (("P1", "P3"), ("P2", "P4"))
    assert set(c["orbit_of_family_divisor"].values()) == {"omit P1", "omit P2", "omit P3", "omit P4"}
    t = PS.torus_cp_test(); assert t["S1"]["real_up_to_rephasing"] and t["S2"]["real_up_to_rephasing"]


def test_s1_iterated_extensions_are_unstable_below_r_five():
    for r in (sp.Rational(1, 2), 1, 2, 3, sp.Rational(49, 10)):
        assert PS.s1_iterated_extension_instability(r)["unstable"]
    assert not PS.s1_iterated_extension_instability(5)["unstable"] and not PS.s1_iterated_extension_instability(6)["unstable"]
    # the Kuenneth vanishings the proof uses, from the module's own cohomology helpers
    assert PS.curve_h0_flux(-3) == (0, None) and PS.torus_cohomology(0) == (1, 1) and PS.torus_cohomology(-3) == (0, 3)
    g = {e["class_bundle"]: e for e in PS.extension_graph("S1", 2)["edges"]}
    assert g[(-3, -1)]["ext1"] == 0 and g[(1, -3)]["ext1_E_leg"] == 3


def test_three_stack_trilemma():
    t = PS.three_stack_trilemma()
    assert t["Y+H+C"] == [] and len(t["Y+H"]) == 4 and all(w[0] * w[1] > 0 for *_, w in t["Y+H"])
    assert all(q == u for q, u, *_ in t["H+C"]) and len(t["Y+C"]) == 8 and all(H[0] * H[1] >= 0 for _, _, H, _ in t["Y+C"])
    assert PS.colour_slope_condition((3, 1), (1, -3))["w"] == (-5, -5) and PS.colour_slope_condition((-3, 1), (1, 3))["w"] == (7, 1)


def test_chirality_parity_on_the_spin_surface():
    t = PS.chirality_parity_theorem(box=4)
    assert t["odd_cases"] == [] and t["untwisted_nets"] == {0}
    assert all(f["vector_like"] and f["index_block"] == f["index_conjugate"] for f in t["families"].values())
    assert PS.adjoint_net_chirality((1, 1), (12, 1))["net_chirality"] == 2 and PS.adjoint_net_chirality((0, 1), (0, 0))["net_chirality"] == -24


def test_rank2_net_index_lemma_and_empty_recombination():
    r = 3
    assert PS.rank2_net_index_lemma([(0, 0)], [(3, -1)], r)["net_index"] == -3                               # a line pair on the ray
    lem = PS.rank2_net_index_lemma([(0, 0)], [(5, -2), (1, 0)], r)                                            # rank 2: m = (3,-1) on the ray, extension block (4,-2)
    assert lem["holds"] and lem["net_index"] == -10
    import itertools, fractions
    box = [(a, b) for a in range(-6, 7) for b in range(-4, 5)]
    found = 0
    for dx, dxp in itertools.combinations(box, 2):                                                           # option A at r = 3, Q = (3, -1)
        if dx[0] + 3 * dx[1] + dxp[0] + 3 * dxp[1] != 0: continue
        e2 = (dx[0] - dxp[0], dx[1] - dxp[1])
        if e2[0] * e2[1] > 0 or dx[0] + 3 * dx[1] == 0: continue                                              # no extension in either direction
        if dx[0] * dx[1] + dxp[0] * dxp[1] == 3: found += 1                                                   # u^c net +3 opposite to Q's -3
    assert found == 0
