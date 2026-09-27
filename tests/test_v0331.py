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
    assert any("g^1_4" in n for n in PS.block_census(2, "S1")["notes"])                       # the degree-4 caveat is declared, not hidden
    assert PS.curve_h0_flux(6) == (1, None) and PS.curve_h0_flux(4)[1] is not None and PS.curve_h0_flux(-2) == (0, None)
