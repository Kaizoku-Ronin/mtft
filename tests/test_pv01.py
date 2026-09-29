"""PV-01 (2026-09-29): the 8D ten-form anomaly gate, the B1 block/index census and scan, the flux-energy landscape, and mtft.parent.
Fast tier exact; the box-3 three-singlet scan is marked slow (~45 s)."""
import pytest
import sympy as sp
from fractions import Fraction
from mtft.research import anomaly8d as A8, vacuum_energy as VE, product_surface as PS
from mtft import parent as P


# ------------------------------------------------------------------ anomaly8d
def test_lambda2_and_sym2_quintic_rules():
    N = A8.rank("A"); s = lambda k: A8.s("A", k)
    assert sp.expand(120 * A8.antisym2("A").c[5] - ((N - 16) * s(5) + 5 * s(1) * s(4) + 10 * s(2) * s(3))) == 0
    assert sp.expand(120 * A8.sym2("A").c[5] - ((N + 16) * s(5) + 5 * s(1) * s(4) + 10 * s(2) * s(3))) == 0
    assert sp.expand(6 * A8.antisym2("A").c[3] - ((N - 4) * s(3) + 3 * s(1) * s(2))) == 0        # the 4D (N - 4) rule
    assert all(A8.adjoint("A").c[k] == 0 for k in (1, 3, 5))                                    # real representation: no anomaly


def test_sun_conditions_and_no_exact_cancellation():
    c = A8.sun_content_conditions()
    nF, nA, nS = sp.symbols("n_F n_A n_S"); N = A8.rank("A")
    assert sp.expand(120 * c["quintic_coefficient"] - (nF + (N - 16) * nA + (N + 16) * nS)) == 0
    assert c["reducible_is_s3_times_X4"]
    assert c["exact_cancellation_solutions"] == [{nA: 0, nF: 0, nS: 0}]


def test_u6_lambda2_plus_ten_fundamentals_is_quintic_free_with_two_gs_fields():
    I = A8.I10([(A8.antisym2("A", 6), 1), (A8.fund("A", 6), 10)])
    d = A8.gs_decomposition(I)
    assert d["cancellable"] and d["green_schwarz_fields_needed"] == 2
    assert d["two_form_factors"][0][0] == A8.s("A", 3)                                           # X_6 = tr F^3
    roots = A8.specialise(I.subs(A8.s("A", 1), 0), {"A": 6})                                     # SU(6): explicit Chern roots
    f = sp.factor(roots); assert f != 0 and len(sp.factor_list(roots)[1]) >= 2                   # factorised (X_4 X_6), nonzero


def test_sun_gs_census_rule():
    rows = A8.sun_gs_census(6, box=1)
    assert next(r["n_F"] for r in rows if (r["n_A"], r["n_S"]) == (1, 0)) == 10
    assert all(r["n_F"] == -(6 - 16) * r["n_A"] - (6 + 16) * r["n_S"] for r in rows)


def test_split_stack_parity_lemma_and_index_table():
    lem = A8.split_stack_parity_lemma(); assert lem["squares_even"]
    for m in ("S1", "S2"):
        t = A8.b1_index_table(m)
        assert all(not r["odd"] for rep in ("antisym2", "sym2") for r in t["blocks"][rep] if r["origin"] == "L_i^2")
        assert all(r["index"] == 0 for r in t["blocks"]["fund"] if r["content"].startswith("c:"))     # trivial colour stack: no coloured families from L_c
        assert {abs(x[3]) for x in t["index_pm3_blocks"]} == {3}


def test_b1_content_blocks_conjugation_and_scan_negatives():
    st = (("c", 3, (0, 1)), ("L", 2, (3, 0)), ("s1", 1, (-1, 2)), ("s2", 1, (-1, 2)), ("s3", 1, (0, 1)))
    b = A8.b1_content_blocks(st, 8)
    assert b[("3", "2", (1, 1, 0, 0, 0))] == 3 and b[("3bar", "1", (-1, 0, -1, 0, 0))] == 3
    assert A8.b1_scan(k_singlets=1, box=3)["n_accepted"] == 0
    assert A8.b1_scan(k_singlets=2, box=3)["n_accepted"] == 0
    assert A8.b1_scan(k_singlets=2, box=3, n_anti=0, n_sym=1)["n_accepted"] == 0


@pytest.mark.slow
def test_b1_scan_three_singlet_stacks_box3():
    r = A8.b1_scan(k_singlets=3, box=3)
    assert r["n_accepted"] == 60 and r["n_fund"] == 8
    for acc in r["accepted"]:
        core = {k: v for k, v in acc["net"].items() if not (k[0] == "1" and k[1] == "1" and k[2] == 0)}
        assert core == A8.SM_NET(3)
    tri = A8.b1_trilemma(r)
    assert (tri["Y+H"], tri["Y+C"], tri["Y+H+C_common_r"], tri["H_and_C_at_different_r"]) == (6, 6, 0, [])


@pytest.mark.slow
def test_b1_scan_three_singlet_stacks_box4():
    r = A8.b1_scan(k_singlets=3, box=4); tri = A8.b1_trilemma(r)
    assert r["n_accepted"] == 156 and tri["Y+H+C_common_r"] == 0
    assert set(tri["H_and_C_at_different_r"]) == {(Fraction(2), Fraction(5, 2)), (Fraction(1, 2), Fraction(2, 5))}


def test_hypercharge_solver():
    blocks = {("3", "2", (1, 1, 0)): 3, ("3bar", "1", (-1, 0, -1)): 3, ("1", "2", (0, 1, 1)): 3}
    fit = A8.hypercharge_assignment(blocks, {("3", "2", (1, 1, 0)): Fraction(1, 6), ("3bar", "1", (-1, 0, -1)): Fraction(-2, 3), ("1", "2", (0, 1, 1)): Fraction(-1, 2)})
    assert fit["pinned"] and fit["y"] == [Fraction(2, 3), Fraction(-1, 2), Fraction(0)]
    assert A8.hypercharge_assignment(blocks, {("3", "2", (1, 1, 0)): Fraction(1, 6), ("3", "2", (1, 1, 0)): Fraction(1, 3)}) is not None   # dict dedups keys


# ------------------------------------------------------------------ vacuum energy
def test_line_bundle_energy_amgm_floor():
    e = VE.line_bundle_energy(3, -1); assert e["identity"] and e["floor_over_4pi2"] == 6 and e["self_dual_ratio"] == 3
    assert sp.simplify(e["E_over_4pi2"].subs(VE.R, 3)) == 6


def test_split_energy_minimiser_and_mirrors():
    s1, s2 = VE.flux_energy(model="S1"), VE.flux_energy(model="S2")
    assert (s1["A"], s1["B"]) == (19, 11) == (s2["A"], s2["B"]) and s1["r_star"] == sp.sqrt(209) / 11
    assert s1["E_star_over_4pi2"] == 2 * sp.sqrt(209) and s1["floor_over_4pi2"] == 18 and s1["polystable_r"] == []
    assert sp.simplify(sp.diff(s1["E_over_4pi2"], VE.R).subs(VE.R, s1["r_star"])) == 0


def test_hym_floor_and_polystable_control():
    h1, h2 = VE.hym_floor(model="S1"), VE.hym_floor(model="S2")
    assert (h1["Delta"], h1["c1"], h1["r_hym"], h1["E_HYM_min_over_4pi2"]) == (14, (-5, -5), 1, sp.Rational(32, 3))
    assert (h2["Delta"], h2["c1"], h2["r_hym"], h2["E_HYM_min_over_4pi2"]) == (50, (7, 1), 7, sp.Rational(32, 3))
    assert sp.Poly(h1["excess_polynomial"], VE.R).discriminant() < 0 and sp.Poly(h2["excess_polynomial"], VE.R).discriminant() < 0   # never touches the floor
    st = (("a", 1, (3, -1)), ("h", 1, (-6, 2)), ("b", 2, (0, 0)))                                                   # the M4 triangle: polystable at r = 3
    assert VE.flux_energy(stacks=st)["polystable_r"] == [3] and VE.hym_floor(3, stacks=st)["excess_split_minus_HYM_at_r"] == 0


def test_einstein_frame_no_volume_minimum():
    ef = VE.einstein_frame_scaling(); vol = sp.Symbol("Vol", positive=True)
    assert ef["volume_stationary_point"] is False and ef["flux_exponent"] == -2 and ef["curvature_exponent"] == sp.Rational(-3, 2)
    assert all(t.is_negative is not False for t in sp.Add.make_args(ef["dV_dVol"]))                      # every term negative: monotone in Vol


# ------------------------------------------------------------------ mtft.parent
def test_parent_scoreboard_matches_ledger_status_column():
    g = P.gates(P.THREE_STACK_8D_S1)
    assert [g[k]["pass"] for k in ("1_odd_chirality", "2_yukawa", "4_colour_slope", "5_no_twisted_colour_destabiliser", "6_tachyon_free")] == [False] * 5
    assert g["7_light_higgs"]["pass"] and g["9_anomaly"]["pass"] and g["10_hypercharge_normalisation"]["pass"]
    assert P.gates(P.M3_E7_6D)["9_anomaly"]["pass"] is False                                             # CC-33: 133 not 0 mod 28


def test_b1_exemplars_exhibit_the_trilemma():
    yh, yc = P.gates(P.B1_U8_YH), P.gates(P.B1_U8_YC)
    for g in (yh, yc):
        assert g["1_odd_chirality"]["pass"] and g["2_yukawa"]["pass"] and g["3_uc_dc_summands"]["pass"] and g["9_anomaly"]["pass"]
        assert g["9_anomaly"]["witness"]["green_schwarz_fields_needed"] == 2 and g["6_tachyon_free"]["pass"] is False
    assert yh["7_light_higgs"]["pass"] and yh["4_colour_slope"]["pass"] is False and yh["7_light_higgs"]["witness"]["slope_free_loci"] == {"u": 2, "d": 2}
    assert yc["4_colour_slope"]["pass"] and yc["4_colour_slope"]["witness"]["r_colour"] == sp.Rational(3, 2) and yc["7_light_higgs"]["pass"] is False
    hc = P.gates(P.B1_U8_HC)
    assert hc["7_light_higgs"]["witness"]["slope_free_loci"] == {"u": 2, "d": 2} and hc["4_colour_slope"]["witness"]["r_colour"] == sp.Rational(5, 2)
    assert "T" in P.ledger_table()["table"]
