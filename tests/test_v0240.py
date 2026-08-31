"""v0.24.0 — Eisenstein/cuspidal release.

Covers: the CC-11..CC-15 correction family, the new exact lattice tools, the
level-generic Manin model, cuspidal subgroups with two-route certification
and external (Mazur) anchors, the AL morphology table, and the Wave-8
lambda_2 congruence conjecture.
"""
from fractions import Fraction as Fr

import numpy as np
import pytest

import mtft.al_morphology as AM
import mtft.cuspidal as CU
import mtft.hecke as HK
import mtft.integral_lattice as IL
import mtft.levels as LV


# ── CC-11 .. CC-14: the truncation family ───────────────────────────
# v0.23.0 coerced inputs with int(), silently zeroing Fractions.  Each of
# these fails against the old behaviour.

def test_cc11_saturate_rejects_rational_input():
    with pytest.raises(IL.InexactInputError):
        IL.saturate([[Fr(1, 2)], [Fr(1)]], [2])


def test_cc11_guard_is_live_not_dead():
    """The guard must fire before coercion, not after it."""
    with pytest.raises(IL.InexactInputError):
        IL._as_obj([[Fr(1, 3)]])


def test_cc12_solve_in_lattice_keeps_fractional_rhs():
    H = IL.hnf(np.array([[2, 0], [0, 2]], dtype=object))
    x = IL.solve_in_lattice(H, [Fr(1), Fr(1)])
    assert x == [Fr(1, 2), Fr(1, 2)]        # v0.23.0 returned [0, 0]


def test_cc13_rational_kernel_accepts_rational_matrix():
    k = IL.rational_kernel([[Fr(1, 2), Fr(-1)]])
    assert len(k) == 1
    a, b = k[0]
    assert Fr(1, 2) * a - b == 0            # v0.23.0 truncated 1/2 to 0


def test_cc14_class_order_scales_exactly():
    H = IL.hnf(np.array([[2, 0], [0, 2]], dtype=object))
    assert IL.class_order(H, [Fr(1), Fr(1)], [1, 2, 4]) == 2


def test_non_integral_float_is_rejected():
    with pytest.raises(IL.InexactInputError):
        IL._exact(1.5)


def test_clear_denominators_is_primitive():
    arr, mult = IL.clear_denominators([[Fr(1, 2), Fr(1, 3)]])
    assert arr.tolist() == [[3, 2]] and mult == [6]


# ── new exact lattice tools ─────────────────────────────────────────

@pytest.mark.parametrize("seed", [3, 11, 29])
def test_snf_transform_properties_and_agreement(seed):
    import random
    random.seed(seed)
    for _ in range(12):
        m, n = random.randint(1, 5), random.randint(1, 5)
        A = [[random.randint(-9, 9) for _ in range(n)] for _ in range(m)]
        U, S, V = IL.snf_transform(A)
        An = np.array(A, dtype=object)
        assert np.array_equal(U @ An @ V, S)
        d = [S[i][i] for i in range(min(m, n)) if S[i][i] != 0]
        assert all(d[i + 1] % d[i] == 0 for i in range(len(d) - 1))
        # cross-check against the independent HNF-ping-pong implementation
        assert [abs(x) for x in d] == IL.smith_invariants(A)
        K = IL.int_kernel(A)
        if K.shape[1]:
            assert not any((An @ K).flat)


# ── level-generic model ─────────────────────────────────────────────

@pytest.mark.parametrize("N,g,index,ncusps", [
    (11, 1, 12, 2), (15, 1, 24, 4), (23, 2, 24, 2),
    (33, 3, 48, 4), (35, 3, 48, 4), (143, 13, 168, 4),
])
def test_level_data(N, g, index, ncusps):
    d = LV.level_data(N)
    assert (d["genus"], d["index"], d["ncusps"]) == (g, index, ncusps)
    assert d["nu2"] == 0 and d["nu3"] == 0


def test_unsupported_levels_raise_with_reason():
    with pytest.raises(LV.UnsupportedLevelError, match="elliptic"):
        LV.check_supported(21)          # two order-3 points
    with pytest.raises(LV.UnsupportedLevelError, match="squarefree"):
        LV.check_supported(45)          # 3^2 * 5
    with pytest.raises(LV.UnsupportedLevelError, match="genus 0"):
        LV.check_supported(6)


def test_generic_model_reproduces_shipped_143_model():
    """The gate that licenses reading level-generic results onto 143."""
    mg, ms = LV.manin_model(143), HK.model()
    for k in ("P1", "sS", "sT", "sR", "iota", "E", "free", "nq",
              "freeK", "D2", "DEL", "cols", "D2r", "K", "Binv"):
        assert mg[k] == ms[k], f"field {k} differs"


@pytest.mark.parametrize("p", [2, 3, 5])
def test_generic_hecke_matches_shipped(p):
    assert LV.hecke_matrix(143, p) == HK.hecke_matrix(p)
    assert LV.cuspidal_hecke(143, p) == HK.cuspidal_hecke(p)


def test_cusp_labels_are_divisors_and_infinity_is_first():
    """CC-15: index 0 is the cusp infinity, label N — not label 1."""
    m = LV.manin_model(143)
    assert LV.cusp_labels(143) == [143, 13, 11, 1]
    assert sorted(m["cusp_label"]) == [1, 11, 13, 143]
    i_inf = m["idx"][m["canon"](0, 1)]
    assert m["cusp_label"][m["cusp_of"][i_inf]] == 143
    # fan length is the cusp width N/d
    for k, d in enumerate(m["cusp_label"]):
        assert len(m["fans"][k]) == 143 // d


# ── cuspidal subgroups ──────────────────────────────────────────────

def test_charpoly_matches_known_h6_factor():
    M = [[0, -1], [1, 0]]
    assert CU.charpoly(M) == [1, 0, 1]          # x^2 + 1


def test_cuspidal_group_143_reproduces_wave8():
    cg = CU.cuspidal_group(143)
    assert cg["structure"] == [10, 420]
    assert cg["order"] == 4200
    assert all(cg["gates"].values())
    assert cg["invariants_cokernel"] == [10, 420]


def test_cuspidal_group_two_routes_are_certified():
    cg = CU.cuspidal_group(143)
    assert cg["gates"]["two_route_agreement"] is True
    assert cg["order"] == cg["order_cokernel"]


def test_cc15_corrected_generator_names():
    """Generators are named by true divisor labels, base cusp = infinity."""
    cg = CU.cuspidal_group(143)
    assert cg["generators"][420] == "1*[C1 - C143]"
    assert cg["per_cusp_orders"] == {
        "[C13 - C143]": 70, "[C11 - C143]": 60, "[C1 - C143]": 420}


@pytest.mark.parametrize("p,expected", [(11, 5), (23, 11), (47, 23)])
def test_mazur_prime_level_anchor(p, expected):
    """External anchor: |C(J_0(p))| = numerator((p-1)/12)  (Mazur)."""
    assert Fr(p - 1, 12).numerator == expected
    assert CU.cuspidal_group(p)["order"] == expected


def test_two_torsion_image_equals_eisenstein_plane_at_143():
    tt = CU.two_torsion_image(143)
    assert tt["dim_E_Eis"] == 2
    assert tt["dim_C2_image"] == 2
    assert tt["equals_eisenstein"] is True


@pytest.mark.parametrize("N", [15, 33, 35])
def test_cross_level_control_c2_equals_eisenstein(N):
    """Level-universality control for the Wave-8 identity C[2] = E_Eis."""
    tt = CU.two_torsion_image(N)
    assert tt["contained_in_eisenstein"] is True
    assert tt["equals_eisenstein"] is True


def test_cross_level_control_table_shape():
    rows = CU.cross_level_control((15, 33))
    assert [r["N"] for r in rows] == [15, 33]
    assert all(r["routes_agree"] for r in rows)


# ── AL morphology ───────────────────────────────────────────────────

def test_al_traces_and_fixed_points():
    m = AM.morphology()
    assert m["traces"] == {11: 2, 13: -2, 143: -18}
    assert m["fixed_points"] == {11: 0, 13: 4, 143: 20}


def test_al_quotient_genus_cascade():
    m = AM.morphology()
    assert m["quotient_genera"] == {11: 7, 13: 6, 143: 2}
    assert m["genus_full_quotient"] == 1
    assert m["methods_agree"] is True       # RH and character formula agree


def test_w11_acts_freely():
    assert AM.morphology()["free_involutions"] == [11]


def test_cusps_form_a_free_torsor():
    t = AM.cusp_torsor()
    assert t["all_free"] is True and t["single_orbit"] is True


# ── Wave-8 open conjecture ──────────────────────────────────────────

def test_lambda2_congruence_small_range():
    r = CU.lambda2_congruence_scan(Lmax=37)
    assert r["disc_h6"] == 194616205        # odd => Z[a_2] is 2-maximal
    assert r["counterexamples"] == []
    assert r["primes_tested"] >= 8


@pytest.mark.slow
def test_lambda2_congruence_full_wave8_range():
    r = CU.lambda2_congruence_scan(Lmax=199)
    assert r["primes_tested"] == 43
    assert r["holds"] is True
