"""v0.24.0 — Arf parity direction geometry, finite Kakeya, and the corrected
cross-level control (CC-16, CC-17).

The cross-level tests here pin BOTH directions: the identity C[2] = E_Eis
holds at the two-prime-factor levels and FAILS at N = 105.  Both are asserted,
so a future change that quietly makes 105 "pass" will break the suite.
"""

import pytest

import mtft
from mtft import kakeya as KY


# ---------------------------------------------------------------- Arf theorem

def test_arf_direction_theorem_x0143():
    r = KY.arf_direction_theorem()
    assert (r["n_characteristics"], r["n_even"], r["n_odd"]) == (128, 96, 32)
    assert (r["dim_V"], r["dim_R"]) == (7, 5)
    assert r["n_directions_odd"] == 31 == r["expected_odd"]
    assert r["n_directions_even"] == 127 == r["expected_even"]
    assert r["delta_odd_equals_radical"]
    assert r["delta_even_equals_V"]


def test_arf_forcing_argument_is_derived_not_assumed():
    """q|_R = 0 is forced by the counts: n_odd == 2^dim R, not 2^dim V / 2."""
    r = KY.arf_direction_theorem()
    assert r["q_restricted_to_radical_is_zero"]
    assert r["odd_locus_is_single_R_coset"]
    assert r["n_odd"] == 2 ** r["dim_R"]
    assert r["n_odd"] != 2 ** r["dim_V"] // 2      # the q|_R != 0 alternative


def test_odd_locus_is_directionally_confined():
    r = KY.arf_direction_theorem()
    assert r["odd_confined"]
    assert r["n_directions_odd"] < r["n_directions_even"]


def test_not_advertised_as_a_kakeya_theorem():
    """Guard against promotion: this is a difference-direction statement."""
    assert KY.arf_direction_theorem()["is_kakeya_theorem"] is False


def test_direction_set_of_an_affine_coset_is_its_direction_space():
    rows = [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1]]
    D = KY.direction_set(rows)
    assert D == {(0, 0, 1), (0, 1, 0), (0, 1, 1)}


# ---------------------------------------------------------------- finite Kakeya

@pytest.mark.parametrize("p", [5, 7, 11, 13])
def test_besicovitch_is_a_kakeya_set_and_respects_dvir(p):
    rep = KY.kakeya_report(p)
    assert rep["is_kakeya"]
    assert rep["size"] >= KY.dvir_bound(p)
    assert rep["excess_over_bound"] == (p - 1) // 2


def test_dvir_bound_values():
    assert KY.dvir_bound(11) == 66
    assert KY.dvir_bound(13) == 91


def test_kakeya_set_is_a_strict_subset_of_the_plane():
    rep = KY.kakeya_report(11)
    assert rep["size"] < rep["ambient"]


# ---------------------------------------------------------------- CRT bridge

def test_crt_direction_bridge_11_13():
    b = KY.crt_direction_bridge(11, 13)
    assert b["size_P1_N"] == 168
    assert (b["size_P1_p"], b["size_P1_q"]) == (12, 14)
    assert b["product"] == 168
    assert b["bijective"]
    assert b["equals_index_Gamma0"]


def test_p1_counts():
    assert len(KY.p1_points(11)) == 12
    assert len(KY.p1_points(13)) == 14
    assert len(KY.p1_points(143)) == 168


# ------------------------------------------------- cross-level control (CC-16/17)

@pytest.mark.slow
def test_cross_level_rows_carry_a_vacuity_flag():
    from mtft import cuspidal as CU
    rows = {r["N"]: r for r in CU.cross_level_control(levels=(15, 33))}
    assert rows[15]["informative"] is False          # genus 1: E_Eis is everything
    assert rows[15]["vacuity_reason"]
    assert rows[33]["informative"] is True


@pytest.mark.slow
def test_c2_equals_eisenstein_at_two_prime_factor_levels():
    from mtft import cuspidal as CU
    for r in CU.cross_level_control(levels=(33, 35, 143)):
        assert r["n_prime_factors"] == 2
        assert r["ncusps"] == 4
        assert r["dim_E_Eis"] == 2
        assert r["C2_equals_E_Eis"], f"expected identity to hold at N={r['N']}"


@pytest.mark.slow
def test_c2_identity_FAILS_at_three_prime_factors():
    """CC-17: the identity is NOT level-universal.

    N = 105 = 3*5*7 has eight cusps and the same genus as 143.  There
    dim E_Eis = 8, dim C[2] = 6, and the identity fails.  This negative is
    asserted so it cannot be silently lost.
    """
    from mtft import cuspidal as CU
    r = CU.cross_level_control(levels=(105,))[0]
    assert r["n_prime_factors"] == 3
    assert r["ncusps"] == 8
    assert r["genus"] == 13
    assert r["informative"] is True
    assert r["dim_E_Eis"] == 8
    assert r["dim_C2"] == 6
    assert r["C2_equals_E_Eis"] is False


def test_version_triple_is_current():
    assert mtft.__version__ == "0.24.0"
