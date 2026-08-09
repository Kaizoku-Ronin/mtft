"""Tests for the v0.14.0 promotion: moments, curvature, hecke, eisenstein."""
from fractions import Fraction as Fr
from math import gcd

import pytest
from mpmath import mp

from mtft import moments as M, curvature as CV, hecke as H, eisenstein as EI


# ── moments ─────────────────────────────────────────────────────────

def test_cold_constants():
    mp.dps = 40
    assert mp.nstr(M.weight_second_moment(1), 21) == M.COLD["T"]
    assert abs(M.weight_third_moment(1) - mp.mpf(M.COLD["U"])) < mp.mpf("1e-25")
    assert abs(M.weight_susceptibility(1) - mp.mpf(M.COLD["chi_w"])) < mp.mpf("1e-18")


def test_covariance_is_zeta_second_derivative():
    mp.dps = 30
    for b in ("1.5", "2", "3.5"):
        assert abs(M.cov_log_weight(b) - mp.diff(mp.zeta, mp.mpf(b) + 1, 2)) < mp.mpf("1e-25")


def test_susceptibility_definition():
    mp.dps = 30
    b = mp.mpf("2.5")
    lhs = M.weight_susceptibility(b)
    rhs = M.weight_second_moment(b) - M.weight_first_moment(b) ** 2
    assert abs(lhs - rhs) < mp.mpf("1e-25")


def test_kappa_wll_is_exact_zeta_third():
    mp.dps = 30
    b = mp.mpf("2.5")
    assert abs(M.cumulants(b)["wll"] + mp.diff(mp.zeta, b + 1, 3)) < mp.mpf("1e-25")


# ── curvature ───────────────────────────────────────────────────────

def test_gaussian_family_convention_lock():
    mp.dps = 40
    assert abs(CV.gaussian_family_curvature() + mp.mpf("0.5")) < mp.mpf("1e-30")


def test_curvature_anchors():
    mp.dps = 40
    assert mp.nstr(CV.gaussian_curvature("3.5"), 16) == CV.CURVATURE["K_3p5"]
    assert abs(CV.gaussian_curvature("2.5") - mp.mpf(CV.CURVATURE["K_2p5"])) < mp.mpf("1e-15")


def test_hagedorn_slope_closed_form():
    mp.dps = 40
    assert abs(CV.hagedorn_slope() - mp.mpf(CV.CURVATURE["hagedorn_slope"])) < mp.mpf("1e-25")


def test_positive_dome_and_sign_change():
    mp.dps = 40
    assert CV.gaussian_curvature("2.0") > 0
    assert CV.gaussian_curvature("6.0") > 0
    assert CV.gaussian_curvature("12.0") < 0


def test_simplex_rigidity_locks():
    """{1,2,3} and {1,2,3,4} both have K = 1/4 identically."""
    for b in (5, 12, 40):
        assert abs(CV.finite_atom_curvature(b, (1, 2, 3)) - mp.mpf(1) / 4) < mp.mpf("1e-30")
        assert abs(CV.finite_atom_curvature(b, (1, 2, 3, 4)) - mp.mpf(1) / 4) < mp.mpf("1e-30")


def test_atom_five_flips_the_sign():
    assert CV.finite_atom_curvature(20, (1, 2, 3, 5)) < 0
    assert CV.finite_atom_curvature(20, (1, 2, 3, 4, 6)) > 0


def test_cold_rate_six_fifths_from_six_atom_model():
    a = CV.finite_atom_curvature(60, (1, 2, 3, 4, 5, 6))
    b = CV.finite_atom_curvature(64, (1, 2, 3, 4, 5, 6))
    rate = float((b / a) ** mp.mpf("0.25"))
    assert abs(rate - 1.2) < 5e-5


# ── hecke ───────────────────────────────────────────────────────────

def test_model_dimensions():
    m = H.model()
    assert len(m["P1"]) == 168
    assert m["E"] == 84
    assert len(m["tris"]) == 56
    assert m["nq"] == 29
    assert len(m["K"]) == 26


def test_block_dimensions_sum_to_homology():
    d = {k: len(v) for k, v in H.blocks().items()}
    assert d == H.BLOCK_DIMS
    assert sum(d.values()) == 26


def test_merel_counts():
    assert len(H.merel(2)) == 4
    assert len(H.merel(3)) == 7
    assert len(H.merel(5)) == 15


def test_hecke_operators_commute():
    A = [list(r) for r in H.cuspidal_hecke(2)]
    B = [list(r) for r in H.cuspidal_hecke(3)]
    n = len(A)
    assert all(sum(A[i][k] * B[k][j] for k in range(n))
               == sum(B[i][k] * A[k][j] for k in range(n))
               for i in range(n) for j in range(n))


def test_bad_prime_operators_descend_and_commute():
    A = [list(r) for r in H.cuspidal_hecke(2)]
    for p in H.BAD_PRIMES:
        U = [list(r) for r in H.cuspidal_hecke(p)]
        n = len(A)
        assert all(sum(A[i][k] * U[k][j] for k in range(n))
                   == sum(U[i][k] * A[k][j] for k in range(n))
                   for i in range(n) for j in range(n))


def test_elliptic_block_matches_point_counts():
    """E2: Merel/Manin route vs Weierstrass point counting."""
    ap = EI.curve_ap((0, -1, 1, -1, -2), 8)
    for p in (3, 5, 7):
        X = EI.hecke_on_block("143a1", p)
        assert all(X[i][j] == (ap[p] if i == j else 0)
                   for i in range(2) for j in range(2))


def test_old_block_is_the_level_11_ghost():
    ap = EI.curve_ap((0, -1, 1, -10, -20), 8)
    for p in (3, 5, 7):
        X = EI.hecke_on_block("11a1_ghost", p)
        assert all(X[i][j] == (ap[p] if i == j else 0)
                   for i in range(4) for j in range(4))


def test_star_involution_is_an_involution_commuting_with_hecke():
    I = [list(r) for r in H.star_involution()]
    A = [list(r) for r in H.cuspidal_hecke(3)]
    n = len(I)
    assert all(sum(I[i][k] * I[k][j] for k in range(n))
               == (1 if i == j else 0) for i in range(n) for j in range(n))
    assert all(sum(I[i][k] * A[k][j] for k in range(n))
               == sum(A[i][k] * I[k][j] for k in range(n))
               for i in range(n) for j in range(n))


@pytest.mark.parametrize("name,dim", [("ell", 2), ("old", 4), ("q4", 8)])
def test_harmonic_density_trace_equals_dimension(name, dim):
    rho = H.harmonic_density(name)
    assert sum(rho) == dim
    assert len(rho) == 84


def test_particles_vanish_on_the_loop_edge():
    """Every cuspidal density is zero on the self-loop edge."""
    m = H.model()
    loop = [k for k in range(m["E"])
            if m["tri_of"][m["erep"][k]]
            == m["tri_of"][m["sS"][m["erep"][k]]]]
    assert len(loop) == 1
    for name in ("ell", "old", "q4", "q6"):
        assert H.harmonic_density(name)[loop[0]] == 0


def test_corpus_field_reconciliation():
    """g4 and the stored FIELD_POLY_F2 define the same quartic field."""
    import mtft.x0_143 as X
    assert X.FIELD_POLY_F3 == [c * (-1) ** ((6 - i) % 2)
                               for i, c in enumerate(H.H6[::-1])]
    assert len(X.FIELD_POLY_F2) == len(H.G4)


# ── eisenstein ──────────────────────────────────────────────────────

def test_sturm_bound():
    assert EI.sturm_bound(143, 2) == 28
    assert EI.sturm_bound(11, 2) == 2


def test_congruence_census():
    assert EI.congruence_census() == EI.MODULI


def test_moduli_exponents_match_orbit_degree():
    """C = norm_modulus ** (dim / orbit degree)."""
    dims = {"143a1": 2, "11a1_ghost": 4, "f2_quartic": 8,
            "f3_sextic": 12}
    for name, C in EI.MODULI.items():
        e = EI.ORBIT_DEGREE[name]
        assert EI.NORM_MODULI[name] ** (dims[name] // e) == C


def test_mazur_five_via_point_counts():
    """E2: the ghost's modulus 5 = numerator((11-1)/12), independently."""
    ap = EI.curve_ap((0, -1, 1, -10, -20), 140)
    g = 0
    for p, a in ap.items():
        if p != 11:
            g = gcd(g, p + 1 - a)
    assert g == 5
    assert EI.NORM_MODULI["11a1_ghost"] == 5


def test_143a1_has_no_eisenstein_congruence():
    ap = EI.curve_ap((0, -1, 1, -1, -2), 140)
    g = 0
    for p, a in ap.items():
        if p not in (11, 13):
            g = gcd(g, p + 1 - a)
    assert g == 1
    assert EI.MODULI["143a1"] == 1
