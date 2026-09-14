"""SM-02 step 5: the three-family Riemann–Roch space constructed and certified (sub-second)."""
import numpy as np
import pytest

from mtft.surface import rrspace as RR


def test_three_family_sections_certified():
    r = RR.three_family_sections()
    assert r["f1_max_rel_at_Z"] < 1e-8                       # div(f1 dz) contains all 20 W143 points
    assert r["rank"] == 33 and r["dimension"] == 3           # 33 independent conditions on 36 unknowns
    GP = r["G_at_W13_classes"]                               # rows: P4, other, W11 P4, W11 other
    assert GP[0].max() < 1e-6                                # vanishes at the excluded point
    assert GP.max(axis=1).max() == 1.0 and (GP.max(axis=1) > 1e-3).sum() == 2   # poles at exactly two points (the same-sign pair)


def test_higgs_target_and_up_yukawa_M1():
    r = RR.up_yukawa_M1()
    assert r["target"]["rank"] == 42 and r["target"]["dimension"] == 18
    assert r["residual"] < 1e-9                                                  # G_i G_j in f_K * H^0(K(2 SigmaP))
    Y = r["Y"]; assert np.allclose(Y, np.transpose(Y, (1, 0, 2)))
    rng = np.random.default_rng(0)
    assert all(np.linalg.matrix_rank(np.einsum("ijk,k->ij", Y, rng.standard_normal(18)), tol=1e-8) == 3 for _ in range(3))
    fK = RR.f_K_series(20); assert fK[0] == 0 and fK[1] == 1                    # f_K = q - 2q^2 - q^3 + ...


def test_w13_texture_of_the_yukawa_tensor():
    r = RR.up_yukawa_M1(); g = RR.w13_grading_and_texture(r)
    assert abs(g["T_squared"] + 13) < 1e-3 and g["closure_residual"] < 1e-6
    assert abs(sum(g["grades_sections"])) == 1 and sorted(g["grades_higgs"]).count(1) == 8   # 1 + 2 split (singleton grade is a branch convention); 8 even / 10 odd
    assert g["max_on_forbidden"] < 1e-6 and g["max_on_allowed"] > 0.5      # selection rule grade_i grade_j grade_k = -1
    assert g["rank_even_higgs"] == 2 and g["rank_odd_higgs"] == 3          # even VEV: one massless family


def test_twisted_sector_and_down_yukawa():
    r = RR.down_yukawa_M1()
    tw = r["twisted"]
    assert (tw[4]["rank"], tw[4]["dimension"]) == (37, 3) and (tw[6]["rank"], tw[6]["dimension"]) == (50, 18)
    assert r["residual"] < 1e-7                                       # Q x d^c products lie in f_K x twisted Higgs space
    rng = np.random.default_rng(0)
    assert all(np.linalg.matrix_rank(np.einsum("ijk,k->ij", r["Y"], rng.standard_normal(18)), tol=1e-8) == 3 for _ in range(3))


def test_cubic_twist_lepton_tensor_lifts_degeneracy():
    lep = RR.lepton_yukawa_M1(); assert lep["L_space"]["dimension"] == 3 and lep["E_space"]["dimension"] == 3
    assert lep["residual"] < 1e-6                                     # chi_133 * chi_56 = chi_13: products land in the chi_13 Higgs space
    Yd = RR.down_yukawa_M1()["Y"]; Yl = lep["Y"]
    assert RR.family_equivalence_residual(Yd, Yd) < 1e-8             # control: self-equivalence
    assert RR.family_equivalence_residual(Yl, Yd) > 1e-3 and RR.family_equivalence_residual(Yl, np.transpose(Yd, (1, 0, 2))) > 1e-3
