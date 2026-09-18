import numpy as np
from fractions import Fraction
from mtft.surface import hym as HY, smflux as SF

def test_cc26_normalisation_is_basis_invariant_and_transpose_is_not():
    d = HY.load_m2_tensors()
    for Y, NA, NB, NH in ((d["Y_up_raw"], d["N_Q"], d["N_u"], d["N_Hu"]), (d["Y_down_raw"], d["N_L"], d["N_L"], d["N_Hd"])):
        assert HY.normalisation_is_basis_invariant(Y, NA, NB, NH)
        chol = lambda N: np.linalg.cholesky((N + N.conj().T) / 2); A = np.linalg.inv(chol(NA))
        assert np.allclose(A @ NA @ A.conj().T, np.eye(3)) and not np.allclose(A.T @ NA @ A.T.conj().T, np.eye(3))
    Yu = HY.normalise_yukawa(d["Y_up_raw"], d["N_Q"], d["N_u"], d["N_Hu"]); assert np.allclose(Yu, d["Yu"])

def test_m2_corrected_ratios_are_order_one():
    d = HY.load_m2_tensors(); rng = np.random.default_rng(0); R = []
    for _ in range(300):
        v = rng.standard_normal(14) + 1j * rng.standard_normal(14); R.append(HY.mass_ratios(d["Yu"], v))
    R = np.array(R); assert 0.3 < np.median(R[:, 1]) < 0.9 and 0.05 < np.median(R[:, 0]) < 0.6     # CC-26: no hierarchy at leading order

def test_hypercharge_normalisation_M1():
    r = SF.hypercharge_normalisation((3, 2, 1, 1, 1), (Fraction(1, 6), 0, Fraction(-1, 2), Fraction(1, 2), Fraction(-1, 2)))
    assert r["k_Y"] == Fraction(5, 3) and r["sin2_thetaW"] == Fraction(3, 8)
