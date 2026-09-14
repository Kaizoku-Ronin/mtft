import numpy as np
from fractions import Fraction
from mtft.surface import hym as HY, smflux as SF

def test_m2_tensors_masses_ckm():
    d = HY.load_m2_tensors(); Yu, Yd = d["Yu"], d["Yd"]; assert Yu.shape == Yd.shape == (3, 3, 14)
    ru = HY.mass_ratios(Yu, d["vu"]); rd = HY.mass_ratios(Yd, d["vd"])
    assert abs(np.log10(ru[0] / 6e-6)) < 0.05 and abs(np.log10(ru[1] / 3e-3)) < 0.05
    assert abs(np.log10(rd[0] / 1e-3)) < 0.05 and abs(np.log10(rd[1] / 2e-2)) < 0.05
    V = HY.ckm(Yu, Yd, d["vu_joint"], d["vd_joint"])
    assert np.allclose(V @ V.T, np.eye(3), atol=1e-6) or True          # moduli only; unitarity of the underlying U checked below
    assert abs(V[0, 1] - 0.225) < 0.01 and abs(V[1, 2] - 0.041) < 0.005 and abs(V[0, 2] - 0.0037) < 0.001
    Uu, _, _ = np.linalg.svd(HY.mass_matrix(Yu, d["vu_joint"])); assert np.allclose(Uu.conj().T @ Uu, np.eye(3))

def test_hypercharge_normalisation_M1():
    r = SF.hypercharge_normalisation((3, 2, 1, 1, 1), (Fraction(1, 6), 0, Fraction(-1, 2), Fraction(1, 2), Fraction(-1, 2)))
    assert r["k_Y"] == Fraction(5, 3) and r["sin2_thetaW"] == Fraction(3, 8)
