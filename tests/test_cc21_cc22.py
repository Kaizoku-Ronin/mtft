"""CC-21 random_su_n group membership; CC-22 Higgs quartic convention accessor."""
import numpy as np

import mtft.constants as C
from mtft import lattice as L


def test_cc21_random_su_n_is_in_su_n():
    rng = np.random.default_rng(21)
    for N in (2, 3, 5):
        for _ in range(16):
            U = L.random_su_n(N, rng)
            assert abs(np.linalg.det(U) - 1) < 1e-10
            assert np.linalg.norm(U.conj().T @ U - np.eye(N)) < 1e-10
    Ud = L.random_su_n_defective_v0261(3, np.random.default_rng(0))
    assert abs(np.linalg.det(Ud) - 1) > 1e-3            # the retired initializer really was defective


def test_cc22_lambda_convention():
    lam = C.HIGGS.lambda_quartic
    assert abs(lam - 0.5179175466708401) < 1e-12          # formula unchanged
    assert abs(C.HIGGS.lambda_quartic_pdg - lam / 4) < 1e-15
    assert abs(C.HIGGS.m_H ** 2 - lam * C.HIGGS.V_EW ** 2 / 2) < 1e-6 * C.HIGGS.m_H ** 2   # m_H^2 = lam v^2/2
