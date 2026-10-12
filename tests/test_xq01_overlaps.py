"""XQ-01: the overlaps X_qq' — family normalisation, cusp Bergman jets of X_0(143), charge-conserving weights."""
import numpy as np
from mtft.matter import overlaps as X


def test_cusp_bergman_jets():
    B = X.cusp_bergman_jets()
    assert np.allclose(B, B.T, atol=1e-12) and np.all(np.linalg.eigvalsh(B) > 0)
    assert np.allclose(np.diag(B), [1.0, 2.0941, 2.9255, 4.1597, 4.9221, 5.2081], atol=2e-4)
    C = X.correlations(B); off = np.abs(C - np.diag(np.diag(C)))
    assert off.max() < 0.32 and abs(C[3, 1]) < 0.002                       # nearly conserves the cusp rotation charge


def test_charge_conserving_weights():
    w = X.antidiagonal_weights(X.cusp_bergman_jets())
    assert np.allclose(w[:3], [0.31618, 0.00540, 2.92550], atol=2e-4)
    assert 150 < np.sqrt(w[0] * w[2]) / abs(w[1]) < 200


def test_normalisation_convention():                                        # CC-41
    assert abs(X.geometric_ratio(0.02, 0.0002, True) - 0.25) < 0.002 and abs(X.geometric_ratio(0.0, 0.0002, True) - 0.25) < 1e-6
    assert abs(X.geometric_ratio(0.02, 0.0002, False) - 1.0) < 1e-6
    r = X.data_rho(True); assert abs(r["up"] - 0.568) < 0.002 and abs(r["down"] - 0.684) < 0.002 and abs(r["leptons"] - 0.408) < 0.002
