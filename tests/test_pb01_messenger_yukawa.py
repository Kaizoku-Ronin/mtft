"""PB-01: the holomorphic messenger route closes; the D-term vacuum has one orientation; the cusp kernel of weight-2 forms."""
import numpy as np
import pytest
from mtft.matter import messenger_yukawa as MY


def test_su_invariants_and_holomorphic_orders():
    assert (MY.su_invariants(2, 2, 0), MY.su_invariants(3, 3, 0), MY.su_invariants(4, 1, 1), MY.su_invariants(4, 4, 0)) == (1, 1, 1, 1)
    h = MY.holomorphic_orders()
    assert h["Sigma^0"] == 1 and h["Sigma^1"] == 0 and h["Sigma^2"] == 0 and h["Sigma^3"] == 0 and h["Sigma Sigmabar"] == 2 and h["Sigma^4"] == 4


def test_d_term_vacuum_has_one_orientation():
    v = MY.d_term_vacuum(); s = MY.d_term_vacuum_scan(n=81)
    assert v["Sigmabar_sq"] == 0.0 and s["Sigmabar_sq"] == 0.0 and abs(v["V"] - s["V"]) < 1e-3 and v["residual_nonabelian_D"] > 0


@pytest.fixture(scope="module")
def ck():
    return MY.cusp_kernel()


def test_cusp_kernel(ck):
    N = ck["sector_weights"]
    assert ck["sector_sizes"] == {(1, 1): 1, (1, -1): 6, (-1, 1): 5, (-1, -1): 1}
    assert abs(N[(1, 1)] - 2.4855) < 2e-3 and abs(N[(1, -1)] - 4.4384) < 2e-3 and abs(N[(-1, 1)] - 4.5016) < 2e-3 and abs(N[(-1, -1)] - 1.0661) < 2e-3
    ev = ck["eigenvalues_four_cusps"]; assert np.allclose(sorted(ev), sorted(4 * np.array(list(N.values()))), atol=1e-9)   # circulant: eigenvalues 4 N_chi
    three = ck["eigenvalues_three_cusps"]; vals = list(three.values())
    assert all(np.allclose(v, vals[0], atol=1e-9) for v in vals) and np.allclose(vals[0], [6.4375, 13.1556, 17.8817], atol=2e-3)
