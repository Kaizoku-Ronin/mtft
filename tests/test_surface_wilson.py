"""WL-01: Wilson-line lifting of the O-block pair; twisted ground state equals the Hodge quadratic form."""
import numpy as np
import pytest

from mtft.surface import spectral as SP


def test_untwisted_limit_and_hermiticity_fast():
    lam0 = SP.wilson_spectrum(11, np.zeros(2), q=1, h=0.25, nx=6, k=3)
    assert abs(lam0[0]) < 1e-6                         # theta = 0: constant zero mode survives
    lam = SP.wilson_spectrum(11, np.array([0.3, 0.0]), q=1, h=0.25, nx=6, k=3)
    assert lam[0] > 1e-3                              # any nontrivial Wilson line lifts it


@pytest.mark.slow
def test_x0143_hodge_prediction_direction_independent():
    e0 = np.eye(26)[0]; e7 = np.eye(26)[7]
    r0 = SP.wilson_line_mass_check(0.05 * e0); r7 = SP.wilson_line_mass_check(0.05 * e7)
    assert abs(r0["ratio"] - 1) < 0.03 and abs(r7["ratio"] - 1) < 0.03
    assert abs(r0["ratio"] - r7["ratio"]) < 0.005     # the SHAPE of the mass form is the frozen Hodge metric
    r3 = SP.wilson_line_mass_check(0.05 * e0, q=3)
    assert abs(r3["lambda_0"] / r0["lambda_0"] - 9) < 0.1
