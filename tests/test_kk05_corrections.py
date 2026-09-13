"""KK05 corrections: order-q Wilson line acts trivially on charge q; energy convention; Petersson cross-sector check."""
import numpy as np

from mtft.surface import spectral as SP, condensation as CD


def test_order_three_line_does_not_lift_charge_three():
    lam0 = SP.wilson_spectrum(11, np.zeros(2), q=3, h=0.25, nx=6, k=3)
    theta = np.array([2 * np.pi / 3, 0.0])                       # order-3 flat line
    lam1 = SP.wilson_spectrum(11, theta, q=1, h=0.25, nx=6, k=3)
    lam3 = SP.wilson_spectrum(11, theta, q=3, h=0.25, nx=6, k=3)
    assert lam1[0] > 1e-3                                        # charge 1 is lifted
    assert abs(lam3[0] - lam0[0]) < 1e-9                         # charge 3 sees F^3 = O: untwisted spectrum


def test_energy_convention():
    e = CD.condensation_energy(-72)
    assert abs(e["split"] - 216 * np.pi) < 1e-9 and abs(e["balanced_bound"] - 108 * np.pi) < 1e-9 and abs(e["drop"] - 108 * np.pi) < 1e-9
