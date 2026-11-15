"""CC-72 (2026-11-15): the Burning Ship fold exponents are a caveat, not a measurement.

Paper 4's off-axis table (periods 48 to 768) is not one period-doubling cascade, so beta_x and
beta_y, fitted to it, are unconfirmed.  The constants stay exposed (append-only).  The control is
the real-axis cascade of x^2 + c, which does shrink by Feigenbaum's delta.
"""
import math

import pytest

from mtft import burning_ship as BS
from mtft import constants as K


def test_paper4_table_is_not_one_cascade():
    r = BS.cc72_cascade_check()
    assert r["gaps"] == (2.42e-4, 1.81e-4, 2.08e-4, 1.34e-4)
    assert abs(r["shrink"] - 1.806) < 0.01                      # 1.81 times over three doublings
    assert abs(r["cascade_shrink"] - 101.79) < 0.01             # a cascade: delta^3
    assert r["cascade_shrink"] / r["shrink"] > 50               # not a matter of tolerance
    assert r["ratios"][1] < 1.0 and not r["monotone"]           # the third gap is larger than the second
    assert r["one_cascade"] is False


def test_control_a_real_cascade_passes_the_same_criterion():
    c = BS.real_axis_cascade(10)
    assert c["c"][0] == 0.0 and c["c"][1] == -1.0
    assert abs(c["c"][2] - (-1.3107026413368329)) < 1e-12
    assert abs(c["ratios"][-1] - 4.6691951) < 1e-6              # delta_10, as the lab's page has it
    assert abs(c["ratios"][-1] - K.FEIGENBAUM_DELTA) < 1e-5
    assert abs(c["c_inf"] - (-1.401155189092)) < 1e-9
    last4 = c["gaps"][-4:]                                      # four gaps, as in Paper 4's table
    ratios = [last4[i] / last4[i + 1] for i in range(3)]
    assert all(abs(r / K.FEIGENBAUM_DELTA - 1.0) < 0.25 for r in ratios)      # the criterion of cc72_cascade_check
    assert abs(last4[0] / last4[-1] - K.FEIGENBAUM_DELTA ** 3) < 0.01         # shrinks 101.8 times, not 1.81


def test_constants_stay_exposed_unchanged():
    assert K.BETA_X == +0.29 and K.BETA_Y == -0.21
    assert math.isclose(K.DELTA_X, K.FEIGENBAUM_DELTA ** 1.29) and math.isclose(K.DELTA_Y, K.FEIGENBAUM_DELTA ** 0.79)
    assert abs(BS.ANISOTROPIC.fold_budget - 0.5) < 1e-12
    assert "CC-72" in BS.AnisotropicScaling.__doc__
    assert "mtft.constants.BETA_X" in BS.CC72_DEPENDENTS and len(BS.cc72_cascade_check()["dependents"]) == 7


def test_real_axis_cascade_rejects_depths_double_precision_cannot_carry():
    for bad in (1, 12):
        with pytest.raises(ValueError):
            BS.real_axis_cascade(bad)
