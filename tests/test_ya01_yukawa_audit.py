"""YA-01: character audit of the Yukawa amplitude for jet families."""
import numpy as np
from mtft.matter import yukawa_audit as Y


def test_holomorphic_higgs_amplitudes_vanish():
    assert Y.forbidden_holomorphic() and min(Y.holomorphic_amplitude_character(m, -2) for m in range(10)) == 4


def test_curve_charge_alone_allows_only_n_three_and_four():
    assert sorted(Y.curve_allowed()) == [(3, 0, 2), (4, 0, 1), (4, 1, 2)]
    assert abs(Y.scale_ratio() - 0.5774) < 1e-4


def test_leading_texture():
    assert np.abs(Y.leading_texture(0.0)).max() < 1e-15
    for t in (0.06, 0.1, 0.14):
        m = Y.masses(Y.leading_texture(t)); assert m[2] < 1e-12 and abs(m[1] / t ** 2 - 1) < 0.08
