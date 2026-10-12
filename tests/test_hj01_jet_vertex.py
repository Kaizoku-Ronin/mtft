"""HJ-01: the family vertex of the jet classes — selection rule, the round blow-up, profiles."""
import numpy as np
from mtft.matter import jet_vertex as J


def test_selection_rule_forbids_holomorphic_messengers():
    assert all(J.curve_charge(n, b) > 0 for n in range(5) for b in range(20))
    assert J.report()["min_curve_charge"] == 1


def test_round_blowup_has_no_hierarchy():
    rng = np.random.default_rng(7)
    for _ in range(6):
        v = rng.normal(size=2) + 1j * rng.normal(size=2); v /= np.linalg.norm(v)
        assert np.allclose(J.singular_ratios(J.round_texture(*v)), [2.0, 1.0, 1.0], atol=1e-9)


def test_profiles():
    a = J.newform_coefficients(); assert a[1:7] == [1, 0, -1, -2, -1, 0]
    p = J.profiles(); assert abs(p["newform"][1] - 0.25) < 1e-12 and abs(p["bergman_peaked"][1] - 0.4306) < 2e-4
    Y = J.taylor_texture([1, 1, 1, 1, 1], [1, 1, 1, 1, 1]); assert Y.shape == (3, 3) and abs(Y[1, 1] - 2) < 1e-12
