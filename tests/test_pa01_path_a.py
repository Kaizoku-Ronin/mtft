"""PA-01: Path A — content of the matter-untwisted entry, one-loop stability of the untwisting Wilson line, messenger masses."""
from mtft.matter import path_a as PA


def test_content_and_charges():
    c = PA.content()
    assert c["net"] == 3 and c["pairs"] == 25 and c["Rbar"]["jets"] == 2
    assert PA.T_charges() == {"adj": 0, "(4,2,1)": 2, "(4bar,1,2)": 2, "(1,2,2)": 0}       # only matter is charged under T


def test_untwisted_point_is_a_maximum():
    for eps in ((0, 1), (1, 0)):
        prof = PA.potential_profile(eps=eps)
        assert prof[0][1] == -1.0 and abs(prof[-1][1] - 1.0) < 1e-12
        assert all(prof[k][1] < prof[k + 1][1] for k in range(len(prof) - 1))                 # monotone from the jet vacuum to the untwisted point
        assert all(abs(prof[k][1] + prof[-1 - k][1]) < 1e-12 for k in range(len(prof)))      # V(1 - s) = -V(s): odd windings only
        h = PA.hessian_sign(eps=eps); assert h["at_jet_vacuum"] > 0 > h["at_untwisted_point"] and abs(h["at_jet_vacuum"] + h["at_untwisted_point"]) < 1e-9


def test_messenger_masses():
    m = PA.messenger_masses()
    assert abs(m["real (P = 1/2)"] - 3.1057) < 1e-3 and abs(m["complex pair (P = tau/2)"] - 3.5370) < 1e-3
