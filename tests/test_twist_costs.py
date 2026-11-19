"""TC-01: the twist table, the closed-form Hessian at the Pati–Salam pattern, the continuous twist and the pull on τ."""
import math

import numpy as np
import pytest

from mtft.matter import twist_costs as TC

E2 = TC.E2


def test_complementarity_on_all_nine_pairs():
    rows = TC.twist_rows()
    assert len(rows) == 9 and TC.complementarity()
    for k, r in rows.items():
        same = r["P"] == r["h"]
        assert r["tree Yukawas"] == same and (r["V"] == 0) == same and r["MA = 0"] == same and r["twelve pinned"] == (not same)
        assert r["Higgs bosons untwisted"] and not r["condensate boson untwisted"]
        assert r["signature"] == ((2, 14, 0) if same else (1, 2, 13))


def test_anchors_of_ff13():
    r = TC.twist_row(E2["p2"], E2["p2"])
    assert np.allclose(r["T_S"], (-1335.483488, -231.1462), rtol=1e-6) and r["twelve"] == (0.0, 0.0) and r["V"] == 0
    r = TC.twist_row(E2["p1"], E2["p2"])
    assert math.isclose(r["V"], -140.88499928066574, rel_tol=1e-12)
    assert np.allclose(r["T_S"], (-628.441487, 1033.507557), rtol=1e-6) and np.allclose(r["twelve"], (21.125436, 1395.401655), rtol=1e-6)
    r = TC.twist_row(E2["p2"], E2["p1"])
    assert math.isclose(r["V"], -50.8490106650271, rel_tol=1e-12) and np.allclose(r["twelve"], (49.701113, 505.467554), rtol=1e-6)


def test_closed_hessian_against_differences():
    for P, h in (("p1", "p2"), ("p2", "p1"), ("p2", "p2")):
        err = TC.compare_routes(TC.weights_spin(E2[P]), E2[h])
        assert err < 1e-8
    W = TC.weights_continuous(0.1)
    assert TC.compare_routes(W, E2["p2"]) < 1e-8
    assert TC.compare_routes(W, E2["p2"]) > 0                     # non-vacuity: the two routes are computed apart


def test_continuous_twist_pin_reproduces_ff12_and_masses():
    cv = TC.pattern_curvatures(TC.weights_continuous(0.05), E2["p2"])
    assert np.allclose(cv["twelve"], (0.00820168913707508, 0.06625313929920104), rtol=1e-9)
    m = TC.continuous_masses(0.05, brute=True)
    for k, v in m["closed"].items():
        assert math.isclose(v, m["brute force"][k], rel_tol=1e-12)
    assert math.isclose(m["closed"]["Phi-type Higgs"], 2 * 0.05 * TC.m_ss(), rel_tol=1e-15)
    assert math.isclose(TC.m_ss(), 0.4942809911220414, rel_tol=1e-12)


def test_hold_pull_jet_vacuum_and_twists():
    lo, hi = TC.jet_vacuum_pull()
    assert math.isclose(lo, 8.730665725945551, rel_tol=1e-6) and math.isclose(hi, 10.628636535933715, rel_tol=1e-6)
    same = TC.hold_pull(TC.pattern_energy_of_tau(E2["p2"], E2["p2"]))
    assert same["V"] == 0 and same["max sinh d"] == 0 and same["held"]
    other = TC.hold_pull(TC.pattern_energy_of_tau(E2["p1"], E2["p3"]))
    assert math.isclose(other["min sinh d"], 619.8887162532019, rel_tol=1e-5) and not other["held"]

    def Vc(t):
        return TC.pattern_energy(TC.weights_continuous(0.05, t=complex(t)), E2["p2"], t=complex(t))
    cont = TC.hold_pull(Vc)
    assert math.isclose(cont["max sinh d"], 0.010981504030410701, rel_tol=1e-5) and cont["held"]


def test_block_status_is_exact():
    st = TC.block_status(E2["p2"], E2["p2"])
    assert all(st[b]["fermions untwisted"] for b in TC.MATTER) and all(st[b]["bosons untwisted"] for b in TC.HIGGS)
    st = TC.block_status(E2["p1"], E2["p2"])
    assert not any(st[b]["fermions untwisted"] for b in TC.MATTER)


def test_raise_paths():
    m, n, _ = TC.windings()
    with pytest.raises(ValueError):
        TC.twist_row((0.25, 0), E2["p2"])
    with pytest.raises(ValueError):
        TC.pattern(E2["O"])
    with pytest.raises(ValueError):
        TC.continuous_masses(0.6)
    with pytest.raises(ValueError):
        TC.hold_pull(lambda t: 0.0, window=())
    with pytest.raises(ValueError):
        TC.windings(0)
    with pytest.raises(ValueError):
        TC.pattern_matrices(np.ones(3), E2["p2"])
    with pytest.raises(ValueError):
        TC.pattern_matrices(-np.ones_like(m), E2["p2"])
    with pytest.raises(ValueError):
        TC.weights_continuous(1.5)
    with pytest.raises(ValueError):
        TC.chi((0.25, 0), m, n)
    with pytest.raises(ValueError):
        TC.hessian_by_differences(TC.weights_spin(E2["p1"]), E2["p2"], [np.zeros((8, 2))])
    with pytest.raises(ValueError):
        TC.twist_row((1, 2, 3), E2["p2"])
