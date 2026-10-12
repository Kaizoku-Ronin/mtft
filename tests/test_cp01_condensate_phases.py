"""CP-01: spin-1 phases of the condensate alignment, the tree-level preference, the sign-flip test."""
import numpy as np
from mtft.matter import condensate_phases as CP


def test_spin1_invariants_and_phases():
    inv = CP.fit_invariants(); assert abs(inv["a"] - 0.2) < 1e-12 and abs(inv["b"] + 1 / 15) < 1e-12 and inv["max_error"] < 1e-12
    q = lambda c: CP.quartic(c) / CP.norm2(c) ** 2
    assert abs(q(CP.coherent(0.4)) - 1.8) < 1e-12 and abs(q(CP.polar(0.4)) - 1.2) < 1e-12
    rng = np.random.default_rng(3)
    assert all(1.2 - 1e-12 <= q(rng.normal(size=3) + 1j * rng.normal(size=3)) <= 1.8 + 1e-12 for _ in range(500))


def test_texture_of_each_phase():
    d1, d2 = CP.D_from_alignment(CP.coherent(0.1)); assert abs(d2 / d1 ** 2 - 0.5) < 1e-9
    assert CP.D_from_alignment(CP.polar(0.0)) == (0.0, 0.0)
    m2q, m1q = CP.hankel_masses(0.02, 0.0002); m2l, m1l = CP.hankel_masses(0.02, -0.0002)
    assert abs(m2l / m2q - 3) < 0.01 and abs(m1l / m1q - 1 / 3) < 0.01                     # rho = 1/2: Georgi–Jarlskog in the small-D1 limit


def test_sign_flip_against_data():
    g = CP.gj_test()
    assert 0.55 < g["rho_down"] < 0.62 and 3.4 < g["pred_mmu_over_ms"] < 3.7 and 0.27 < g["pred_me_over_md"] < 0.30
    assert abs(g["pred_mmu_over_ms"] * g["pred_me_over_md"] - 1) < 0.02 and 0.85 < g["data_product"] < 0.95


def test_sector_phases_are_a_zero_flip_artefact():                                              # CC-40, CC-61
    from mtft.matter import leptons as LP
    import math
    sols = LP.all_alignments(tuple(LP.LEPTON_POLE))
    reading = lambda s: abs(CP.theta((s[0], s[1] / math.sqrt(2), s[2]))) / (s[0] ** 2 + s[1] ** 2 / 2 + s[2] ** 2)
    assert max(reading(s) for s in sols) > 0.9 and min(reading(s) for s in sols) < 0.1               # the same masses read as either phase
    s = next(x for x in sols if 3 < x[1] < 4); f = CP.zero_flip(s)
    assert np.allclose(CP.raw_masses(s), CP.raw_masses(f), rtol=1e-9) and abs(f[1] - 0.205) < 0.005 and abs(f[2] + 0.028) < 0.002
    assert any(abs(x[1] - f[1]) < 1e-3 and abs(x[2] - f[2]) < 1e-3 for x in sols)                    # the flip is itself in the solution set
    assert set(CP.sector_phases()) == {"up", "down", "leptons"}


def test_rho_branches_and_geometric_masses():
    from mtft.matter.leptons import UP_MZ, DOWN_MZ, LEPTON_POLE
    b = {k: CP.rho_branches(*m) for k, m in (("up", UP_MZ), ("down", DOWN_MZ), ("leptons", LEPTON_POLE))}
    assert abs(b["up"][0] - 0.453) < 0.002 and abs(b["up"][1] + 4.80) < 0.02 and abs(b["down"][0] - 0.577) < 0.002 and abs(b["down"][1] - 3.73) < 0.01
    assert abs(b["leptons"][0] - 0.302) < 0.002 and abs(b["leptons"][1] + 0.765) < 0.002
    g = CP.geometric_test(); assert abs(g["up"] - 0.567) < 0.002 and abs(g["down"] - 2.553) < 0.002 and abs(g["leptons"] - 0.0813) < 0.0005


def test_landau_ground_state_is_parity_even():
    for kappa in (3.0, 1.5, 0.5):
        c = CP.landau_ground(kappa); assert abs(c[1]) < 1e-6
    assert abs(CP.landau_ground(3.0)[0] - 1) < 1e-9 and abs(CP.landau_ground(1.0)[2]) > 0.4


def test_zeeman_vacuum_is_coherent_with_rho_one_half():
    for beta in (0.1, 0.29, 0.6):
        c, e = CP.zeeman_ground(beta); d1, rho = CP.texture_invariants(c)
        assert abs(rho - 0.5) < 1e-5 and abs(e + 1) < 1e-9 and abs(CP.theta(c)) < 1e-6
    assert abs(CP.texture_invariants(CP.zeeman_ground(0.29)[0])[0] - 0.2022) < 0.001


def test_quartic_threshold_and_runaway():
    c, e = CP.zeeman_ground(0.29, 0.45); assert abs(CP.theta(c)) < 1e-5 and abs(e + 1) < 1e-9          # still exactly coherent below 1/2
    c, e = CP.zeeman_ground(0.29, 0.55); assert abs(CP.theta(c)) > 0.3 and abs(CP.texture_invariants(c)[1] - 0.5) > 2   # past 1/2: rho runs away
    d = CP.rho_departures(); assert abs(d["up"] + 0.047) < 0.003 and abs(d["down"] - 0.077) < 0.003 and d["leptons"] < -0.19
