"""Tests for the v0.23.0 modules: homology, thetachar, thetafun, liealg.

Runnable under pytest or directly (python tests/test_v0230.py).
Set MTFT_SLOW=1 to additionally recompute the full theta census
(a few minutes) and compare against the packaged certificate.
"""
import os

import numpy as np


def test_homology():
    from mtft import homology as H
    m = H.matrices()                      # all structural asserts fire here
    U, ops = H.symplectic_frame()
    J = H.standard_J()
    assert np.array_equal(U.T @ m["P"] @ U, J)
    per = H.periods_frame_ops()
    J2 = H.mod2(per["J"])
    g = 13
    Jstd = np.zeros((26, 26), np.uint8)
    Jstd[:g, g:] = np.eye(g, dtype=np.uint8)
    Jstd[g:, :g] = np.eye(g, dtype=np.uint8)
    assert np.array_equal(J2, Jstd)
    for k in ("W11", "W13", "STAR"):
        M2 = H.mod2(per[k])
        assert np.array_equal((M2.T @ Jstd @ M2) % 2, Jstd)


def test_thetachar():
    from mtft import thetachar as TC
    ref = TC.census()
    for make in (TC.x0143_periods_frame, TC.x0143_gp_frame):
        act = make()
        TC.verify_action(act, samples=25)
        t0, B = act.joint_fixed_locus()
        assert B.shape[1] == ref["joint"]["affine_dim"] == 7
        ev, od = act.parity_count(t0, B)
        assert (ev, od) == (ref["joint"]["even"], ref["joint"]["odd"]) \
            == (96, 32)
        r, rad = act.parity_polarization(t0, B)
        assert (r, rad) == (2, 5)
        pts, par = act.invariant_characteristics()
        assert pts.shape == (128, 26)
        assert int((par == 0).sum()) == 96 and int((par == 1).sum()) == 32
    if os.environ.get("MTFT_SLOW"):
        act = TC.x0143_periods_frame()
        full = act.full_census()
        assert full["element_fixed"] == ref["element_fixed"]
        assert full["burnside"] == ref["burnside"]
        assert full["orbit_sizes"] == ref["orbit_sizes"]


def _brute(a, b, tau, N):
    import itertools
    g = len(a)
    tot = 0j
    grad = np.zeros(g, complex)
    for n in itertools.product(range(-N, N + 1), repeat=g):
        x = np.array(n, float) + a
        t = np.exp(1j * np.pi * (x @ tau @ x) + 2j * np.pi * (x @ b))
        tot += t
        grad += 2j * np.pi * x * t
    return tot, grad


def test_thetafun_g1():
    import mpmath as mp
    from mtft import thetafun as TF
    tau1 = np.array([[0.3 + 1.1j]])
    rd = TF.siegel_ready(tau1)
    q = mp.e ** (1j * mp.pi * mp.mpc(tau1[0, 0]))
    v, vb = TF.theta_char(*TF.char_to_ab(np.array([0, 0])), rd, tol=1e-12)
    assert abs(v - complex(mp.jtheta(3, 0, q))) < 1e-10
    for t2 in ([0, 1], [1, 0], [1, 1]):
        a, b = TF.char_to_ab(np.array(t2))
        v, vb = TF.theta_char(a, b, rd, tol=1e-12)
        bv, bg = _brute(a, b, tau1, 30)
        assert abs(v - bv) < 1e-10 + vb
    a, b = TF.char_to_ab(np.array([1, 1]))
    gr, v, vb, gb = TF.theta_grad(a, b, rd, tol=1e-12)
    bv, bg = _brute(a, b, tau1, 30)
    assert abs(v - bv) < 1e-9 and np.abs(gr - bg).max() < 1e-9


def test_thetafun_g2():
    from mtft import thetafun as TF
    # diagonal tau: exact factorization into g = 1 thetas
    tauA = np.array([[0.2 + 1.3j, 0], [0, -0.4 + 0.9j]])
    rdA = TF.siegel_ready(tauA)
    for t4 in ([0, 0, 0, 0], [1, 0, 0, 1], [1, 1, 1, 1], [0, 1, 1, 0]):
        a, b = TF.char_to_ab(np.array(t4))
        v, vb = TF.theta_char(a, b, rdA, tol=1e-12)
        prod = 1 + 0j
        for i in (0, 1):
            rd1 = TF.siegel_ready(tauA[i:i + 1, i:i + 1])
            vi, _ = TF.theta_char(a[i:i + 1], b[i:i + 1], rd1, tol=1e-13)
            prod *= vi
        assert abs(v - prod) < 1e-9
    # mildly skewed tau: full reduction path vs direct brute
    tauB = np.array([[0.4 + 1.2j, 0.3 - 0.5j], [0.3 - 0.5j, -0.6 + 1.5j]])
    rdB = TF.siegel_ready(tauB)
    for t4 in ([1, 0, 1, 1], [0, 0, 1, 1], [1, 1, 0, 1]):
        a, b = TF.char_to_ab(np.array(t4))
        gr, v, vb, gb = TF.theta_grad(a, b, rdB, tol=1e-11)
        bv, bg = _brute(a, b, tauB, 12)
        assert abs(v - bv) < 1e-8 + vb
        assert np.abs(gr - bg).max() < 1e-7 + gb


def test_thetafun_split():
    from mtft import thetafun as TF
    rng = np.random.default_rng(5)
    B = rng.normal(0, .35, (8, 8))
    Y = B @ B.T + 1.1 * np.eye(8)
    X = np.round(rng.normal(0, .5, (8, 8)), 2)
    tau = (X + X.T) / 2 + 1j * Y
    rd = TF.siegel_ready(tau)
    odd = np.array([1] + [0] * 7 + [1] + [0] * 7)
    even = np.zeros(16, np.int64)
    for tv in (even, odd):
        a, b = TF.char_to_ab(tv)
        a2, b2, _ = TF.reduce_char(a, b, rd)
        s1 = TF.theta_reduced(a2, b2, rd, tol=1e-10, deriv=True,
                              split=4, stats=True)
        s0 = TF.theta_reduced(a2, b2, rd, tol=1e-10, deriv=True,
                              split=0, stats=True)
        assert min(s1["npoints"], s0["npoints"]) > 1000
        assert abs(s1["value"] - s0["value"]) < 1e-9
        assert np.abs(s1["grad"] - s0["grad"]).max() < 1e-8
    a, b = TF.char_to_ab(odd)
    gr, v, vb, gb = TF.theta_grad(a, b, rd, tol=1e-10)
    assert abs(v) < 1e-9 + vb            # odd theta null vanishes
    assert np.linalg.norm(gr) > 1e-2     # with a genuinely nonzero gradient


def test_liealg_d4():
    from mtft import liealg as L
    rep = L.d4_report(screen=True)
    assert rep["dim"] == 28
    assert rep["closure"]["growth"][:4] == [3, 6, 17, 28]
    s = rep["structure"]
    assert (s["center_dim"], s["derived_dim"], s["rank"]) == (0, 28, 4)
    assert s["killing_signature"] == [28, 0, 0]
    r = rep["representation"]
    assert (r["common_fixed_dim"], r["active_dim"]) == (5, 8)
    assert r["active_commutant_dim"] == 1 and r["u_n_commutant_dim"] == 26
    assert rep["normalizer"]["kernel_dim"] == 54
    assert rep["normalizer"]["gap"] > 1e6
    assert rep["roots"]["count"] == 24
    assert rep["roots"]["length_ratio"] < 1 + 1e-10
    assert rep["roots"]["plus_minus_paired"] and rep["roots"]["zero_modes"] == 4
    assert rep["roots_cosine_dev"] < 1e-12
    sc = rep["symmetry_screen"]
    assert sc["STAR"]["kind"] == "antilinear"
    assert sc["STAR"]["identity_deviation_on_g"] < 1e-5
    for k in ("W11", "W13", "W143"):
        assert sc[k]["normalization_residual"] > 0.1


if __name__ == "__main__":
    import time
    for fn in (test_homology, test_thetachar, test_thetafun_g1,
               test_thetafun_g2, test_thetafun_split, test_liealg_d4):
        t0 = time.time()
        fn()
        print("%-22s PASS  %.1fs" % (fn.__name__, time.time() - t0))
    print("ALL v0.23.0 TESTS PASS")
