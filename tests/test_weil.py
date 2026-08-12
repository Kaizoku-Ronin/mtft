"""Tests for the candidate mtft.weil module (W1). Fast tier runs in seconds;
the E2 test uses the cached zero file zeros_gamma_T100.npy when present
(371 ordinates, gamma <= 640.7, mpmath dps 15) and skips otherwise.
Sympy-free by house rule."""
import os
import numpy as np
import pytest
from mtft import weil


def test_w_series_cc02():
    r = weil.w_series_check(s=3, N=150_000)
    assert r["diff_correct"] < 1e-9            # -zeta(s)zeta'(s+1) is the series
    assert r["diff_paper1"] > 1e-3             # printed Paper 1 value excluded
    assert r["diff_ag"] > 1e-3                 # printed AG value excluded


def test_w_convolution_e2():
    # w = Lambda_1 * 1 (Lambda_1 = log/id) equals the direct divisor sieve
    N = 20_000
    w1 = np.zeros(N + 1)
    for dd in range(2, N + 1):
        w1[dd::dd] += np.log(dd) / dd
    lam = np.zeros(N + 1)
    comp = np.zeros(N + 1, dtype=bool)
    for p in range(2, N + 1):
        if not comp[p]:
            comp[p * p::p] = True
            pk = p
            while pk <= N:
                lam[pk] = np.log(p)
                pk *= p
    sig = np.zeros(N + 1)
    for dd in range(1, N + 1):
        sig[dd::dd] += 1.0 / dd
    w2 = np.zeros(N + 1)
    for e in range(2, N + 1):
        if lam[e] > 0:
            m = np.arange(1, N // e + 1)
            w2[e * m] += (lam[e] / e) * sig[m]
    assert np.max(np.abs(w1[1:] - w2[1:])) < 1e-12


def test_rank_trace_random_and_equality():
    a = weil.audit_rank_trace(trials=20_000, seed=143)
    assert a["violations"] == 0
    d, r, b, c = 10, 3, 4, 2.0
    P = np.zeros((d, d)); P[:r, :r] = (c / 2) * np.eye(r)
    Q = np.zeros((d, d)); Q[r:r + b, r:r + b] = c * np.eye(b)
    assert abs(weil.rank_trace_gap(P, Q, r, b, c)) < 1e-12


def test_mt_constant():
    assert abs(weil.mt_constant(1.0) - 0.753296067856) < 1e-9
    assert abs((2 - 1 / weil.mt_constant(1.0)) - 0.672500703679) < 1e-9


def test_taper_constants_and_poisson():
    g = weil.gabor(100.0, 1.0)
    win = weil.Window(g["L"], eta=0.2, n_u=3001)
    assert abs(win.a - (1 - 0.603 * 0.2)) < 2e-4
    assert abs(win.b - (1 - 0.688 * 0.2)) < 2e-4
    # Lemma 2.2 at modest lattice truncation
    K = 1500
    kk = np.arange(-K, K + g["d"])
    tk = 100.0 + g["h"] * kk
    rng = np.random.default_rng(7)
    for _ in range(3):
        ta, tb = rng.uniform(90, 210, 2)
        lhs = float(np.sum(win.phihat(ta - tk) * win.phihat(tb - tk)))
        rhs = float(g["L"] * win.Phi(ta - tb))
        assert abs(lhs - rhs) / abs(rhs) < 5e-4


def test_synthetic_inertia_inequalities():
    g = weil.gabor(100.0, 1.0)
    win = weil.Window(g["L"], eta=0.2, n_u=3001)
    rng = np.random.default_rng(11)
    ords = np.sort(rng.uniform(95, 205, 30))
    for deps in (np.full(30, 0.2),
                 np.where(np.arange(30) % 2 == 0, 0.05, 0.0),
                 rng.uniform(1e-4, 0.45, 30) * (rng.random(30) < 0.5)):
        G, cnt = weil.G_zero(ords, 100.0, eta=0.2, depths=deps, mirror=False, win=win)
        ev = np.linalg.eigvalsh(G / g["L"])
        assert int((ev > 1e-9).sum()) <= cnt["s1"] + cnt["s2"] + cnt["p"]
        Ah = G / (win.a * g["L"] ** 2)
        cert = 4 * np.trace(Ah) - 2 * cnt["N"] - np.sum(Ah * Ah)
        assert cert <= cnt["s1"] + 1e-9


def test_e2_prime_vs_zero_cached():
    path = os.path.join(os.path.dirname(__file__), "zeros_gamma_T100.npy")
    if not os.path.exists(path):
        pytest.skip("zero cache not present; run w1_study_driver.py to build it")
    gammas = np.load(path)
    Gp, g, win = weil.G_prime(100.0, 1.0, eta=0.2, pad=250.0, dtau=0.01)
    Gz, _ = weil.G_zero(gammas, 100.0, 1.0, eta=0.2, win=win)
    rel = np.max(np.abs(Gp - Gz)) / np.max(np.abs(Gp))
    assert rel < 5e-5      # certified 3.4e-6 at dtau=0.005; slack for coarse grid
