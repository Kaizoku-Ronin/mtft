"""mtft.maass: the K-Bessel routine, the reduction with its characters, and Hejhal's solver.

Gates of the v1.1.0 plan (Tool 1).  Outside references: LMFDB's Maass form 13.1 (R = 0.97081541)
and the classical first eigenvalues of SL(2, Z) (9.5336952614, 13.7797513519).
"""
import math

import mpmath as mp
import numpy as np
import pytest

from mtft.maass import domain as D
from mtft.maass import hejhal as H
from mtft.maass import kbessel as KB
from mtft.maass import scattering as S


def test_kbessel_fast_route_against_mpmath():
    worst = 0.0
    for R in (0.05, 0.37, 1.0, 2.2, 3.0):
        for x in (0.02, 0.1, 0.7, 2.5, 9.0, 30.0):
            a, b = float(KB.kir(R, np.array([x]))[0]), float(mp.besselk(1j * R, x).real)
            worst = max(worst, abs(a - b) / abs(b))
    assert 0 < worst < 1e-12                                     # the lab's record: 2.1e-14
    assert KB.kir(1.0, np.array([61.0]))[0] == 0.0               # cut off above x = 60
    with pytest.raises(ValueError):
        KB.kir(3.5, np.array([1.0]))                             # outside the validated range
    assert abs(KB.k_bessel(3.5, np.array([1.0]))[0] - float(mp.besselk(3.5j, 1.0).real)) < 1e-15


def test_reduction_and_atkin_lehner_characters_on_eisenstein_series():
    """F_eps(z) = chi F_eps(z*) for the four Eisenstein combinations at s = 2: a route that uses no Maass form."""
    mp.mp.dps = 20
    rng = np.random.default_rng(1)
    worst, low, moved = 0.0, 1.0, 0
    for x in rng.uniform(-0.5, 0.5, 12):
        y = float(rng.uniform(0.004, 0.03))
        xs, ys, par, _ = D.reduce_point(float(x), y, 143)
        low, moved = min(low, ys), moved + (ys > y * 1.0001)
        Ez, Es = S.eisenstein(mp.mpc(float(x), y), 2), S.eisenstein(mp.mpc(xs, ys), 2)
        for e11 in (1, -1):
            for e13 in (1, -1):
                Fz, Fs = (v["inf"] + e11 * e13 * v["0"] + e13 * v["1/11"] + e11 * v["1/13"] for v in (Ez, Es))
                chi = (e11 ** par[0]) * (e13 ** par[1])
                worst = max(worst, float(abs(Fz - chi * Fs) / max(abs(Fz), abs(Fs))))
    assert moved == 12 and low > 0.0098                          # every point was carried up, into the domain
    assert worst < 1e-8                                          # the lab's record: 2.4e-13


def _level1_word(x, y):
    """z carried into the domain of SL(2, Z) by T and S, with the sign of the word: chi(T) = chi(S) = -1 for nu_eta^12."""
    sign = 1
    for _ in range(100000):
        k = math.floor(x + 0.5)
        x -= k
        if k & 1:
            sign = -sign
        if x * x + y * y >= 1.0:
            return x, y, sign
        d = x * x + y * y
        x, y, sign = -x / d, y / d, -sign
    raise RuntimeError("no convergence")


def test_eta_twist_sign_is_the_character_of_eta_to_the_twelfth():
    """The sign carried through the reduction equals the sign of the motion as a word in S and T.

    nu_eta^12 is -1 on both generators of SL(2, Z), so it is (-1)^(length of the word).  The reduction of
    Gamma_0(143) (Atkin-Lehner steps off) moves z to z* by an element of SL(2, Z); carrying both to the domain
    of SL(2, Z) by S and T gives that element as a word, and its sign must be the one reduce_point returns.
    """
    mats = [(1, 0, 0, 1), (1, 1, 0, 1), (1, 0, 1, 1), (0, 1, 1, 0), (0, 1, 1, 1), (1, 1, 1, 0)]      # GL(2, F_2)
    mul = lambda A, B: (A[0] * B[0] + A[1] * B[2], A[0] * B[1] + A[1] * B[3], A[2] * B[0] + A[3] * B[2], A[2] * B[1] + A[3] * B[3])
    assert all(D.sign_mod2(*mul(A, B)) == D.sign_mod2(*A) * D.sign_mod2(*B) for A in mats for B in mats)     # a homomorphism
    assert sorted(D.sign_mod2(*A) for A in mats) == [-1, -1, -1, 1, 1, 1]
    assert D.sign_mod2(1, 1, 0, 1) == -1 and D.sign_mod2(0, -1, 1, 0) == -1                                # T and S
    with pytest.raises(ValueError):
        D.sign_mod2(2, 0, 0, 1)                                  # an even determinant has no sign
    rng = np.random.default_rng(7)
    seen = {1: 0, -1: 0}
    for x in rng.uniform(-0.5, 0.5, 60):
        y = float(rng.uniform(0.0008, 0.006))
        xs, ys, _, sign = D.reduce_point(float(x), y, 143, use_al=False)
        if ys <= 1.0001 * y:                                     # under no circle of Gamma_0(143): not moved
            continue
        a, b = _level1_word(float(x), y), _level1_word(xs, ys)
        assert abs(a[0] - b[0]) < 1e-7 and abs(a[1] - b[1]) < 1e-7          # one orbit of SL(2, Z)
        assert a[2] == sign * b[2]
        seen[sign] += 1
    assert min(seen.values()) >= 8 and sum(seen.values()) >= 30  # points were moved, with both signs: the check is not empty


def test_level_13_form_agrees_with_lmfdb():
    st = H.setup(13)
    r = H.refine(0.97081541, 1, (-1,), False, st, half_width=1e-6)
    assert abs(r["R"] - 0.97081541) < 5e-8                       # LMFDB 13.1 gives eight digits
    assert r["agreement"] < 1e-9 and r["hecke"] < 1e-9
    assert abs(r["coefficients"][12] * math.sqrt(13) - 1.0) < 1e-9          # c_13 sqrt 13 = -eps_13


def test_level_11_forms_of_the_record():
    st = H.setup(11)
    for R0, parity, eps in ((2.48359105955, 0, 1), (2.96820576382, 1, 1), (2.03309099386, 1, -1)):
        r = H.refine(R0, parity, (eps,), False, st, half_width=1e-6)
        assert abs(r["R"] - R0) < 1e-9 and r["agreement"] < 1e-9 and r["hecke"] < 1e-9
        assert abs(r["coefficients"][10] * math.sqrt(11) + eps) < 1e-8      # c_11 sqrt 11 = -eps_11


def test_level_143_three_eigenvalues_found_again_on_two_horocycles():
    st = H.setup()
    assert st.M == 632 and st.horocycles[0].Q == 788             # the records' setup
    assert min(h.ys.min() for h in st.horocycles) > 0.00985      # the reduced points are in the domain
    for R0, parity, eps, twist in ((0.3738518157, 1, (-1, -1), False), (0.4449803610, 0, (-1, -1), False),
                                   (0.4627415219, 1, (-1, -1), True)):
        r = H.refine(R0, parity, eps, twist, st)
        assert abs(r["R"] - R0) < 1e-9 and r["agreement"] < 1e-9 and r["hecke"] < 1e-9, (R0, r["R"], r["agreement"], r["hecke"])


def test_control_a_wrong_sector_or_parity_does_not_pass():
    R1 = 0.373851815749                                           # odd, (-, -)
    good = H.solve_pair(R1, 1, (-1, -1))
    assert H.agreement(*good) < 1e-9
    for parity, eps in ((1, (1, 1)), (1, (1, -1)), (0, (-1, -1))):
        assert H.agreement(*H.solve_pair(R1, parity, eps)) > 1e-3
    with pytest.raises(ValueError):
        H.refine(R1, 1, (1, 1))                                   # no sign change in the wrong sector
    off = H.solve_pair(R1 + 0.01, 1, (-1, -1))                    # the right sector, the wrong R
    assert H.agreement(*off) > 1e-3


def test_detector_dips_at_the_eigenvalue():
    hc, fr = H.detector()
    assert fr.size == 300 and hc.Q == 381
    chi = hc.chi((-1, -1))
    at = H.sigma_min(0.373851815749, hc, fr, 1, chi)
    near = [H.sigma_min(0.373851815749 + d, hc, fr, 1, chi) for d in (-0.01, 0.01)]
    assert at < 1e-9 and min(near) > 1e-6 and min(near) / at > 1e4
    wrong = H.sigma_min(0.373851815749, hc, fr, 1, hc.chi((1, 1)))
    assert wrong / at > 1e4


def test_error_paths():
    with pytest.raises(ValueError):
        D.reduce_point(0.1, -1.0)
    with pytest.raises(ValueError):
        D.reduce_point(0.1, 0.01, 12)                             # not squarefree
    with pytest.raises(ValueError):
        D.prime_factors(0)
    with pytest.raises(ValueError):
        D.Horocycle(143, 0.05, 20)                                # not below the domain
    hc = H.setup(11).horocycles[0]
    with pytest.raises(ValueError):
        hc.chi((1, 1))                                            # one prime, two signs
    with pytest.raises(ValueError):
        H.setup(2).horocycles[0].chi((1,), twist=True)             # the eta-twist needs an odd level
    with pytest.raises(ValueError):
        H.hecke_residual(np.ones(10))
    with pytest.raises(ValueError):
        H.detector(11)


@pytest.mark.slow
def test_level_1_first_eigenvalues_with_mpmath_bessel():
    st = H.setup(1, 12, (0.72, 0.62))
    for R0, parity, known in ((9.5336952614, 1, 9.5336952613535575), (13.7797513519, 0, 13.779751351890738)):
        r = H.refine(R0, parity, (), False, st, half_width=1e-6, upto=4)
        assert abs(r["R"] - known) < 1e-9 and r["agreement"] < 1e-9
        c = r["coefficients"]
        assert abs(c[1] * c[2] - c[5]) < 1e-8                     # c_2 c_3 = c_6
