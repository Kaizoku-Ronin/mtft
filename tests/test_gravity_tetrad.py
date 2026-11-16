"""mtft.gravity.tetrad: a tetrad for the squared hopfion, frozen into its null flow, and the two helicities.

Expected numbers are the lab's recorded ones (MP23, 14 November 2026: studies/gravity_2026nov/lab/mp23) or are computed in the
test by another route, and the test says which.  The lab took its derivatives by nested differences at 40 digits; here they
are exact (second-order jets), which is a routine of another kind."""
import json
import math
import random
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from mtft.gravity import fields as FD
from mtft.gravity import knotted as KN
from mtft.gravity import tetrad as TE

MP23 = json.loads((Path(__file__).resolve().parents[1] / "studies" / "gravity_2026nov" / "lab" / "mp23" / "mp23_numbers.json").read_text(encoding="utf-8"))


def _events(n, seed):
    rng = random.Random(seed)
    return [(rng.uniform(-1, 2), rng.gauss(0, 1.2), rng.gauss(0, 1.2), rng.gauss(0, 1.2)) for _ in range(n)]


def test_the_tetrad_is_frozen_and_its_curvature_is_the_squared_hopfion():
    """MP23's Q3, which is T-G1's criterion: 60 events with t from −1 to 2."""
    worst = {}
    for x in _events(60, 77):
        for key, v in TE.exhibit(x).items():
            worst[key] = max(worst.get(key, 0.0), v)
    assert set(worst) == {"frozen", "E", "B", "ricci", "Rk"} and max(worst.values()) < 1e-12, worst
    rec = MP23["Q3"]                                                                 # the lab's figures were of the same size
    assert max(rec["|i_k F| / |F|"], rec["E - Re W"], rec["B + Im W (B_ij = (1/2) eps_jkl R_0ikl)"], rec["Ricci"]) < 1e-12
    # a control, so that the zeros mean something: with the other sign of the magnetic part the difference is of the field's own size
    x = (0.4, 0.6, -0.3, 0.5)
    c = TE.curvature(x)
    scale = np.max(np.abs(c["W"]))
    assert np.max(np.abs(c["E"] - c["W"].real)) < 1e-13 * scale and np.max(np.abs(c["B"] + c["W"].imag)) < 1e-13 * scale
    assert np.max(np.abs(c["B"] - c["W"].imag)) > 0.1 * scale                        # MP23: with the other sign there is no solution
    v = FD.flow(x[0], np.array(x[1:]))                                               # MP23's Q0: W is trace-free and across the flow
    assert abs(np.trace(c["W"])) < 1e-13 * scale and np.max(np.abs(c["W"] @ v)) < 1e-13 * scale
    assert np.max(np.abs(c["W"] - FD.squared(x[0], np.array(x[1:])))) < 1e-13 * scale    # two codes of the squared hopfion
    assert np.max(np.abs(c["ricci"])) < 1e-13 * scale and np.max(np.abs(c["riemann"])) > scale / 10


def test_the_same_at_forty_digits_and_over_the_rationals():
    with mp.workdps(40):
        for x in _events(8, 5):
            got = TE.exhibit(tuple(mp.mpf(v) for v in x), "mp")
            assert max(got.values()) < mp.mpf("1e-30"), got
    rng = random.Random(9)
    for _ in range(10):                                                              # exact: every difference is the rational number 0
        x = (Fr(rng.randint(-8, 16), 8), Fr(rng.randint(-12, 12), 7), Fr(rng.randint(-12, 12), 5), Fr(rng.randint(-12, 12), 9))
        assert TE.exhibit(x, "exact")["exact"] is True
    a, b = TE.GQ(Fr(1, 2), Fr(-3, 4)), TE.GQ(2, 5)                                   # the exact numbers themselves
    assert a * b / b == a and (a - a) == 0 and 1 / (1 / a) == a and (a * a.conjugate()).imag == 0 and a + 1 - 1 == a


def test_the_certificate_for_every_event():
    """Exact zeros at ten random rational events, and what they prove: a false identity would pass with chance below 3e-21."""
    cert = TE.certificate()
    assert cert["holds"] is True and cert["points"] == 10 and cert["values_per_coordinate"] == 4001 and cert["differences"] == 16 + 9 + 9 + 16 + 64
    assert cert["degree_bound"] == 22 and cert["false_pass_below"] == pytest.approx(114 * (22 / 4001) ** 10) and cert["false_pass_below"] < 3e-21
    # a control in the same exact arithmetic: the other sign of the magnetic part is not zero
    x = (Fr(3, 10), Fr(-7, 20), Fr(11, 25), Fr(1, 2))
    _, _, d2B, s = TE._parts(x, "exact")
    R, W = TE._riemann(d2B), TE._fields(s)[1]
    Bm = sum(e * R[0][1][k + 1][m + 1] for (j, k, m), e in TE._EPS3.items() if j == 1) / 2
    assert Bm + W[0][1].imag == 0 and Bm - W[0][1].imag != 0


def test_the_field_strength_is_the_light_two_form_times_a_weight():
    for x in _events(20, 3):
        F = TE.field_strength(x)
        assert np.max(np.abs(F + np.transpose(F, (2, 1, 0)))) < 1e-15 * (1 + np.max(np.abs(F)))   # a 2-form for each leg
        assert np.max(np.abs(F - TE.strength_from_weight(x))) < 1e-12 * np.max(np.abs(F))         # 𝓕_μ = Re[G_μ dα ∧ dβ]
        al, be = FD.labels(x[0], np.array(x[1:]))
        da, db = np.array([0, 1, -1j, 0]), np.array([1, 0, 0, -1])
        assert TE.weight(x) == pytest.approx(-0.5 * (1 - al) * da + 0.5j * be * db, abs=1e-12)    # linear in the labels (MP23's fitted form)
        k = TE.null_flow(x)
        assert k[0] == 1 and k[1:] == pytest.approx(FD.flow(x[0], np.array(x[1:])), abs=1e-13) and -k[0] ** 2 + k[1:] @ k[1:] == pytest.approx(0, abs=1e-13)
        h = TE.metric(x)
        assert np.max(np.abs(h - h.T)) == 0 and np.max(np.abs(TE.potential(x) - TE.potential(x).T)) > 1e-6   # the frame's turn is not zero


def test_the_two_helicities():
    """MP23's Q4, which is T-G2: (2/3)π² and π², on the label sphere and on two slices of space."""
    rec = MP23["Q4"]
    sphere = TE.helicities_on_sphere()
    assert sphere["H_em"] == pytest.approx(math.pi ** 2, rel=1e-12) and sphere["H_g"] == pytest.approx(2 * math.pi ** 2 / 3, rel=1e-12)
    assert sphere["ratio"] == pytest.approx(2 / 3, abs=1e-12) and TE.H_G / TE.H_EM == pytest.approx(2 / 3)
    # the lab's file holds the same two numbers with a minus sign: its order of the sphere's coordinates, (η, ξ₁, ξ₂), is the
    # orientation opposite to space's, and its test took absolute values.  Here the sphere is oriented as space orients it.
    assert sphere["H_g"] == pytest.approx(-rec["on the label sphere"]["H_g"], rel=1e-13) and sphere["H_em"] == pytest.approx(-rec["on the label sphere"]["H_em"], rel=1e-13)
    for t, key in ((0.0, "in space, t = 0"), (0.8, "in space, t = 0.8")):
        got = TE.helicities_on_slice(t)
        assert got["H_g"] == pytest.approx(rec[key]["H_g"], rel=1e-12) and got["H_em"] == pytest.approx(rec[key]["H_em"], rel=1e-12)
        assert got["H_g"] == pytest.approx(TE.H_G, rel=1e-8) and got["H_em"] == pytest.approx(TE.H_EM, rel=1e-8)   # the same in time, and the sign of the sphere's
        assert got["ratio"] == pytest.approx(2 / 3, abs=1e-8)
    assert rec["3 pi^2 / 4"] == pytest.approx(0.75 * math.pi ** 2) and abs(sphere["H_g"] - rec["3 pi^2 / 4"]) > 0.8   # MP23's prediction of 3/4: a recorded miss
    assert KN.helicity(1, 1) == pytest.approx(TE.H_EM) and KN.helicity_on_sphere(1, 1) == pytest.approx(sphere["H_em"], rel=1e-12)   # two codes, one number


def test_the_orientation_of_the_label_sphere():
    """Found in this build: (ξ₁, η, ξ₂) is the orientation that space gives the label sphere."""
    rng = np.random.default_rng(3)
    for X in rng.normal(0, 1, (20, 3)):
        def hopf(Y):
            al, be = FD.labels(0.3, Y)
            return np.array([np.angle(al), math.acos(min(1.0, abs(al))), np.angle(be)])
        J = np.array([(hopf(X + 1e-6 * e) - hopf(X - 1e-6 * e)) / 2e-6 for e in np.eye(3)]).T
        assert np.linalg.det(J) > 0
