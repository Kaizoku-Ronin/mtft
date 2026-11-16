"""mtft.gravity.knotted and mtft.gravity.fields: the tube round a core line, its rotation number and its shear.

Expected numbers are the lab's recorded ones (MP09 and MP21, 10 and 14 November 2026: studies/gravity_2026nov/lab) or are
computed in the test by another route, and the test says which."""
import json
import math
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import pytest

import mtft.gravity as GR
from mtft.gravity import fields as FD
from mtft.gravity import knotted as KN
from mtft.plane import kernel as PK

LAB = Path(__file__).resolve().parents[1] / "studies" / "gravity_2026nov" / "lab"
MP09 = json.loads((LAB / "mp09" / "mp81_numbers.json").read_text(encoding="utf-8"))
MP21 = json.loads((LAB / "mp21" / "mp21_numbers.json").read_text(encoding="utf-8"))


def _pair(key):
    return tuple(int(v) for v in key.split(","))


def test_the_front_page_says_what_this_is_not():
    doc = GR.__doc__
    assert "weak wave on flat" in doc and "Newton's constant is not derived" in doc and "is not the same as new" in doc
    assert set(GR.OPEN) == {"T-G3", "H-G7"} and set(GR.NOT_DERIVED) == {"N-G1"}
    for name in ("Kedia", "Swearngin", "Garnier", "Huang", "Asenjo"):                # whose work is whose
        assert name in doc


def test_the_rotation_number_and_the_pairs_where_it_is_whole():
    assert KN.whole_pairs(60) == [(1, 1), (1, 8), (1, 49), (32, 49)]                # MP09's P6, for coprime pairs
    assert KN.pell(4) == [1, 8, 49, 288]                                             # MP09's P6: Pell's equation for p = 1
    every = KN.whole_pairs(60, coprime=False)                                        # found here: without coprimality there are more
    assert set(KN.whole_pairs(60)) < set(every) and {(4, 4), (9, 9), (18, 32)} <= set(every)
    assert [p for p, q in every if p == q] == [k * k for k in range(1, 8)]           # (p, p) exactly when p is a square: 4p³
    for p, q in every:
        assert math.isqrt(KN.nu_squared(p, q)) ** 2 == 2 * p * q * (p + q)
    for key, v in MP09["nu"].items():
        assert KN.nu(*_pair(key)) == pytest.approx(v, rel=1e-15)
    assert KN.seifert_winding(1, 8) == 4.0                                           # MP09: "ν = 12 and winding 4 against Seifert"


def test_the_matrix_of_one_circuit():
    for key, rec in MP09["M"].items():                                               # the lab's closed form, entry by entry
        M = KN.tube_matrix(*_pair(key))
        assert [M[0][0], M[0][1], M[1][0], M[1][1]] == pytest.approx(rec, abs=2e-14)
    for p, q in ((2, 3), (3, 2), (2, 5), (3, 5)):
        L, M = KN.generator(p, q), KN.tube_matrix(p, q)
        E = PK.expm(L, 2 * math.pi * (p + q))                                        # a second route: the plane's own exponential
        assert max(abs(E[i][j] - M[i][j]) for i in range(2) for j in range(2)) < 1e-12
        assert L[0][0] * L[1][1] - L[0][1] * L[1][0] == pytest.approx(2 * p * q / (p + q))
        assert M[0][0] * M[1][1] - M[0][1] * M[1][0] == pytest.approx(1.0) and PK.read(M)["kind"] == "ell"
    # MP09's result: the traces by the formula, -0.0506786 and 0.9871991 (its plan had printed 0.9868 for the second, by hand)
    assert KN.tube_matrix(2, 3)[0][0] * 2 == pytest.approx(-0.0506786, abs=5e-8)
    assert KN.tube_matrix(2, 5)[0][0] * 2 == pytest.approx(0.9871991, abs=5e-8)
    assert KN.tube_matrix(1, 8) == ((1.0, 0.0), (0.0, 1.0)) and KN.tube_matrix(1, 1) == ((1.0, 0.0), (0.0, 1.0))


def test_the_shear_table_of_mp21():
    for key, rec in MP21["table"].items():
        p, q = _pair(key)
        assert KN.nu(p, q) == pytest.approx(rec["nu"], rel=1e-14) and KN.shear(p, q) == pytest.approx(rec["d iota/du"], rel=1e-13, abs=1e-15)
        assert KN.h_max(p, q) == pytest.approx(rec["h_max"], rel=1e-14) and KN.shear_dh(p, q) == pytest.approx(rec["d iota/dh"], rel=1e-13, abs=1e-15)
    assert KN.shear(1, 8) == pytest.approx(35 / 36) and 12 * KN.shear_factor(1, 8) == Fr(35, 36) and KN.shear(1, 1) == 0.0   # MP21: 35/36 for (1, 8)
    assert [round(KN.shear(p, q), 6) for p, q in ((1, 1), (1, 2), (2, 3), (2, 5))] == [0.0, 0.192450, 1.118862, 1.859339]   # MP21's table, as printed
    for q in range(1, 40):                                                           # MP21: for p = 1 the bracket is (q + 2)(q − 1)
        assert KN.shear_factor(1, q) == Fr((q + 2) * (q - 1), 12 * q * (q + 1))
    assert all(KN.shear_factor(p, q) > 0 for p in range(1, 30) for q in range(1, 30) if (p, q) != (1, 1))   # zero for the hopfion only
    for p, q in ((2, 3), (1, 8)):                                                    # the three coefficients are the derivatives of λ at the core
        a, b, c = KN.area_coefficients(p, q)
        assert KN.nu(p, q) == pytest.approx(2 * q * math.sqrt(a))
        K = Fr(15, 32) * b ** 2 / a ** 3 - Fr(3, 8) * c / a ** 2 - Fr(1, 16)
        assert -4 * K * 2 * q * math.sqrt(a) == pytest.approx(KN.shear(p, q), rel=1e-13)


def test_the_derivation_is_exact():
    checks = KN.exact_checks()
    assert len(checks) == 7 and all(checks.values()), checks


def test_the_shear_from_the_profile_integral():
    """MP21's P1, the second route: the integral knows nothing of the closed form."""
    worst = 0.0
    for key, rec in MP21["P1"]["pairs"].items():
        p, q = _pair(key)
        s = KN.shear_from_profile(p, q)
        assert s == pytest.approx(rec["extrapolated"], abs=1e-10)                    # the lab's own extrapolated value
        worst = max(worst, abs(s - KN.shear(p, q)) / KN.nu(p, q))
    assert worst < 1e-8 and worst == pytest.approx(MP21["P1"]["worst difference / nu"], rel=1e-2)
    assert abs(KN.iota(1, 1, "1e-3") - 2) < 1e-20                                    # the hopfion: 2 on every torus
    assert abs(KN.iota(1, 1, 0.5) - 2) < 1e-20
    # toward the edge the winding tends to 2pq (MP09's P5 recorded 11.9994 and 19.9798 at h = 1e-6 h_max)
    assert float(KN.iota(2, 3, 1 - 1e-6)) == pytest.approx(11.9994, abs=5e-5) and float(KN.iota(2, 5, 1 - 1e-6)) == pytest.approx(19.9798, abs=5e-5)
    assert KN.iota_edge(2, 3) == 12 and KN.iota_edge(2, 5) == 20
    levels = [float(KN.iota(2, 3, u)) for u in (0.01, 0.1, 0.5, 0.9, 0.999)]
    assert levels == sorted(levels) and KN.nu(2, 3) < levels[0] and levels[-1] < 12  # it rises from ν to 2pq without a turn


def test_the_shear_from_the_section_as_a_differential_equation():
    """The third route, the package's own: the period of the cross-section's Hamiltonian system, integrated."""
    for p, q in ((2, 3), (1, 8), (2, 5)):
        S, v = KN.shear(p, q), KN.nu(p, q)
        for u in (0.004, 0.016):
            ode, quad = KN.iota_ode(p, q, u), float(KN.iota(p, q, u))
            assert ode == pytest.approx(quad, rel=1e-11)
            assert abs(ode - v - S * u) < 3 * S * u * u                              # MP21's F1: the next term is of ordinary size
    assert KN.iota_ode(1, 1, 0.3) == pytest.approx(2.0, abs=1e-11)


def test_the_fields_in_space():
    rng = np.random.default_rng(5)
    for t in (0.0, 0.8):
        X = rng.normal(0, 1.2, (200, 3))
        al, be = FD.labels(t, X)
        ga, gb = FD.grad_labels(t, X)
        assert np.max(np.abs(np.abs(al) ** 2 + np.abs(be) ** 2 - 1)) < 1e-13         # the labels lie on the sphere
        F11 = FD.F(t, X)
        assert np.max(np.abs(np.cross(ga, gb) - F11)) < 1e-13 * np.max(np.abs(F11))  # ∇α × ∇β = n/d³
        assert np.max(np.abs(np.sum(F11 * F11, -1))) < 1e-12 * np.max(np.abs(F11)) ** 2   # a null field: E·B = 0 and |E| = |B|
        assert np.max(np.abs(FD.chart(t, al, be) - X)) < 1e-11                       # the chart undoes the labels
        v = FD.flow(t, X)
        assert np.max(np.abs(np.linalg.norm(v, axis=1) - 1)) < 1e-13
        dt = 1e-6                                                                    # the labels ride the rays: (∂_t + v·∇)α = 0
        assert np.max(np.abs((FD.labels(t + dt, X)[0] - FD.labels(t - dt, X)[0]) / (2 * dt) + np.sum(v * ga, -1))) < 1e-7
        eta, xi1, xi2 = np.arccos(np.minimum(1, np.abs(al))), np.angle(al), np.angle(be)
        for p, q in ((1, 1), (2, 3), (1, 8)):
            B = FD.F(t, X, p, q).imag
            gH = np.real((p * al ** (p - 1) * be ** q)[:, None] * ga + (q * al ** p * be ** (q - 1))[:, None] * gb)
            assert np.max(np.abs(np.sum(B * gH, -1)) / (np.linalg.norm(B, axis=1) * np.linalg.norm(gH, axis=1))) < 1e-12   # Re(α^p β^q) is kept
            assert KN.H(al, be, p, q) == pytest.approx(np.cos(eta) ** p * np.sin(eta) ** q * np.cos(p * xi1 + q * xi2), abs=1e-13)
            # MP09's line field on the label sphere: η′ = −sin Θ sin η cos η, ξ₁′ = cos Θ sin²η, ξ₂′ = −cos Θ cos²η, up to a factor of one
            # sign.  Found here: with B = Im F_pq the factor is negative (MP09 wrote "positive"); no result depends on the sense.
            Th = p * xi1 + q * xi2
            de, d1, d2 = -np.sin(Th) * np.sin(eta) * np.cos(eta), np.cos(Th) * np.sin(eta) ** 2, -np.cos(Th) * np.cos(eta) ** 2
            want = np.stack([(-np.sin(eta) * de + 1j * np.cos(eta) * d1) * np.exp(1j * xi1), (np.cos(eta) * de + 1j * np.sin(eta) * d2) * np.exp(1j * xi2)], -1)
            got = np.stack([np.sum(ga * B, -1), np.sum(gb * B, -1)], -1)             # the velocity of the labels along B
            w4, g4 = np.concatenate([want.real, want.imag], -1), np.concatenate([got.real, got.imag], -1)
            cos = np.sum(w4 * g4, -1) / (np.linalg.norm(w4, axis=1) * np.linalg.norm(g4, axis=1))
            assert np.max(cos) < -1 + 1e-12


def test_helicity_and_energy():
    """Out of sample: neither closed form is in the lab's records.  Their ratio 1/(p + q) is Kedia et al.'s (2013)."""
    for p, q in ((1, 1), (1, 2), (2, 3), (2, 5), (1, 8)):
        closed = 2 * math.pi ** 2 * p * q * math.factorial(p) * math.factorial(q) / math.factorial(p + q)
        assert KN.helicity(p, q) == pytest.approx(closed) and KN.energy(p, q) == pytest.approx((p + q) * closed)
        assert KN.helicity_on_sphere(p, q) == pytest.approx(closed, rel=1e-12)       # the quadrature knows nothing of the closed form
    for p, q, t in ((1, 1, 0.0), (2, 3, 0.0), (2, 3, 0.8), (1, 8, 0.0)):
        assert KN.energy_in_space(p, q, t, n=72) == pytest.approx(KN.energy(p, q), rel=1e-9)   # and it does not change in time
    assert KN.helicity(1, 1) == pytest.approx(math.pi ** 2) and KN.energy(1, 1) == pytest.approx(2 * math.pi ** 2)
    # found here, against my first expectation: the density does not depend on the angles, so three points in each are enough
    assert KN.helicity_on_sphere(1, 8, n_xi=3) == pytest.approx(KN.helicity(1, 8), rel=1e-12)


def test_the_integrator_on_a_known_orbit():
    """The package's own integrator against an exact solution: the Kepler ellipse of eccentricity 0.6, one period."""
    def f(y):
        r3 = (y[0] * y[0] + y[1] * y[1]) ** 1.5
        return np.array([y[2], y[3], -y[0] / r3, -y[1] / r3])
    e = 0.6
    y0 = [1 - e, 0.0, 0.0, math.sqrt((1 + e) / (1 - e))]
    run = FD.integrate(f, y0, 2 * math.pi, 2 * math.pi / 100)
    assert np.max(np.abs(run["y"][-1] - np.array(y0))) < 1e-11 and run["error"] < 1e-9
    again = FD.integrate(f, y0, 4 * math.pi, 2 * math.pi / 200, event=lambda y: y[1], max_events=1)   # the next upward crossing of y = 0
    assert len(again["events"]) == 1 and again["events"][0][0] == pytest.approx(2 * math.pi, abs=1e-11)   # the start is not counted
