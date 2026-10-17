"""CC-69 (2026-10-17): the repaired lattice action — the arithmetic tower on Sym^n characters on plaquettes — and what follows from it.
Every number in docs/SM/CC69_REPAIRED_ACTION.md is tied to the code here."""
import math

import numpy as np
import pytest
from scipy.special import iv

from mtft.legend import REGISTRY
from mtft.matter import puzzle as Z
from mtft.research import repaired_action as RA


@pytest.mark.parametrize("N", [2, 3])
def test_weyl_integration_is_orthonormal_on_characters(N):
    G = RA.integrator(N); reps = G.reps(5)
    gram = np.array([[np.real(G.integrate(G.chi(a) * np.conj(G.chi(b)))) for b in reps] for a in reps])
    assert np.abs(gram - np.eye(len(reps))).max() < 1e-6


@pytest.mark.parametrize("N", [2, 3])
def test_sym_recursion_dimensions(N):
    h = [1.0, float(N)]
    for n in range(2, 9):
        h.append(N * h[n - 1] - h[n - 2] if N == 2 else N * h[n - 1] - N * h[n - 2] + (h[n - 3] if n >= 3 else 0))
    assert all(abs(h[n] - RA.sym_dim(N, n)) < 1e-9 for n in range(1, 9))


def test_su2_wilson_coefficients_are_bessel():
    c = RA.coefficients(2, 2.0, 0.0, top=3)
    assert abs(c[1] / c[0] - iv(2, 2.0) / iv(1, 2.0)) < 1e-8 and abs(c[2] / c[0] - iv(3, 2.0) / iv(1, 2.0)) < 1e-8


@pytest.mark.parametrize("N", [2, 3])
def test_every_character_coefficient_is_nonnegative(N):            # positive-type weight: the Osterwalder–Seiler hypothesis
    for beta in (0.0, 2.5, 6.0):
        for kappa in (0.0, 1.0, 100.0, 1000.0):
            assert min(RA.coefficients(N, beta, kappa).values()) > -1e-12


def test_tower_is_a_per_cent_correction_at_kappa_one():
    for N, beta, lo, hi in ((2, 1.0, 0.015, 0.025), (3, 5.7, 0.005, 0.015)):
        shift = RA.u_fundamental(N, beta, 1.0) / RA.u_fundamental(N, beta, 0.0) - 1
        assert lo < shift < hi
    assert abs(RA.beta_equivalent(3, 5.7, 1.0) - 5.7577) < 2e-3 and abs(RA.beta_equivalent(2, 2.3, 1.0) - 2.3509) < 2e-3


def test_tower_alone_is_second_order_in_kappa():
    u = {k: RA.u_fundamental(2, 0.0, k) for k in (0.5, 1.0, 2.0)}
    assert abs(u[1.0] / u[0.5] - 4) < 0.25 and abs(u[2.0] / u[1.0] - 4) < 0.25 and u[1.0] < 1e-4


def test_strong_coupling_observables():
    s = RA.strong_coupling(3, 5.7, 1.0)
    assert abs(s["u_f"] - 0.40658) < 1e-4 and abs(s["sigma_a2"] + math.log(s["u_f"])) < 1e-12 and abs(s["glueball_a"] - 4 * s["sigma_a2"]) < 1e-12


def test_centre_symmetry_and_nality():
    for N in (2, 3):
        assert RA.centre_cost(N, "power") > 1e-3 and RA.centre_cost(N, "sym") > RA.centre_cost(N, "power")
    assert RA.nality(3, (1, 0)) == 1 and RA.nality(3, (1, 1)) == 2 and RA.nality(3, (2, 1)) == 0 and RA.nality(3, (2, 0)) == 2 and RA.nality(2, 3) == 1
    cls = RA.tower_nality_classes(3)
    assert set(cls) == {0, 1, 2} and all(v > 0 for v in cls.values())
    ov = [RA.vortex_overlap(3, b, 1.0) for b in (0.5, 2.0, 5.7)]
    assert ov[0] > 0.95 and ov[0] > ov[1] > ov[2] > 0


def test_legend_and_map():
    e = REGISTRY["cc69_repaired_action"]
    assert "Osterwalder" in e.nature and "cc35_reflection_positivity" in e.upstream and e.ref == "docs/SM/CC69_REPAIRED_ACTION.md"
    sub = Z.subs()["I.5"][1]
    assert sub.status == "answered" and "CC-69" in sub.answer and "CC69_REPAIRED_ACTION.md" in sub.docs and "continuum limit" in sub.answer
