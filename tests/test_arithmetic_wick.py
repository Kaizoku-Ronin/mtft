"""
tests/test_arithmetic_wick.py
=============================
First test suite for mtft.arithmetic_wick (new in v0.7.0): weight parity
with mtft.arithmetic, thermodynamic self-consistency of both ensembles,
the Wick rotation bridge, and the v0.7.0 fixes (canonical Y_C, the
underflow-branch TypeError).
"""
import math

import pytest

from mtft.arithmetic import weight
from mtft.arithmetic_wick import (
    Y_C,
    Y_S1,
    Y_S2,
    compute_skeleton_weights,
    compute_weights,
    dirichlet_ensemble,
    laplace_ensemble,
    sieve_primes,
    wick_at_critical_depths,
    wick_rotate,
)
from mtft.constants import CriticalDepths


# ── canonical constants (audit lock) ─────────────────────────────

def test_critical_depths_are_canonical():
    assert Y_C == CriticalDepths.y_conf == 0.18174
    assert Y_S1 == CriticalDepths.y_s1
    assert Y_S2 == CriticalDepths.y_s2


# ── weights ──────────────────────────────────────────────────────

def test_weights_match_arithmetic_module():
    w = compute_weights(300)
    for n in range(1, 301):
        assert w[n] == pytest.approx(weight(n), abs=1e-12), f"n={n}"


def test_skeleton_weights_supported_on_primes():
    n_max = 200
    ws = compute_skeleton_weights(n_max)
    primes = set(sieve_primes(n_max))
    for n in range(2, n_max + 1):
        if n in primes:
            assert ws[n] == pytest.approx(math.log(n) / n)
        else:
            assert ws[n] == 0.0


# ── Laplace ensemble ─────────────────────────────────────────────

def test_laplace_ensemble_thermodynamics():
    y, n_max = CriticalDepths.y_conf, 800
    ens = laplace_ensemble(y, n_max)
    w = compute_weights(n_max)
    two_pi_y = 2 * math.pi * y

    Z = sum(w[n] * math.exp(-two_pi_y * n) for n in range(2, n_max + 1))
    assert ens.partition_fn == pytest.approx(Z, rel=1e-12)
    assert ens.free_energy == pytest.approx(math.log(Z), rel=1e-12)

    mean = sum(n * w[n] * math.exp(-two_pi_y * n) for n in range(2, n_max + 1)) / Z
    assert ens.mean_energy == pytest.approx(mean, rel=1e-10)

    # C = (2πy)² Var and S = 2πy⟨n⟩ + log Z
    assert ens.specific_heat == pytest.approx(
        (2 * math.pi * y) ** 2 * ens.energy_variance, rel=1e-10
    )
    assert ens.entropy == pytest.approx(
        two_pi_y * ens.mean_energy + ens.free_energy, rel=1e-10
    )
    assert ens.energy_variance > 0
    assert ens.stiffness_N3 > 0


def test_laplace_underflow_branch_regression():
    # v0.7.0 fix: this used to raise TypeError (field passed as beta=).
    ens = laplace_ensemble(80.0)
    assert ens.y == 80.0
    assert ens.partition_fn == 0
    assert ens.free_energy == float("-inf")
    assert ens.mean_energy == 0


# ── Dirichlet ensemble ───────────────────────────────────────────

def test_dirichlet_partition_function():
    beta, n_max = 2.5, 400
    ens = dirichlet_ensemble(beta, n_max)
    w = compute_weights(n_max)
    Z = sum(w[n] * n**-beta for n in range(2, n_max + 1))
    assert ens.partition_fn == pytest.approx(Z, rel=1e-12)
    assert ens.mean_energy > 0  # ⟨log n⟩ > 0


def test_dirichlet_partition_decreases_with_beta():
    z2 = dirichlet_ensemble(2.0, 400).partition_fn
    z3 = dirichlet_ensemble(3.0, 400).partition_fn
    assert z3 < z2


# ── Wick rotation ────────────────────────────────────────────────

def test_wick_rotate_consistency():
    y = CriticalDepths.y_conf
    beta = 2 * math.pi * y
    wr = wick_rotate(y, beta, n_max=600)
    assert wr.y == y
    assert wr.beta == beta
    assert wr.curvature_ratio == pytest.approx(
        wr.dirichlet.energy_variance / wr.laplace.energy_variance, rel=1e-12
    )
    assert wr.entropy_difference == pytest.approx(
        wr.dirichlet.entropy - wr.laplace.entropy, rel=1e-10
    )
    assert wr.mean_energy_ratio == pytest.approx(
        wr.dirichlet.mean_energy / wr.laplace.mean_energy, rel=1e-10
    )


def test_wick_at_critical_depths_structure():
    res = wick_at_critical_depths(n_max=500)
    assert set(res) == {"y_s1", "y_c", "y_s2"}
    assert res["y_c"].y == pytest.approx(0.18174)
    assert res["y_s1"].y == pytest.approx(0.1236)
    assert res["y_s2"].y == pytest.approx(0.2106)
    for wr in res.values():
        assert wr.beta == pytest.approx(2 * math.pi * wr.y)
