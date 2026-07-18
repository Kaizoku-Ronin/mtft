"""
tests/test_jacobian.py
======================
Regression suite for the Jacobian stiffness engine (Paper 30, ported to
verified X₀(143) data in v0.7.0).  Eigenvalue pins were computed from the
audit-verified ORBIT_TRACE_F1/F2/F3 tables; the pre-audit Paper 30 values
(1.52628, 6.42705, 45.46824 at y = 0.1812) are documented in the module
docstring for provenance.
"""
import math

import numpy as np
import pytest

from mtft.arithmetic import weight
from mtft.constants import CriticalDepths, DELTA_X_MEASURED, FEIGENBAUM_DELTA
from mtft.jacobian import N_VERIFIED, SECTORS, JacobianStiffness
from mtft.x0_143 import ORBIT_TRACE_F1, ORBIT_TRACE_F2, ORBIT_TRACE_F3


@pytest.fixture(scope="module")
def eng():
    return JacobianStiffness()


# ── data provenance ──────────────────────────────────────────────

def test_traces_are_the_verified_tables(eng):
    assert N_VERIFIED == 50
    assert np.array_equal(eng._traces[0], np.array(ORBIT_TRACE_F1, dtype=float))
    assert np.array_equal(eng._traces[1], np.array(ORBIT_TRACE_F2, dtype=float))
    assert np.array_equal(eng._traces[2], np.array(ORBIT_TRACE_F3, dtype=float))


def test_weights_come_from_arithmetic_module(eng):
    for n in (1, 2, 3, 6, 12, 30, 143):
        assert eng._w_full[n - 1] == pytest.approx(weight(n), abs=1e-14)


def test_skeleton_plus_lambertization_is_full(eng):
    assert np.allclose(eng._w_skel + eng._w_lamb, eng._w_full, atol=1e-15)
    # skeleton is supported exactly on the primes
    assert eng._w_skel[2 - 1] == pytest.approx(math.log(2) / 2)
    assert eng._w_skel[4 - 1] == 0.0
    assert eng._w_skel[143 - 1] == 0.0  # 143 = 11 × 13 is composite


# ── eigenvalue pins ──────────────────────────────────────────────

def test_eigenvalues_at_paper30_depth(eng):
    _, evals, _ = eng.jacobian_matrix(0.1812)
    assert evals == pytest.approx([1.531802, 6.591154, 45.466583], rel=1e-5)


def test_eigenvalues_at_canonical_depth(eng):
    _, evals, _ = eng.jacobian_matrix(CriticalDepths.y_conf)
    assert evals == pytest.approx([1.515739, 6.507311, 45.102141], rel=1e-5)


def test_matrix_symmetric_positive_definite(eng):
    M, evals, evecs = eng.jacobian_matrix(CriticalDepths.y_conf)
    assert np.allclose(M, M.T)
    assert np.all(evals > 0)
    # eigh reconstruction
    assert np.allclose(evecs @ np.diag(evals) @ evecs.T, M, atol=1e-10)


def test_feigenbaum_ratios(eng):
    f = eng.feigenbaum_ratios(0.1812)
    assert f["lambda_3/lambda_2"] == pytest.approx(6.89812, rel=1e-4)
    assert f["lambda_2/lambda_1"] == pytest.approx(4.30288, rel=1e-4)
    assert f["delta_x_measured"] == DELTA_X_MEASURED
    assert f["delta_F"] == FEIGENBAUM_DELTA
    # documented post-audit deviations: ~4.8% and ~7.9%
    assert f["discrepancy_x"] == pytest.approx(0.0475, abs=0.002)
    assert f["discrepancy_F"] == pytest.approx(0.0785, abs=0.002)


def test_couplings(eng):
    c = eng.couplings(CriticalDepths.y_conf)
    assert set(c) == {"e-mu", "e-tau", "mu-tau"}
    assert c["e-mu"] == pytest.approx(0.0222, abs=5e-4)
    assert c["e-tau"] == pytest.approx(0.1881, abs=5e-4)
    assert c["mu-tau"] == pytest.approx(0.7172, abs=5e-4)


# ── projections and budget ───────────────────────────────────────

def test_sector_projection(eng):
    c = eng.sector_projection(CriticalDepths.y_conf)
    assert c == pytest.approx([-0.3993826, 1.272535, 1.398059], rel=1e-5)


def test_skeleton_budget(eng):
    b = eng.skeleton_budget(CriticalDepths.y_conf)
    assert tuple(b) == SECTORS
    sf = {k: v[0] for k, v in b.items()}
    assert sf["electron"] == pytest.approx(0.1351, abs=1e-3)
    assert sf["muon"] == pytest.approx(0.5324, abs=1e-3)
    assert sf["tau"] == pytest.approx(0.0478, abs=1e-3)
    for skel, lamb in b.values():
        assert skel + lamb == pytest.approx(1.0)


def test_lambertization_angle(eng):
    assert eng.lambertization_angle(CriticalDepths.y_conf) == pytest.approx(
        42.222, abs=0.01
    )


# ── stiffness sum ────────────────────────────────────────────────

def test_stiffness_value_and_kinds(eng):
    y = CriticalDepths.y_conf
    assert eng.stiffness(y) == pytest.approx(0.443043137, rel=1e-8)
    full = eng.stiffness(y, weight="full")
    skel = eng.stiffness(y, weight="skeleton")
    lamb = eng.stiffness(y, weight="lambertization")
    assert skel + lamb == pytest.approx(full, rel=1e-12)
    with pytest.raises(ValueError):
        eng.stiffness(y, weight="bulk")


# ── guards ───────────────────────────────────────────────────────

def test_shallow_depth_warns(eng):
    with pytest.warns(RuntimeWarning):
        eng.jacobian_matrix(0.05)


def test_n_max_floor():
    with pytest.raises(ValueError):
        JacobianStiffness(n_max=10)
