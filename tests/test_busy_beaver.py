"""
tests/test_busy_beaver.py
=========================
First test suite for mtft.busy_beaver (new in v0.7.0).  Locks the
audit-canonical data sources (HECKE_TRACES prefix = TRACE_TOTALS_50,
ORBIT_TRACES = the verified v0.6.1 table) and pins the exact BB values
that the module itself computes quickly and deterministically.
"""
import pytest

from mtft.busy_beaver import (
    GENUS,
    HECKE_TRACES,
    INDEX,
    LEVEL,
    ORBIT_TRACES,
    HeckeSign,
    bb_fatou,
    bb_genus,
    bb_sample,
    count_atms,
    dominant_sector,
    enumerate_atms,
    faulhaber_decompose,
    hecke_constraint_density,
    hecke_sign,
    verify_telescoping,
)
from mtft.x0_143 import ORBIT_TRACES_VERIFIED, TRACE_TOTALS_50
from mtft import x0_143


# ── canonical data sources ───────────────────────────────────────

def test_structural_constants_match_x0_143():
    assert LEVEL == x0_143.LEVEL == 143
    assert GENUS == x0_143.GENUS == 13
    assert INDEX == x0_143.INDEX == 168


def test_hecke_traces_prefix_is_verified_totals():
    assert len(HECKE_TRACES) == 200
    assert list(HECKE_TRACES[:50]) == list(TRACE_TOTALS_50)


def test_orbit_traces_is_the_verified_table():
    # v0.7.0 audit fix: the module re-exports the verified table itself.
    assert ORBIT_TRACES is ORBIT_TRACES_VERIFIED


def test_orbit_traces_columns_sum_to_totals():
    # Internal consistency that the pre-audit table failed at 7 of 9 primes.
    for p, (t1, t2, t3) in ORBIT_TRACES.items():
        assert t1 + t2 + t3 == HECKE_TRACES[p - 1], f"p={p}"


# ── Hecke sign oracle ────────────────────────────────────────────

def test_hecke_sign_values():
    assert hecke_sign(1) == HeckeSign.BOSONIC      # trace 11 > 0
    assert hecke_sign(5) == HeckeSign.FREE          # trace 0
    assert hecke_sign(6) == HeckeSign.FERMIONIC     # trace -4 < 0


def test_hecke_sign_bounds():
    with pytest.raises(ValueError):
        hecke_sign(0)
    with pytest.raises(ValueError):
        hecke_sign(len(HECKE_TRACES) + 1)


def test_hecke_constraint_density():
    d = hecke_constraint_density()
    assert d["fermionic"] == pytest.approx(0.485)
    assert d["bosonic"] == pytest.approx(0.465)
    assert d["free"] == pytest.approx(0.05)


def test_dominant_sector_from_verified_traces():
    # Follows directly from ORBIT_TRACES_VERIFIED (argmax of |trace|).
    assert dominant_sector(2) == 1    # (0, 3, 0)   → muon
    assert dominant_sector(3) == 2    # (-1, 0, 3)  → tau
    assert dominant_sector(17) == 1   # (-4, 6, 0)  → muon
    assert dominant_sector(19) == 2   # (2, 8, -10) → tau
    assert dominant_sector(23) == 2   # (7, -4, 11) → tau
    assert dominant_sector(29) is None  # not in the verified table


# ── enumeration ──────────────────────────────────────────────────

def test_count_and_enumerate_atms():
    assert count_atms(2) == (1296, 1296)
    machines = enumerate_atms(2)
    assert len(machines) == 1296


def test_enumeration_size_guard():
    with pytest.raises(ValueError):
        enumerate_atms(5)


# ── exact BB pins (fast, deterministic) ──────────────────────────

def test_bb_hecke_1_and_2():
    r1 = bb_genus(1)
    assert (r1.bb_value, r1.total_machines, r1.halting_count) == (0, 16, 16)
    r2 = bb_genus(2)
    assert (r2.bb_value, r2.total_machines, r2.halting_count, r2.cycle_count) == (
        0, 1296, 1224, 72,
    )


def test_bb_unconstrained_1_and_2():
    u1 = bb_genus(1, hecke_constrained=False)
    assert (u1.bb_value, u1.total_machines) == (7, 64)
    u2 = bb_genus(2, hecke_constrained=False)
    assert (u2.bb_value, u2.total_machines, u2.halting_count) == (13, 20736, 20306)


def test_bb_fatou_agrees_with_plain_run():
    f2 = bb_fatou(2)
    assert f2.bb_value == 0
    assert f2.halting_count == 1224


def test_bb_sample_deterministic_seed_143():
    s1 = bb_sample(4, sample_size=2000)
    s2 = bb_sample(4, sample_size=2000)
    assert s1.bb_value == s2.bb_value == 9
    assert s1.champion_steps == s2.champion_steps == 11
    assert s1.halting_count == s2.halting_count


# ── Faulhaber structure ──────────────────────────────────────────

def test_faulhaber_decomposition():
    res = {1: bb_genus(1), 2: bb_genus(2)}
    fd = faulhaber_decompose(res)
    assert fd.n_values == [1, 2]
    assert fd.bb_values == [0, 0]
    assert fd.naive_dominant == 8192          # 2^13, genus tape window
    assert fd.corrections == [8192, 8192]
    assert fd.degree_bound == 24              # canonical degree
    assert verify_telescoping(fd) is True
