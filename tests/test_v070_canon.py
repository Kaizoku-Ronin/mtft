"""
tests/test_v070_canon.py
========================
Cross-module audit locks for the v0.7.0 additions.  Every module that
mentions X₀(143) data or the critical depths must source it from the
verified single-source-of-truth tables (mtft.x0_143, mtft.constants) —
these tests fail if a local pre-audit copy ever reappears.
"""
import numpy as np

from mtft import arithmetic_machine, arithmetic_wick, busy_beaver, music, x0_143
from mtft.constants import CriticalDepths
from mtft.jacobian import JacobianStiffness
from mtft.x0_143 import (
    ORBIT_TRACE_F1,
    ORBIT_TRACE_F2,
    ORBIT_TRACE_F3,
    ORBIT_TRACES_VERIFIED,
    TRACE_TOTALS_50,
)


def test_orbit_traces_single_source_of_truth():
    # identity, not equality: the modules must re-export the verified dict
    assert arithmetic_machine.ORBIT_TRACES is ORBIT_TRACES_VERIFIED
    assert busy_beaver.ORBIT_TRACES is ORBIT_TRACES_VERIFIED


def test_hecke_traces_agree_across_modules():
    assert list(arithmetic_machine.HECKE_TRACES) == list(busy_beaver.HECKE_TRACES)
    assert list(busy_beaver.HECKE_TRACES[:50]) == list(TRACE_TOTALS_50)
    assert len(busy_beaver.HECKE_TRACES) == 200


def test_orbit_totals_consistency():
    for p, (t1, t2, t3) in ORBIT_TRACES_VERIFIED.items():
        assert t1 + t2 + t3 == TRACE_TOTALS_50[p - 1], f"p={p}"


def test_structural_constants_agree():
    for mod in (arithmetic_machine, busy_beaver):
        assert mod.LEVEL == x0_143.LEVEL == 143
        assert mod.GENUS == x0_143.GENUS == 13
        assert mod.INDEX == x0_143.INDEX == 168


def test_critical_depths_canonical_everywhere():
    assert arithmetic_wick.Y_C == CriticalDepths.y_conf == 0.18174
    assert arithmetic_wick.Y_S1 == CriticalDepths.y_s1
    assert arithmetic_wick.Y_S2 == CriticalDepths.y_s2
    assert music.Y_C == CriticalDepths.y_conf


def test_jacobian_engine_uses_verified_tables():
    eng = JacobianStiffness(n_max=100)
    assert np.array_equal(eng._traces[0], np.array(ORBIT_TRACE_F1, dtype=float))
    assert np.array_equal(eng._traces[1], np.array(ORBIT_TRACE_F2, dtype=float))
    assert np.array_equal(eng._traces[2], np.array(ORBIT_TRACE_F3, dtype=float))


def test_no_phantom_polynomial_roots():
    # The v0.6.1 canon: all six roots of the f3 T2 charpoly are real.
    roots = np.roots(x0_143.hecke_polynomial_f3_T2())
    assert np.all(np.abs(roots.imag) < 1e-9)
