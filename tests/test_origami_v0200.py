"""v0.20.0 — origami tier and the Hardy-Ramanujan benchmark.

Fast tier runs in seconds; the symbolic-heavy and solver gates are marked slow.
"""
import numpy as np
import pytest
import sympy as sp

from mtft.origami import gates as G
from mtft.origami import (DimerGraph, ensemble_conservation, prism_36,
                          PRISM_C, PRISM_LAMBDA0)
from mtft import hardy_ramanujan as HR


# ------------------------------------------------------------------ fast
def test_version_triple():
    """Genuinely three-way: __init__, pyproject.toml and CITATION.cff.

    Through v0.23.0 this gate only read ``mtft.__version__``, so a stale
    pyproject or CITATION could ship undetected.  It now compares all three
    sources when they are present (they are in a source tree / sdist, and
    absent from a wheel install, which is skipped).
    """
    import pathlib
    import re

    import mtft

    expected = "0.26.0"
    assert mtft.__version__ == expected

    root = pathlib.Path(__file__).resolve().parents[1]
    pyproject = root / "pyproject.toml"
    citation = root / "CITATION.cff"
    if not (pyproject.exists() and citation.exists()):
        import pytest as _pytest
        _pytest.skip("not a source tree; version triple not checkable")

    pv = re.search(r'^version\s*=\s*"([^"]+)"', pyproject.read_text(),
                   re.M).group(1)
    cv = re.search(r'^version:\s*(\S+)', citation.read_text(), re.M).group(1)
    assert pv == expected, f"pyproject.toml says {pv}"
    assert cv == expected, f"CITATION.cff says {cv}"


def test_boundary_measurement_24():
    D = G.gate_boundary_measurement_24()
    p, q, r, s = sp.symbols("p q r s", positive=True)
    assert D[(0, 2)] == sp.expand(p * r + q * s)


def test_mandelstam_channels_24():
    """S(1,3)=pr and S(2,4)=qs: the two terms of the Plucker exchange."""
    G.gate_mandelstam_24(draws=3)


def test_ensemble_conservation_catches_the_bug():
    """The (0,0) class of section B is (1+q)(1+s), not 1+q+s."""
    clsA, clsB = G.gate_ensemble_conservation_24()
    q, s = sp.symbols("q s", positive=True)
    assert sp.factor(clsB[(0, 0)]) == sp.factor((1 + q) * (1 + s))
    assert sp.expand(clsB[(0, 0)] - (1 + q + s)) == sp.expand(q * s)


def test_prism_top_cell_and_symmetry():
    P, cnt = G.gate_prism_top_cell()
    assert len(P) == 20 and cnt == 30   # 6 * C(5,4) three-term relations
    G.gate_prism_symmetry()


def test_prism_partition_function():
    assert prism_36().partition_function() == 280
    assert prism_36(heavy=1).partition_function() == 91


def test_square_mechanism():
    fac, balance = G.gate_square_mechanism()
    w = sp.Symbol("w", positive=True)
    assert fac == sp.factor((2 * w + 3) * (w ** 2 + 3))
    assert balance == [0, 2]


def test_hr_modularity_is_an_identity():
    """eta modularity holds as an identity, not an asymptotic."""
    for b in ("0.1", "0.02"):
        assert HR.modularity_residual(b, dps=30) < 1e-25


def test_hr_prefactor_identity():
    expr, resid = HR.prefactor_identity()
    assert resid == 0
    n = sp.Symbol("n", positive=True)
    assert sp.simplify(expr - 1 / (4 * sp.sqrt(3) * n)) == 0


def test_dimer_graph_is_general():
    """The machinery is not hard-wired to the two instances."""
    color = {"a": "w", "b": "b", "u0": "b", "u1": "w"}
    g = DimerGraph(color, [("a", "b", 2), ("a", "u0", 1), ("u1", "b", 1)],
                   ["u0", "u1"])
    assert g.partition_function() == sp.Integer(3)
    ensemble_conservation(g, [0])


# ------------------------------------------------------------------ slow
@pytest.mark.slow
def test_simplex_curvature_exact():
    assert G.gate_simplex_curvature() == sp.Rational(1, 4)


@pytest.mark.slow
def test_closed_forms_both_routes():
    G.gate_closed_forms_B()


@pytest.mark.slow
def test_spinflip_parity_and_flatness_by_symmetry():
    rho = G.gate_spinflip_parity()
    c = sp.Symbol("c", positive=True)
    assert sp.simplify(rho.subs(c, sp.Rational(9, 4))) == sp.Rational(1, 5)


@pytest.mark.slow
def test_curvature_bounds_strict():
    lo, hi = G.gate_curvature_bounds()
    c = sp.Symbol("c", positive=True)
    assert sp.simplify(lo.subs(c, sp.Rational(9, 4))) == sp.Rational(-5, 16)
    assert sp.simplify(hi.subs(c, sp.Rational(9, 4))) == sp.Rational(5, 36)


@pytest.mark.slow
def test_fisher_quadric_clifford_torus():
    G.gate_fisher_quadric()


@pytest.mark.slow
def test_path_independence_gate():
    out = G.gate_path_independence()
    assert len(out) == 3
    assert all(res < 1e-20 for _v, _t, res in out)


@pytest.mark.slow
def test_prism_lambda0_exact():
    G.gate_prism_lambda0()


@pytest.mark.slow
def test_prism_four_branches_1_plus_3():
    br, fixed, cycles = G.gate_prism_branches(n_starts=800)
    assert len(br) == 4 and len(fixed) == 1 and len(cycles[0]) == 3


@pytest.mark.slow
def test_hr_saddle_beats_closed_form():
    from sympy.functions.combinatorial.numbers import partition
    for n in (50, 500):
        val, _bs, _lay = HR.saddle_partition(n, dps=25)
        ex = int(partition(n))
        hr = HR.hardy_ramanujan_asymptotic(n)
        assert abs(val / ex - 1) < abs(hr / ex - 1)
        assert abs(val / ex - 1) < 0.03
