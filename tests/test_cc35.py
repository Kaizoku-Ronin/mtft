"""CC-35 (2026-09-29): reflection positivity of the MTFT lattice action retracted from PROVEN to OPEN.  Every number in
docs/SM/CC35_REFLECTION_POSITIVITY.md is tied to the code here."""
import itertools
import numpy as np
import sympy as sp
import mtft
from mtft.research import reflection_positivity as RP


def test_power_sums_are_virtual_characters_exact():
    x = sp.symbols("x1:4")
    def schur(lam):                                     # bialternant formula in three variables
        lam = list(lam) + [0] * (3 - len(lam))
        num = sp.Matrix(3, 3, lambda i, j: x[i] ** (lam[j] + 2 - j)); den = sp.Matrix(3, 3, lambda i, j: x[i] ** (2 - j))
        return sp.cancel(num.det() / den.det())
    for n in range(1, 6):
        p = sum(xi ** n for xi in x)
        rhs = sum(sgn * schur(part) for part, sgn in RP.hook_expansion(n) if len(part) <= 3)
        assert sp.expand(p - rhs) == 0
    assert RP.power_sum_characters(2, 3) == {(2, 0): 1, (1, 1): -1}
    assert RP.power_sum_characters(3, 3) == {(3, 0): 1, (2, 1): -1, (0, 0): 1}
    assert RP.power_sum_characters(3, 2) == {(3,): 1, (1,): -1}


def test_first_order_coefficients_exact_and_negative():
    a = RP.a_n()
    fo3 = RP.first_order_coefficients(3); fo2 = RP.first_order_coefficients(2)
    assert np.isclose(fo3[(1, 1)], (a[3] - a[1]) / 6) and fo3[(1, 1)] < 0
    assert np.isclose(fo3[(2, 1)], -a[2] / 3) and fo3[(2, 1)] < 0
    assert np.isclose(fo2[(1,)], -a[2] / 2) and fo2[(1,)] < 0
    assert float(mtft.weight_array(1)[0]) == 0.0                       # w_1 = 0: the n = 1 term is absent
    assert RP.first_order_coefficients(3, y=0.05)[(1, 1)] > 0 > RP.first_order_coefficients(3, y=0.06)[(1, 1)]   # 3bar sign flips at ln(sqrt 2)/2pi


def test_counterexample_at_package_defaults():
    from mtft.lattice import MTFTAction
    act = MTFTAction()
    assert (act.kappa, act.y) == (RP.KAPPA_DEFAULT, RP.Y_DEFAULT)
    c3 = RP.kernel_coefficients(3, reps=[(1, 1), (2, 1)]); c2 = RP.kernel_coefficients(2, reps=[(1,)])
    assert np.isclose(c3[(1, 1)], -0.004698, atol=2e-6) and np.isclose(c3[(2, 1)], -0.003969, atol=2e-6) and np.isclose(c2[(1,)], -0.005944, atol=2e-6)
    assert RP.site_reflection_pairing(3, (2, 1)) < 0 and RP.site_reflection_pairing(3, (1, 1)) < 0 and RP.site_reflection_pairing(2, (1,)) < 0
    assert all(RP.min_coefficient(3, 1.0, y)[0] < 0 for y in (0.05, 0.18174, 0.3, 1.0))


def test_grid_convergence():
    a = RP.kernel_coefficients(3, reps=[(1, 1)], M=120)[(1, 1)]; b = RP.kernel_coefficients(3, reps=[(1, 1)], M=240)[(1, 1)]
    assert abs(a - b) < 1e-12


def test_positivity_thresholds_numerical():
    lo2, hi2 = RP.positivity_bracket(2, iters=10); assert 58.5 < lo2 < hi2 < 60.0
    lo3, hi3 = RP.positivity_bracket(3, iters=6); assert 650 < lo3 < hi3 < 690


def test_repair_is_positive_definite():
    for N in (2, 3):
        assert all(RP.min_coefficient(N, k, weights="sym")[0] > -1e-12 for k in (0.01, 0.1, 1, 10, 100, 1000))


def test_docstrings_downgraded():
    import mtft.lattice as L, mtft.arithmetic as A, mtft.tower as T
    assert "UNPROVEN" in L.__doc__ and "proven via Osterwalder-Seiler (Paper 24)" not in L.__doc__ and "CC-35" in L.MTFTAction.mass_gap_bound.__doc__
    assert "CC-35" in A.mass_gap_stiffness.__doc__ and "CC-35" in T.__doc__
    from mtft.legend import REGISTRY
    assert "cc35_reflection_positivity" in REGISTRY and "CC-35" in REGISTRY["mu_stiffness"].nature
