"""CC-67 and CC-68 (2026-10-16): the hypotheses capture under Hubble friction needs, and the divisor series' second pole and summatory remainder.

CC-67: an integrable H leaves energy behind (Gronwall), a fast-growing H freezes the field off the floor (Codex's exact solution), and a bounded H
of bounded variation with a divergent integral captures (the theorem; the radiation and matter eras decay like t^-1.5 and t^-2).  CC-68: Z_D has a
double pole at s = 0 (two routes), S(x) = sum_{n<=x} w_n lies in -(1/2) log^2 x - log x - 2 <= S + zeta'(2) x <= 0, jumps by the unbounded w_n, and
with a smooth weight the s = 0 terms are exact up to O(1/x)."""
import math

import numpy as np
import pytest

from mtft import divisor_gas as D
from mtft.legend import REGISTRY
from mtft.matter import ledger as L
from mtft.matter import puzzle as Z
from mtft.matter import vacuum_atlas as VA


@pytest.fixture(scope="module")
def capture():
    return VA.capture_conditions()


def test_integrable_damping_leaves_energy(capture):
    r = capture["integrable"]
    floor = math.exp(-6 * 0.1 * (1 - math.exp(-200.0)))                          # Gronwall: E(t)/E(0) >= exp(-6 int_0^t H), H = 0.1 e^-t
    assert r["bound_holds"] and r["fraction_left"] > floor and r["fraction_left"] == pytest.approx(0.7682, abs=2e-3)   # no capture: 77% of E stays


def test_fast_growing_damping_freezes_the_field(capture):
    r = capture["overdamped"]
    assert r["exact_solution"] and r["x_at_12"] == pytest.approx(1.0 + 0.5 * math.exp(-12), abs=1e-6)          # phi = a + b e^-t -> a = 1, not 0
    assert 0.3 < r["energy_left"] < 0.4                                                                        # (cosh 1 - 1)/E(0) stays


def test_bounded_nonintegrable_damping_captures(capture):
    c = capture["captured"]
    assert c["constant H = 0.2"]["E_end_over_E0"] < 1e-12
    for name, p in (("radiation H = 1/(2t)", 0.5), ("matter H = 2/(3t)", 2 / 3)):
        assert c[name]["E_end_over_E0"] < 1e-4 and c[name]["decay_exponent"] == pytest.approx(-3 * p, rel=0.02)  # E ~ t^(-3p)


def test_int_H_E_converges_as_the_proof_requires():
    x = np.linspace(-30, 30, 600001)
    assert np.all(x * np.sinh(x) - (np.cosh(x) - 1) >= -1e-9)                    # the convexity inequality x V'(x) >= V(x) the proof uses
    from scipy.integrate import solve_ivp                                        # its conclusion: int H E dt converges in the radiation era
    sol = solve_ivp(lambda t, y: [y[1], -1.5 / t * y[1] - math.sinh(y[0])], (1.25, 2000.0), [1.5, 0.0], rtol=1e-10, atol=1e-12, dense_output=True)
    ts = np.linspace(1.25, 2000.0, 400001); y = sol.sol(ts); E = 0.5 * y[1] ** 2 + np.cosh(y[0]) - 1
    HE = E / (2 * ts); cum = np.concatenate([[0.0], np.cumsum((HE[1:] + HE[:-1]) / 2 * np.diff(ts))])
    i1 = np.searchsorted(ts, 1000.0)
    assert cum[-1] - cum[i1] < 2e-3 * cum[-1]


def test_divisor_series_pole_at_zero():
    p = D.pole_at_zero()
    assert p["double_pole_minus_half"] and p["simple_part_minus_half_log_2pi"] and p["two_routes_agree"] and p["no_triple_pole"]
    assert p["c0"] == pytest.approx(-0.9667703052124541, abs=1e-12) and p["residue_at_1"] == pytest.approx(0.9375482543158438, abs=1e-12)


def test_summatory_function_bounds_and_jumps():
    b = D.summatory_bounds(200000)
    assert b["upper_bound_holds"] and b["lower_bound_holds"]                                                     # -(1/2)log^2 x - log x - 2 <= S + zeta'(2)x <= 0
    assert b["jumps_exceed_bound"] and b["jumps_grow"] and b["jumps"][-1]["w_m_factorial"] > 9                 # w_{12!} = 9.43: unbounded jumps
    assert abs(b["mean_over_x"]) < 1e-3 < b["riesz_bound"]                                                       # the s = 0 terms are the mean: O(X^-1/8)
    assert -6 < b["remainder_min"] < 0 < b["remainder_max"] < 6                                                   # pointwise it wanders, by the jumps


def test_smoothed_sum_is_exact_to_order_one_over_x():
    for x in (500.0, 2000.0):
        s = D.smoothed_sum(x)
        assert s["relative_difference"] < 1e-12                                  # -zeta'(2) x + Res_0 + (ln 2 pi)/(24 x), error O(x^-3)


def test_record_legend_and_map():
    e = {x["id"]: x for x in L.E}
    for cc, mod in (("CC-67", "mtft.matter.vacuum_atlas"), ("CC-68", "mtft.divisor_gas")):
        assert e[cc]["status"] == "active" and e[cc]["module"] == mod and e[cc]["test"] == "test_pg01_corrections.py" and e[cc]["class"] == "EXACT"
    assert e["PG-5"]["status"] == "caveat" and e["PG-5"]["note"].startswith("CC-67") and "for any H(t) > 0" in e["PG-5"]["statement"]   # append-only
    assert e["CC-65"]["status"] == "caveat" and e["CC-65"]["note"].startswith("CC-68") and "+ O(1)" in e["CC-65"]["statement"]
    assert "PG-5" in e["CC-67"]["statement"] and "CC-65" in e["CC-68"]["statement"] and "double pole at s = 0" in e["CC-68"]["statement"]
    a = Z.assignment(); assert a["CC-67"] == "IV" and a["CC-68"] == "IV" and all(Z.check().values()) and all(L.check().values())
    assert "CC-67" in Z.subs()["IV.7"][1].claims and "bounded variation" in Z.subs()["IV.7"][1].answer
    two, mean = REGISTRY["two_primon_energies"].nature, REGISTRY["mean_weight_glaisher"].nature
    assert "+ O(1)" not in two and "O(log^2 x)" in two and "CC-68" in two and "keeps every pole" in two      # the claim is gone, the bound is stated
    assert "not proved" in mean and "CC-68" in mean and "-1/zeta(-1)" in mean
    assert "entire except at s = 1" not in D.__doc__.replace("'entire except at s = 1'", "") and "whatever H(t) > 0" not in VA.energy_argument.__doc__
