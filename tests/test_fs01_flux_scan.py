"""FS-01: stack-level flux scan — curve degrees suffice, TK-1 lifted, the mod-4 Higgs-charge rule, the ratio-3 / charged-nu^c separation."""
from fractions import Fraction
import pytest
from mtft.matter import flux_scan as F


@pytest.fixture(scope="module")
def rows():
    return F.scan(3, 3)


def test_curve_degrees_alone_suffice():
    s = F.StackFlux((0, 1, 1, -1, -2), (0, 0, 0, 0, 0)).summary()
    assert s and s["eps"] == Fraction(1, 3) and s["kappa"] == 0 and s["pair"] == ("A1", "A2") and s["fermions"] == 288
    s4 = F.StackFlux((1, 1, 1, -1, -5), (0, 0, 0, 0, 0)).summary(); assert s4 and s4["su4"] and s4["pair"] == ("4", "A2")
    for s in (s, s4):
        mu = s["mu"] if "mu" in s else None
    k = F.StackFlux((0, 1, 1, -1, -2), (0, 0, 0, 0, 0)).solve_kahler(); assert k["mu"]["L"] == 0 and k["mu"]["B"] == 0 and k["mu"]["c"] == 0 and k["mu"]["A1"] == -k["mu"]["A2"] != 0


def test_tk1_lifted_and_families(rows):
    assert any(r["su4"] and abs(r["charges"]["nu_c"]) == 1 and r["charges"]["delta_d"] == 0 for r in rows)
    fam = F.families(rows=rows)
    assert fam["minimal (curve degrees only)"]["fermions"] == 288
    assert fam["ratio 3 (SU(4) kept, Higgs charges 1 and 3)"]["fermions"] == 672
    assert fam["neutrino-optimal, both Higgs localised"]["fermions"] == 416 and fam["neutrino-optimal, up Higgs global"]["fermions"] == 1056


def test_higgs_charge_rule(rows):
    assert all((r["charges"]["delta_u"] + r["charges"]["delta_d"]) % 4 == 0 for r in rows if r["su4"])
    rule = F.higgs_charge_rule()
    assert rule[Fraction(3)] == (1, 3) and rule[Fraction(5, 3)] == (3, 5) and rule[Fraction(2)] == (4, 8) and rule[Fraction(7, 3)] == (6, 14)
    both = {round(r["eps_ratio"], 6) for r in rows if r["su4"] and r["charges"]["delta_u"] and r["charges"]["delta_d"]}
    assert 3.0 in both and 2.0 not in both


def test_ratio_three_with_charged_nu_c(rows):
    assert not any(r["su4"] and abs(r["charges"]["nu_c"]) == 1 and abs(r["eps_ratio"] - 3) < 1e-9 for r in rows)   # not in the box
    w = F.ratio_three_with_charged_nu_c()
    assert w and w["su4"] and abs(w["charges"]["nu_c"]) == 1 and w["eps_ratio"] == 3.0 and w["fermions"] == 3104
