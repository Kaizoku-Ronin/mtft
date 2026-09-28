"""The parent-action ledger (docs/SM/PARENT_ACTION_REQUIREMENTS.md) and CC-34: every gate the ledger names exists and returns the verdict the
ledger records, so the document cannot drift from the code.  All checks exact except the CC-34 arithmetic (mpmath, 20 digits)."""
import re
from pathlib import Path
import mpmath as mp
import sympy as sp
from mtft.research import product_surface as PS, gravitational_anomaly as GA, unified_parent as UP

ROOT = Path(__file__).resolve().parents[1]
LEDGER = (ROOT / "docs/SM/PARENT_ACTION_REQUIREMENTS.md").read_text(encoding="utf-8")


def test_every_ledger_gate_exists():
    names = set(re.findall(r"`([a-z_0-9]+)`", LEDGER.split("## B.")[0]))
    modules = (PS, GA, UP)
    missing = [n for n in names if n not in {"cc26_normalisation", "hypercharge_assignment"}
               and not any(hasattr(m, n) for m in modules)]
    assert missing == [], missing


def test_recorded_verdicts():
    assert PS.chirality_parity_theorem(box=3)["untwisted_nets"] == {0}                              # row 1
    assert PS.three_stack_trilemma()["Y+H+C"] == []                                                   # rows 2, 4
    assert PS.block_census(sp.Rational(1, 2), "S2")["summary"]["harmonic_negative"] == 90              # row 6
    assert PS.s1_iterated_extension_instability(2)["unstable"]                                        # row 5
    assert PS.higgs_slope_mass((2, -4), (0, 1))["massless_locus"] == [PS.A_E / 2]                     # row 7
    assert PS.torus_cp_test()["S2"]["real_up_to_rephasing"]                                           # row 8
    assert GA.m1_tensor_survivors()["survivors"] == [] and not UP.gaugino_p2_integrality("E8")["integral"]  # row 9 (CC-33)


def test_cc34_arithmetic():
    mp.mp.dps = 20
    delta = mp.mpf("4.669201609102990671853"); ainv = mp.mpf("137.035999084")
    e = mp.exp(-2 * ainv); f = delta ** -6 * e
    assert abs(e / mp.mpf("9.3766e-120") - 1) < 1e-4 and abs(e / mp.mpf("1.342e-119") - 1) > 0.3        # the papers' 1.342e-119 is a slip
    assert abs(f / mp.mpf("9.0487e-124") - 1) < 1e-4                                                  # the formula's value
