"""MA-01 (2026-11-15): the first rung of the curve's tower, by two independent routes; the ledger entries.

Route A is Hejhal's method (mtft.maass, frozen with provenance and found again by the solver in test_maass_solver.py).
Route B is SPEC-01's finite elements (mtft.surface.spectral), which share no step with it.
"""
import math

from mtft.maass import spectrum as SP
from mtft.matter import ledger as L
from mtft.matter import puzzle as Z
from mtft.surface import spectral as FE


def test_finite_elements_reproduce_the_first_nine_eigenvalues():
    hejhal = sorted(f["lambda"] for f in SP.forms(twist=False))[:9]
    assert abs(hejhal[0] - 0.3897651801) < 1e-9
    r = FE.maass_candidates(143, h=0.22, nx=7, k=14)
    fem = sorted(r["stable"])[:9]
    assert r["class"] == "DIAGNOSTIC" and all(r["gates"].values()) and len(fem) == 9
    assert abs(r["lambda_1"] - hejhal[0]) < 0.02                              # the forecast's window; found 0.0027
    for a, b in zip(hejhal, fem):
        assert 0.0 < b / a - 1.0 < 0.015, (a, b)                              # 0.70 % to 1.05 % above, level by level
    assert abs(hejhal[7] - hejhal[8]) < 1e-9 and abs(fem[7] - fem[8]) < 0.01  # the double eigenvalue, seen by both
    assert hejhal[3] - hejhal[2] > 0.3 and fem[3] - fem[2] > 0.3              # and the same wide step after the third


def test_ledger_entries_and_the_map():
    E = {e["id"]: e for e in L.E}
    assert E["MA-1"]["class"] == "NUMERICAL" and E["MA-1"]["status"] == "active" and E["MA-1"]["module"] == "mtft.maass.spectrum"
    assert E["MA-2"]["status"] == "caveat" and "not yet in the package" in E["MA-2"]["note"]
    assert E["MA-3"]["class"] == "EXACT" and E["MA-3"]["module"] == "mtft.maass.scattering"
    assert "sqrt(0.39)" in E["VE-5"]["statement"] and "DIAGNOSTIC" in E["VE-5"]["statement"]      # VE-5's own text stands
    assert Z.WAVE_DATES["MA-01"] == "2026-11-15" and Z.assignment()["MA-1"] == Z.assignment()["VE-5"] == "IV"
    assert all(L.check().values()) and all(Z.check().values())


def test_the_rung_in_ve5_units():
    g = SP.gap()
    assert abs(math.sqrt(g["lambda_1"]) - 0.62431) < 1e-5
    lo, hi = 0.063 * 2 * math.sqrt(g["lambda_1"]), 0.070 * 2 * math.sqrt(g["lambda_1"])       # 1/(2 R_C) = 0.063-0.070 M_KK
    assert abs(lo - 0.079) < 1e-3 and abs(hi - 0.087) < 1e-3                  # VE-5's 0.079-0.087 M_KK, unchanged
