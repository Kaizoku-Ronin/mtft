"""GB-01 (2026-11-16): the lattice comparison's ledger entries, its negative on the puzzle map, the seal, and the lab's
record of needed shifts, which is kept in the study folder and must agree with the package."""
import json
from pathlib import Path

import pytest

import mtft.maass as MA
from mtft.maass import glueballs as G
from mtft.matter import ledger as L
from mtft.matter import puzzle as Z

RECORD = Path(__file__).resolve().parents[1] / "studies" / "maass_2026nov" / "lab" / "needed_shifts.json"


def test_ledger_entries_and_the_map():
    E = {e["id"]: e for e in L.E}
    ids = [e["id"] for e in L.E if e["wave"] == "GB-01"]
    assert ids == ["GB-1", "GB-2", "GB-3", "GB-4", "GB-5"] and all(E[i]["module"] == "mtft.maass.glueballs" for i in ids)
    assert [E[i]["class"] for i in ids] == ["MODEL", "MODEL", "MODEL", "NUMERICAL", "MODEL"]
    assert [E[i]["status"] for i in ids] == ["active", "caveat", "active", "active", "caveat"]
    assert "All five predictions miss" in E["GB-1"]["statement"] and "registered miss" in E["GB-3"]["statement"]
    assert "not evidence" in E["GB-5"]["statement"] and "lookup has not been made" in E["GB-5"]["note"] and "GB-S1" in E["GB-5"]["statement"]
    assert Z.WAVE_DATES["GB-01"] == "2026-11-16" and Z.WAVE_REPORT["GB-01"] == "GB01_GLUEBALL_COMPARISON.md"
    assert {Z.assignment()[i] for i in ids} == {"IV"} == {Z.assignment()["MA-1"]}                  # beside the spectrum it compares
    n24 = [n for n in Z.NEGATIVES if n["id"] == "N24"]
    assert len(n24) == 1 and n24[0]["claims"] == ("GB-1", "GB-2", "GB-3", "GB-4") and "GB-5" not in n24[0]["claims"]      # the open one is not a negative
    assert all(L.check().values()) and all(Z.check().values())
    assert "glueballs" in MA.__all__


def test_the_ledger_numbers_are_the_module_numbers():
    E = {e["id"]: e["statement"] for e in L.E}
    rb = G.reading_b()
    for k in ("N1", "N2", "N3", "N4"):
        assert f"{rb[k]['value']:.3f}" in E["GB-1"]
    need = G.needed_shifts()
    assert f"{need['A1--']['curve_ratio']:.3f}" in E["GB-3"] and f"{need['A1+-']['curve_ratio']:.3f}" in E["GB-3"]
    assert f"{abs(need['A1--']['sigma']):.1f}" in E["GB-3"] and f"{abs(need['A1+-']['sigma']):.1f}" in E["GB-3"]
    assert f"{G.self_interaction_fit()['chi2']:.1f} on 2" in E["GB-4"] and f"{G.one_coupling_fit()['chi2']:.0f} on 4" in E["GB-4"]
    assert f"({G.one_coupling_fit(correlated=True)['chi2']:.0f} with the shared error" in E["GB-4"]
    tw = G.twist_energy()
    assert f"{tw['mu']:.3f}({round(tw['error'] * 1000)})" in E["GB-5"] and f"{tw['chi2']:.2f} on 1" in E["GB-5"]
    assert f"{tw['second_odd']:.2f}(7)" in E["GB-5"] and f"{tw['second_even']:.2f}(7)" in E["GB-5"]


def test_the_lab_record_of_needed_shifts_agrees():
    lab = json.loads(RECORD.read_text(encoding="utf-8"))
    need = G.needed_shifts()
    assert set(lab) == set(need)
    for state, rec in lab.items():
        d = need[state]
        assert d["lattice_ratio"] == pytest.approx(rec["lattice ratio"], abs=1e-12) and d["error"] == pytest.approx(rec["error"], abs=1e-12)
        assert d["curve_ratio"] == pytest.approx(rec["curve ratio"], abs=1e-8)                    # the package's own eigenvalues
        assert d["shift"] == pytest.approx(rec["needed shift in lambda"], abs=1e-9) and d["sigma"] == pytest.approx(rec["sigma"], abs=1e-6)
        assert d["percent"] == pytest.approx(rec["needed shift, percent of the curve value"], abs=1e-6)


def test_the_seal_is_intact_and_the_table_is_as_it_was_sealed():
    """Until E5 (16 November 2026, later the same day) the outcome was "pending"; the lookup's result is tested in test_gb02."""
    assert G.verify_seals() == {"GB-S1": True}
    seal = G.SEALS[0]
    assert seal["sha256"] == "e75106afbfc3e3bd92c29ebcfeea81db48a0c67ab715c50fdc29dad1735500e7"   # as written into the build record at sealing
    assert seal["outcome"] == "pass in one channel, no data in the other" and seal["evaluated"] == "2026-11-16" and "E5" in seal["result"]
    assert set(G.LATTICE) == {"0++", "0++*", "0-+", "0-+*", "A1+-", "A1--"}                        # no excited C = −1 state is in the table
