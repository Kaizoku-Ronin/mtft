"""GB-02 (2026-11-16): E5, the lookup of the two sealed forecasts of the constant twist energy, and the check of the lab's
transcription against the paper.  The numbers "as read" were read from arXiv:2007.06422v2 by automated extraction on that day
(mtft.maass.glueballs.E5_READ); everything computed here is computed from them."""
from pathlib import Path

import pytest

from mtft.maass import glueballs as G
from mtft.matter import ledger as L
from mtft.matter import puzzle as Z

ROOT = Path(__file__).resolve().parents[1]


def test_ledger_entries_and_the_map():
    E = {e["id"]: e for e in L.E}
    assert [e["id"] for e in L.E if e["wave"] == "GB-02"] == ["GB-6", "GB-7"]
    assert all(E[i]["module"] == "mtft.maass.glueballs" and E[i]["test"] == "test_gb02_e5_sealed_forecasts.py" for i in ("GB-6", "GB-7"))
    assert E["GB-6"]["status"] == "caveat" and "weak support and no more" in E["GB-6"]["statement"] and "little power" in E["GB-6"]["note"]
    assert E["GB-7"]["status"] == "active" and "4+- ground state" in E["GB-7"]["statement"] and "does not help H24" in E["GB-7"]["note"]
    assert "GB-6 and GB-7" in E["GB-5"]["note"] and E["GB-5"]["status"] == "caveat"       # GB-5 is kept as written, with a pointer
    assert Z.WAVE_DATES["GB-02"] == "2026-11-16" and Z.WAVE_REPORT["GB-02"] == "GB02_E5_SEALED_FORECASTS.md"
    assert (ROOT / "docs" / "SM" / Z.WAVE_REPORT["GB-02"]).exists()
    n24 = [n for n in Z.NEGATIVES if n["id"] == "N24"]
    assert n24[0]["claims"] == ("GB-1", "GB-2", "GB-3", "GB-4")                           # a pass, however weak, is not a negative


def test_the_sealed_rule_on_the_value_read():
    seal = G.SEALS[0]
    assert G.verify_seals() == {"GB-S1": True} and seal["sha256"] == "e75106afbfc3e3bd92c29ebcfeea81db48a0c67ab715c50fdc29dad1735500e7"
    assert G.EXCITED == {"A1--": (10.39, 0.37), "A1+-": None}
    got = G.e5()
    res = got["sealed"]["result"]
    odd, even = res["A1-- first excited (second odd)"], res["A1+- first excited (second even)"]
    assert odd["verdict"] == "PASS" and odd["sigma"] == pytest.approx(1.4121, abs=1e-4)
    assert odd["lattice_ratio"] == pytest.approx(10.39 / 3.405) and odd["lattice_error"] == pytest.approx(0.1103, abs=1e-4) and odd["forecast"] == 2.8658
    assert even == {"forecast": 3.2399, "verdict": "NO DATA"} and res["overall"] == "PASS (one channel)"
    assert seal["outcome"] == "pass in one channel, no data in the other" and seal["evaluated"] == "2026-11-16"
    # the overall verdict does not depend on which C = -1 column the excited entry stands in
    assert got["if the column were the other"] == "PASS (one channel)" and got["sigma if the column were the other"] == pytest.approx(1.4342, abs=1e-4)
    assert G.summary()["e5"] == seal["outcome"]
    by = {r["id"]: r for r in G.register()}
    assert "1.41 sigma" in by["GB-S1"]["result"] and "little power" in by["GB-S1"]["result"]


def test_the_six_transcribed_values_agree_with_the_paper():
    check = G.transcription_check()
    assert set(check) == {*G.LATTICE, "table 17, J = 0"} and all(check.values())
    assert [len(row) for row in G.TABLE14_A1] == [4, 3, 2, 2]                              # as read: the short rows carry no columns
    c_minus = [v for row in G.TABLE14_A1[1:] for v in row if v[0] > 9.5]                   # beyond the two ground states, one C = -1 entry
    assert c_minus == [(10.39, 0.37)] and "4+-" in G.TRANSCRIPTION_CHECK
    assert set(G.LATTICE) == {"0++", "0++*", "0-+", "0-+*", "A1+-", "A1--"}                 # the frozen table is not edited


def test_the_refit_identifies_the_column():
    """The paper's continuum extrapolation, redone from its own finite-spacing tables: it reproduces the two controls, and
    it gives 10.39(37) for the P=-,C=- states listed first at the four finest spacings and for nothing under P=+,C=-."""
    fits = G.identify_column()
    a, b = fits["control A1++ gs (paper: 3.405(21))"], fits["control A1+- gs, beta >= 5.99 (paper: 9.32(28))"]
    assert a["value"] == pytest.approx(3.405, abs=0.003) and round(a["error"], 3) == 0.021
    assert round(b["value"], 2) == 9.32 and round(b["error"], 2) == 0.28
    hit = fits["P=-,C=- listed first, four finest"]
    assert round(hit["value"], 2) == 10.39 and hit["error"] == pytest.approx(0.37, abs=0.01) and hit["chi2"] / hit["dof"] < 2.5
    second = fits["P=-,C=- listed second, four finest"]                                    # the other A1-- state, beside the paper's 9.78(46)
    assert second["value"] == pytest.approx(9.78, abs=0.2) and round(second["error"], 2) == 0.46
    for key in ("P=+,C=- listed second, four finest", "P=+,C=- second lightest, four finest"):
        assert fits[key]["chi2"] / fits[key]["dof"] > 2.5                                  # the paper stars such fits; 10.39(37) has no star
    assert fits["P=+,C=- listed second, four finest"]["value"] > 11.2
    assert fits["matches 10.39(37) with an unstarred fit"] == ["P=-,C=- listed first, four finest"] and fits["column"] == "P=-,C=-"
    # the same authors' later SU(3) table, whose cells are kept: the three other columns digit for digit, a second entry only under P=-,C=-
    later = G.LATER["A1 block"]
    assert later[0][:3] == G.TABLE14_A1[0][:3] and later[1][:2] == G.TABLE14_A1[1][:2] and later[1][2] is None and later[1][3] == (10.03, 0.47)
    assert all(row[2] is None and row[3] is None for row in later[2:])
    assert {later[2][0], later[3][0]} == {G.TABLE14_A1[2][0], G.TABLE14_A1[3][0]}


def test_the_column_as_read_from_the_pdf():
    """Added after the lookup: the PDF's page was read (Roger Tano's screenshot of Table 14, and the file's text with its layout
    kept).  The columns of the A1 block, and that they are the flattened rows of the first reading and the column identified."""
    cols = G.TABLE14_A1_COLUMNS
    assert list(cols) == ["P=+,C=+", "P=-,C=+", "P=+,C=-", "P=-,C=-"] and [len(v) for v in cols.values()] == [4, 4, 1, 2]
    rows = tuple(tuple(c[k] for c in cols.values() if len(c) > k) for k in range(4))
    assert rows == G.TABLE14_A1                                                            # the flattened text was the table read row by row
    assert cols["P=-,C=-"][1] == G.EXCITED["A1--"] and cols["P=+,C=-"] == (G.LATTICE["A1+-"][:2],) and G.EXCITED["A1+-"] is None
    assert G.identify_column()["column"] == "P=-,C=-" and "P=-,C=-" in G.EXCITED_READ and "confirmed" in G.EXCITED_READ
    assert G.e5()["column read"] == G.EXCITED_READ
    assert G.A1_PLUS_MINUS_ALL[6.235] == ((1.431, 0.020), (1.482, 0.025), (1.77, 0.04))    # the third state, which the flattened text left without a column


def test_what_the_pass_is_worth():
    got = G.e5()
    lo, hi = got["passing window at this error"]
    assert (lo, hi) == pytest.approx((8.866, 10.653), abs=0.002) and lo < got["ground state"][0] < hi   # only a state above 10.65 would have failed
    assert got["outside the seal: the authors' later value"]["value"] == (10.03, 0.47)
    assert got["outside the seal: the authors' later value"]["sigma"] == pytest.approx(0.5103, abs=1e-4)
    spin = G.SPIN_OF_A1["A1+- gs"]
    assert spin["assigned"] == "4+- gs" and spin["stars"] == 2 and G.SPIN_OF_A1["A1-- gs"] is None and G.SPIN_OF_A1["A1-- ex1"] is None
    assert spin["members"][0] == G.LATTICE["A1+-"][:2]                                     # the state the lab used as J = 0, C = -1
    assert round(got["mean of the four 4+- members"], 2) == spin["table 17"][0] == 9.02    # the paper's 4+- is the mean of the four
    assert "weak support at most" in got["worth"] and "spin 4" in got["worth"]
    tw = G.twist_energy()                                                                   # what the fit had used: both ground states as J = 0
    assert tw["chi2"] == pytest.approx(1.5697, abs=1e-3)
