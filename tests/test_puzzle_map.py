"""The programme map (mtft.matter.puzzle): integrity against the ledger, the modules and the reports."""
import contextlib
import importlib
import io
import json
import re
from pathlib import Path

import pytest

from mtft.matter import ledger as L
from mtft.matter import puzzle as Z
from mtft.matter.__main__ import main as matter_main

ROOT = Path(__file__).resolve().parents[1]


def test_map_integrity_and_non_vacuity():
    c = Z.check()
    assert all(c.values()), {k: v for k, v in c.items() if not v}
    assert {"docs_exist", "every_doc_placed"} <= set(c)                      # the source-tree checks ran
    assert len(Z.PIECES) == 10 and sum(len(p.subs) for p in Z.PIECES) >= 50
    assert len(Z.PREREG) >= 6 and len(Z.NEGATIVES) >= 19


def test_every_claim_in_exactly_one_piece():
    a = Z.assignment()
    assert len(a) == len(L.E) == len({e["id"] for e in L.E})
    counts = {pid: len(Z.claims_of(pid)) for pid in Z.ORDER}
    assert sum(counts.values()) == len(L.E)
    assert counts["I"] == 0 and all(counts[p] > 0 for p in Z.ORDER if p != "I")    # the parent lives in the registers, not the ledger
    assert a["LEP-4"] == "X" and a["CC-62"] == "IX" and a["AG-4"] == "IV" and a["AG-1"] == "II" and a["YB-2"] == "IX"


def test_generated_map_in_sync():
    assert (ROOT / "docs" / "SM" / "MATTER_PUZZLE_MAP.md").read_text(encoding="utf-8") == Z.to_markdown()


def test_wave_dates_match_module_docstrings():
    checked = 0
    for wave, date in Z.WAVE_DATES.items():
        for mod in {e["module"] for e in L.E if e["wave"] == wave}:
            doc = importlib.import_module(mod).__doc__ or ""
            m = re.search(r"\d{4}-\d{2}-\d{2}", doc[:200])
            if m:
                assert m.group(0) == date, (wave, mod, m.group(0)); checked += 1
    assert checked >= 50                                                      # non-vacuity: nearly every wave's module is dated


def test_statuses_follow_the_record():
    P = Z.pieces()
    assert [p.status for p in Z.ordered()] == ["open", "stands", "stands", "partial", "partial", "open", "partial", "open", "closed", "closed"]
    for pid in ("IX", "X"):                                                   # closed pieces: only 'reopen' questions stay open
        assert all(s.reopen for s in P[pid].subs if s.status == "open") and any(s.status == "failed" for s in P[pid].subs)
    assert not any(s.status == "failed" for pid in ("II", "III") for s in P[pid].subs)
    E = {e["id"]: e for e in L.E}
    assert all(E[c]["status"] in ("active", "caveat") for n in Z.NEGATIVES for c in n["claims"])
    assert {r["id"] for r in Z.to_json()["retracted"]} == {e["id"] for e in L.E if e["status"] == "retracted"}


def test_preregistered_outcomes_agree_with_the_ledger():
    E = {e["id"]: e["statement"] for e in L.E}
    words = {"T1": "No canonical prescription survives", "T2": "outside", "T3": "PASS", "T4": "FAIL", "T5": "FAIL"}
    T = {t["id"]: t for t in Z.PREREG}
    for tid, w in words.items():
        assert any(w in E[c] for c in T[tid]["claims"]), tid
    assert [t["outcome"] for t in Z.PREREG] == ["mixed", "no hit", "fail", "pass", "fail", "fail", "pending", "pass", "mixed"]
    used = {s.prereg for p in Z.PIECES for s in p.subs if s.prereg}
    assert used == {"T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8"}
    assert "pending" in E["TR-8"] and T["T6"]["outcome"] == "pending"                               # sealed before any candidate: PR-S2 (PR-S1 superseded, CC-63)
    assert "PASS" in E["GF-3"] and T["T7"]["outcome"] == "pass" and "FAIL" in E["GF-2"] and T["T8"]["outcome"] == "mixed"   # GF-01: PR-S4, PR-S3


def test_layers_and_dependencies():
    lay = Z.layers()
    assert lay["I"] == lay["II"] == 0 and lay["IX"] == 5 and lay["X"] == 6
    assert "digraph" in Z.to_dot() and Z.to_dot().count("->") == sum(len(p.depends_on) for p in Z.PIECES)


def test_raise_paths(monkeypatch):
    monkeypatch.setattr(Z, "_module_piece", lambda: {"m": {"II", "IV"}})
    with pytest.raises(ValueError):
        Z.piece_of({"id": "CC-99", "wave": "correction", "module": "m"})
    monkeypatch.undo()
    a = Z.Piece("I", "a", "q", "open", "v", (), (), ("II",), (), ())
    b = Z.Piece("II", "b", "q", "open", "v", (), (), ("I",), (), ())
    monkeypatch.setattr(Z, "PIECES", (a, b)); monkeypatch.setattr(Z, "ORDER", ("I", "II"))
    with pytest.raises(ValueError):
        Z.layers()
    assert Z.check()["dependencies_acyclic"] is False


def test_command_line():
    def run(*argv):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = matter_main(["puzzle", *argv])
        return code, buf.getvalue()
    code, out = run("--json"); d = json.loads(out)
    assert code == 0 and len(d["pieces"]) == 10 and d["ledger_size"] == len(L.E)
    code, out = run("--piece", "IX"); assert code == 0 and "Quark flavour" in out and "IX.9" in out
    code, out = run("--piece", "P7"); assert code == 2 and "unknown piece" in out
    code, out = run("--open"); assert code == 0 and "[answered]" not in out and "[open]" in out
    code, out = run("--markdown"); assert out == Z.to_markdown()
