"""PR-01: seals verify against the code as it stands, evaluation refuses a broken seal, the trials counter is complete and honest."""
import contextlib
import io
import json

import pytest

from mtft.matter import ledger as L
from mtft.matter import prereg as P
from mtft.matter import puzzle as Z
from mtft.matter import star_quotients as S
from mtft.matter.__main__ import main as matter_main


def test_seals_verify_and_are_dated_before_evaluation():
    assert P.verify() == {"PR-S1": "superseded", "PR-S2": True, "PR-S3": True, "PR-S4": "superseded", "PR-S5": True}
    s1, s2, s3, s4, s5 = P.SEALS
    assert s3["sealed"] == s4["sealed"] == s5["sealed"] == "2026-10-11" and s3["evaluated"] == s5["evaluated"] == "2026-10-11"   # GF-01: sealed, then run the same day
    assert s3["outcome"] == "mixed" and s5["outcome"] == "pass" and P.seal_hash(s3) == s3["sha256"] and P.seal_hash(s5) == s5["sha256"]
    assert s4["status"].startswith("superseded") and "CC-64" in s4["status"] and P.seal_hash(s4) != s4["sha256"]   # the constant changed: the old seal no longer matches, by design
    with pytest.raises(RuntimeError):
        P.evaluate("PR-S4")                                                                              # superseded: refused
    assert s1["status"].startswith("superseded") and s1["outcome"] == "superseded" and s1["evaluated"] is None and "CC-63" in s1["status"]
    assert s2["outcome"] == "pending" and s2["evaluated"] is None and s2["sealed"] == "2026-10-10" and s2["status"] == "sealed"
    assert P.seal_hash(s2) == s2["sha256"] and len(s2["sha256"]) == 64 and s1["sha256"] != s2["sha256"]
    assert P.seal_hash(s1) != s1["sha256"]                                                          # the rule changed: the old seal no longer matches, by design
    assert P.seal_hash(dict(s2, params={"tolerance": 1e-3})) != s2["sha256"]                       # any change to the terms changes the seal
    assert P.LEGACY == tuple(t["id"] for t in Z.PREREG if t["id"] not in ("T6", "T7", "T8"))      # the unsealed ones are the prose-only tests
    with pytest.raises(ValueError):
        P.verify("PR-S9")


def test_evaluate_runs_the_sealed_rule_and_refuses_a_broken_seal(monkeypatch):
    r = P.evaluate("PR-S2", ("cm", (1, 0, 143)), det_phase=0.0)
    assert r["result"]["verdict"] == "PASS" and r["sealed"] == "2026-10-10"
    assert P.evaluate("PR-S2", ("cm", (1, 1, 36)), det_phase=0.0)["result"]["verdict"] == "FAIL: theta-bar = pi"
    assert P.evaluate("PR-S2", ("cm", (1, 1, 36)), det_phase=3.141592653589793)["result"]["verdict"] == "PASS"
    with pytest.raises(RuntimeError):
        P.evaluate("PR-S1", ("cm", (1, 0, 143)))                                                   # superseded: refused
    monkeypatch.setitem(P.SEALS[1], "sha256", "0" * 64)
    with pytest.raises(RuntimeError):
        P.evaluate("PR-S2", ("cm", (1, 0, 143)))
    with pytest.raises(ValueError):
        P.evaluate("PR-S9")


def test_trials_counter():
    t = {e["id"]: e for e in P.trials()}
    assert len(t) == 17 and [e["id"] for e in P.trials()] == sorted(t)                              # append-only, ordered
    assert t["X17"]["selects"] == [118, 141, 142, 143] and t["X17"]["source"] == "GF-01" and t["X17"]["sealed"] is None
    assert [i for i, e in t.items() if e["singles_out_143"]] == ["X01", "X04", "X06", "X11"]
    assert t["X10"]["selects"] == [57] and t["X05"]["selects"] == [] and t["X02"]["selects"] == [77, 143] and t["X07"]["selects"] == [123, 143]
    assert t["X13"]["selects"] == [57, 58, 65, 77, 91, 118, 123, 143, 155] and t["X11"]["selects"] == [143]
    assert t["X12"]["selects"] == [57, 77, 91, 123, 141, 143, 155]                                   # a_2 even (CC-63): 57a1 and 91a1 have a_2 = -2
    assert len(t["X06"]["levels_unique_in_their_own_version"]) == 8                                 # eight of the thirteen have a unique orbit pattern
    c = P.counter()
    assert c["trials"] == 17 and c["pre_registered"] == 0 and c["single_level"] == 5 and c["single_out_143"] == 4
    assert c["include_143"] == 15 and c["expected_include_143_by_chance"] == 5.77 and c["expected_single_out_143_by_chance"] == 0.38
    assert "not evidence" in c["reading"]
    # every distinction in star_quotients is logged as a trial (or is the arena itself)
    for key in S.distinctions(143):
        assert key in P.DISTINCTION_TRIALS and (P.DISTINCTION_TRIALS[key] is None or P.DISTINCTION_TRIALS[key] in t)
    assert S.distinctions(143)["only_twin_torus"] and t["X01"]["singles_out_143"]


def test_command_line():
    def run(*argv):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = matter_main(["prereg", *argv])
        return code, buf.getvalue()
    code, out = run(); assert code == 0 and "PR-S1" in out and "17 tried (0 pre-registered)" in out and "X17" in out
    code, out = run("--json"); d = json.loads(out); assert code == 0 and d["counter"]["single_out_143"] == 4 and d["seals"][1]["verifies"] is True
    code, out = run("--verify"); assert code == 0 and "PR-S2 True" in out and "PR-S1 superseded" in out and "PR-S5 True" in out and "PR-S4 superseded" in out


def test_the_record():
    pr = [e for e in L.E if e["wave"] == "PR-01"]
    assert [e["id"] for e in pr] == ["PR-1", "PR-2", "PR-3"] and all(e["module"] == "mtft.matter.prereg" for e in pr)
    cc = next(e for e in L.E if e["id"] == "CC-63")
    assert cc["wave"] == "correction" and "PR-S2" in cc["statement"] and "X12" in cc["statement"]
    a = Z.assignment()
    assert all(a[f"PR-{i}"] == "IV" for i in (1, 2, 3))
    sub = Z.subs()["IV.6"][1]
    assert "PR-3" in sub.claims and sub.status == "open"
