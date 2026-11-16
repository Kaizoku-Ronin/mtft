"""mtft.gravity.ledger (2026-11-16): the gravity work's entries, the map's tests and negatives in its own words, the misses."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import mtft.gravity as GR
from mtft.gravity import ledger as GL

ROOT = Path(__file__).resolve().parents[1]


def test_the_ledger_checks_itself():
    check = GL.check(tests_dir=ROOT / "tests")
    assert len(check) == 11 and all(v is True for v in check.values()), check
    assert all(v is not False for v in GL.check().values())                         # by default: the tests folder of a checkout, if there is one
    assert GL.check(tests_dir=ROOT / "no_such_folder")["tests_exist"] is None      # an installed copy has no tests folder: neither a pass nor a failure


def test_the_check_can_fail(monkeypatch):
    """Each line of the check is shown to fail when what it guards is broken: a check that cannot fail checks nothing."""
    tests = ROOT / "tests"
    assert all(v is True for v in GL.check(tests_dir=tests).values())
    with monkeypatch.context() as m:                                                 # the data file is not the one the ledger was written against
        m.setattr(GL, "DATA_SHA1", "0" * 40)
        assert GL.check(tests_dir=tests)["map_as_parsed"] is False
    with monkeypatch.context() as m:
        m.setattr(GL, "MAP_SHA1", "0" * 40)
        assert GL.check(tests_dir=tests)["map_as_parsed"] is False
    for module in ("mtft.gravity.no_such_module", "mtft.no_such_package.module"):    # a cited module that does not exist; the second raises inside importlib
        with monkeypatch.context() as m:
            m.setitem(GL.E[0], "module", module)
            assert GL.check(tests_dir=tests)["modules_exist"] is False
    with monkeypatch.context() as m:                                                 # a cited test file that does not exist
        m.setitem(GL.E[0], "test", "test_no_such_file.py")
        assert GL.check(tests_dir=tests)["tests_exist"] is False
    with monkeypatch.context() as m:                                                 # a class outside the matter ledger's five
        m.setitem(GL.E[0], "class", "PROVED")
        assert GL.check(tests_dir=tests)["classes_and_statuses"] is False
    with monkeypatch.context() as m:                                                 # a retraction without its pointer
        m.setitem(GL.E[9], "note", "")
        assert GL.E[9]["status"] == "retracted" and GL.check(tests_dir=tests)["retracted_have_pointers"] is False
    with monkeypatch.context() as m:                                                 # a duplicated id
        m.setitem(GL.E[1], "id", "GR-1")
        assert GL.check(tests_dir=tests)["ids_unique"] is False
    with monkeypatch.context() as m:                                                 # a miss said to be among the negatives that is not
        m.setattr(GL, "MISSES_AMONG_NEGATIVES", ("N-G7", "N-G99"))
        assert GL.check(tests_dir=tests)["misses_point_to_negatives"] is False
    with monkeypatch.context() as m:                                                 # an open test dropped from the list
        m.setattr(GL, "OPEN", {"H-G7": "x"})
        assert GL.check(tests_dir=tests)["open_is_listed"] is False
    real = GL._map()                                                                 # the map's records themselves
    for line, breaker in (("sixteen_negatives", lambda m: m["negatives"].pop(3)),
                          ("ten_tests_one_pending", lambda m: m["tests"][2].update(status="done")),
                          ("eleven_pieces", lambda m: m["pieces"].pop())):
        broken = copy.deepcopy(real)
        breaker(broken)
        with monkeypatch.context() as m:
            m.setattr(GL, "_map", lambda broken=broken: broken)
            assert GL.check(tests_dir=tests)[line] is False
    assert len(GL.check(tests_dir=tests)) == 11 and all(v is True for v in GL.check(tests_dir=tests).values())   # all eleven lines shown; nothing left broken


def test_entries():
    E = {e["id"]: e for e in GL.E}
    assert list(E) == [f"GR-{k}" for k in range(1, 13)] and {e["wave"] for e in GL.E} == set(GL.WAVES)
    assert [e["id"] for e in GL.entries(status="retracted")] == ["GR-10", "GR-11", "GR-12"]   # retracted statements stay exposed
    assert "63.9" in E["GR-10"]["note"] and "4.7" in E["GR-11"]["note"] and "69.8" in E["GR-11"]["note"] and "2/3" in E["GR-12"]["note"]
    lab22 = (ROOT / "studies" / "gravity_2026nov" / "lab" / "mp22" / "MP22_TG8_Hopfion_Plan.md").read_text(encoding="utf-8")
    assert "It reaches 4.7." in lab22 and "one turn in 63.9 circuits" in lab22 and "at 0.43 of the range" in lab22 and "69.8" in lab22   # the lab's words
    lab21 = (ROOT / "studies" / "gravity_2026nov" / "lab" / "mp21" / "MP21_TG7_Shear_Plan.md").read_text(encoding="utf-8")
    assert all(x in lab21 for x in ("−15.2404", "−15.0920", "−4.4244", "−4.6719", "−6.018516", "−6.018519"))                          # M-G2's numbers
    assert E["GR-9"]["status"] == "caveat" and "belongs to the frame" in E["GR-9"]["note"] and "AG-D5" in E["GR-9"]["note"]
    for i in ("GR-1", "GR-5", "GR-8"):                                               # whose work is whose, and "not found is not new"
        assert "unreviewed" in E[i]["note"] and "not found in the sources read" in E[i]["note"]
    assert "Kedia" in E["GR-1"]["note"] and "Swearngin" in E["GR-5"]["note"] and "Garnier" in E["GR-5"]["note"] and "not a symbolic proof" in E["GR-8"]["note"]
    assert all((ROOT / "tests" / e["test"]).exists() and (ROOT / GL.REPORT).exists() for e in GL.E)
    entries = GL.entries(wave="GR-02")
    entries[0]["statement"] = "changed"                                              # copies: the ledger cannot be edited through them
    assert GL.E[4]["statement"] != "changed"


def test_the_map_in_its_own_words():
    neg, tests, pieces = GL.negatives(), GL.tests(), GL.pieces()
    assert [n["id"] for n in neg] == [f"N-G{k}" for k in range(1, 17)] and [t["id"] for t in tests] == [f"T-G{k}" for k in range(1, 11)]
    assert [t["id"] for t in tests if t["status"] == "pending"] == ["T-G3"] and sum(t["status"] == "done" for t in tests) == 9
    assert len(pieces) == 11 and {p["status"] for p in pieces} == {"open", "partial"} and [p["piece"] for p in pieces if p["status"] == "partial"] == ["III", "IV"]
    assert neg[0]["statement"].startswith("G is not derived") and neg[0]["pieces"] == "I" and set(GR.NOT_DERIVED) == {"N-G1"}
    assert "63.9" in neg[14]["statement"] and "2/3" in neg[15]["statement"] and "|W|" in neg[15]["statement"]   # a bare | inside a cell is kept
    # the data file is the lab's map as parsed, and the map is in the study folder
    data = ROOT / "src" / "mtft" / "gravity" / "_data" / "gravity_map.json"
    source = ROOT / "studies" / "gravity_2026nov" / "lab" / "MTFT_Gravity_Puzzle_Map.md"
    assert hashlib.sha1(source.read_bytes()).hexdigest() == GL.MAP_SHA1 == json.loads(data.read_text(encoding="utf-8"))["source_sha1"]
    assert hashlib.sha1(data.read_bytes()).hexdigest() == GL.DATA_SHA1 and GL.DATA_SHA1 in (data.parent / "PROVENANCE.txt").read_text(encoding="utf-8")
    text = source.read_text(encoding="utf-8")
    assert all(n["statement"] in text for n in neg) and all(t["criterion"] in text for t in tests)        # carried, not retyped
    # re-derived: the parser run again on the map gives the data file, record for record
    spec = importlib.util.spec_from_file_location("build_gravity_map", ROOT / "studies" / "gravity_2026nov" / "build_gravity_map.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    assert builder.build() == json.loads(data.read_text(encoding="utf-8"))
    assert builder.cells("| a | b |c| | d |") == ["a", "b |c|", "d"]                                      # only ' | ' separates cells
    assert "10 pre-registered tests (1 pending, T-G3; the other nine done), 16 recorded negatives" in text  # the map's own count


def test_misses_and_open():
    assert [m["id"] for m in GL.MISSES] == [f"M-G{k}" for k in range(1, 7)] and sum(m["kind"] == "found in the 1.1.0 build" for m in GL.MISSES) == 3
    neg = {n["id"]: n["statement"] for n in GL.negatives()}
    for i in GL.MISSES_AMONG_NEGATIVES:                                              # each of these says it was wrong or missed
        assert any(w in neg[i] for w in ("wrong", "missed", "Corrected")), i
    assert {"T-G3", "H-G7"} <= set(GL.OPEN) and set(GR.OPEN) <= set(GL.OPEN)


def test_the_document_is_the_ledger():
    doc = (ROOT / "docs" / "gravity" / "GRAVITY_LEDGER.md").read_text(encoding="utf-8")
    assert doc == GL.to_markdown()
    assert "12 entries (8 active, 1 with a caveat, 3 retracted and kept); 10 pre-registered tests (9 done, 1 pending); 16 negatives; 6 further recorded misses." in doc
