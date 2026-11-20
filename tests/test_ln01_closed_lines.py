"""LN-01 (2026-11-16): the closed lines of the curve's unit tangent bundle - the wave's ledger entries, the recorded miss,
the frozen record's provenance, and one pass through the tool from a bare matrix to its thirty integers and its crossings."""
import pytest

import mtft.lines as LN
from mtft.lines import closed as C
from mtft.lines import crossings as X
from mtft.lines import involutions as IV
from mtft.lines import lorenz as LZ
from mtft.lines import periods as P
from mtft.lines import table as TB
from mtft.matter import ledger as L
from mtft.matter import puzzle as Z


def test_ledger_entries_and_the_map():
    E = {e["id"]: e for e in L.E}
    ids = [e["id"] for e in L.E if e["wave"] == "LN-01"]
    assert ids == [f"LN-{k}" for k in range(1, 9)]
    assert all(E[i]["class"] == "EXACT" and E[i]["status"] == "active" for i in ("LN-1", "LN-2", "LN-3", "LN-5", "LN-6", "LN-7"))
    assert E["LN-4"]["class"] == "NUMERICAL" and E["LN-8"]["class"] == "STRUCTURAL" and E["LN-8"]["status"] == "caveat"
    assert "not checked" in E["LN-8"]["note"] and "gravity" in E["LN-8"]["note"]
    assert "recorded miss" in E["LN-6"]["statement"] and "no mechanism" in E["LN-6"]["note"]
    assert "trace 4" in E["M05-3"]["statement"] and E["M05-3"]["status"] == "active"              # the earlier claim stands as it was
    assert Z.WAVE_DATES["LN-01"] == "2026-11-16" and Z.WAVE_REPORT["LN-01"] == "LN01_CLOSED_LINES.md"
    assert {Z.assignment()[i] for i in ids} == {"II"} == {Z.assignment()["M05-3"]}
    assert all(L.check().values()) and all(Z.check().values())


def test_the_miss_is_among_the_recorded_negatives():
    n23 = [n for n in Z.NEGATIVES if n["id"] == "N23"]
    assert len(n23) == 1 and n23[0]["piece"] == "II" and n23[0]["claims"] == ("LN-6",) and "LN01_CLOSED_LINES.md" in n23[0]["docs"]
    assert IV.conjecture_p13()["holds"] is False and IV.conjecture_p13()["negated_by_none"] == [22, 24]
    assert len({n["id"] for n in Z.NEGATIVES}) == len(Z.NEGATIVES)


def test_the_frozen_record_and_the_front_page():
    assert TB.verify_provenance() == {"lines_143.json": True}
    meta = TB.metadata()
    assert meta["level"] == 143 and "MP12" in meta["source"] and "MP14" in meta["source"]
    assert set(LN.__all__) == {"closed", "crossings", "cutting", "darts", "genus", "involutions", "lorenz", "periods", "rademacher", "table",
                               "twists", "intersection", "walls", "twist_energy"}   # genus: wave 7; twists, intersection, walls: wave 8; twist_energy: wave 9
    assert set(TB.motions()) == set(C.MOTIONS) | set(C.MORE_MOTIONS)
    assert all(TB.motions()[k] == A for k, A in {**C.MOTIONS, **C.MORE_MOTIONS}.items())
    rec = TB.record()
    rec["lines"][0]["psi"] = 99                                                    # a copy: changing it does not touch the record
    assert TB.record("lines")[0]["psi"] == 0


def test_from_a_matrix_to_its_integers_and_its_crossings():
    """The tool on a motion no record holds: U^2 V^3 = [[7, 2], [3, 1]] (trace 8, Lorenz word xxyyy, Psi = -1)."""
    A = ((7, 2), (3, 1))
    assert LZ.lorenz_word(A) == "xxyyy" and C.transfer(A)["psi"] == -1
    with pytest.raises(ValueError):
        X.Crossings(143, {"square": ((5, 3), (3, 2))})                             # the cat map squared: a proper power is refused
    lines = C.lifts(A, 143, "new")
    assert sum(l.ell for l in lines) == 168 and all(l.psi == l.ell * lines[0].psi // lines[0].ell for l in lines)
    t = C.transfer(A)
    assert not any(t["sum_cls"]) and t["sum_psi"] == 168 * t["psi"] and t["sum_windings"] == (0, 0, 0)
    for l in lines:
        v = P.thirty(l.gamma)
        assert v["vector"][0] == l.psi and v["cls"] == l.cls and len(v["vector"]) == 30
        assert v["psi"] % 6 == l.chi
    Zc = X.crossings(143, {"cat": C.MOTIONS["cat"], "new": A})
    assert Zc.check_conventions() and len(Zc.geodesics) >= 5
    assert Zc.total("cat", "new") == 168 * LZ.transfer_sum(C.MOTIONS["cat"], A)    # the transfer, on a pair no record holds
    assert all(v is None or v <= 0 for row in Zc.linking() for v in row)
