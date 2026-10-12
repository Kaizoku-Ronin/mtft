"""GF-01: the eight-dot engine reproduces the programme, the sealed rules give their verdicts, hypercharge by two routes, the torus geometry,
the three parent columns (PV-02 recovered; SO(16) and E8 closed by parity), and the clock of twelve."""
import contextlib
import io
import json
import math
from fractions import Fraction as Fr

import pytest

from mtft.matter import forces as F
from mtft.matter import ledger as L
from mtft.matter import prereg as P
from mtft.matter import puzzle as Z
from mtft.matter.__main__ import main as matter_main


# ------------------------------------------------------------------------------------------------------------------------------ 1. the engine
def test_engine_reproduces_the_programme():
    r = F.reproduce()
    assert r["ps01_hypercharges"] and r["ps01_branching"] and r["ds01_touch"] and r["yukawa_triangles_close"]
    assert r["lambda2_blocks"] == 28 and r["lambda2_block_count"] == 12 and r["groupings_of_eight"] == 22
    assert r["PS"]["total"] == 24 and r["PS"]["lines_between"] == 40 and r["PS"]["weyl_order"] == 96 and r["PS"]["check_64"]
    assert r["SM"]["total"] == 16 and r["SM"]["lines_between"] == 48 and r["SM"]["weyl_order"] == 12 and r["SM"]["nonabelian_dim"] == 11
    assert F.carriers(F.GROUPINGS["U8"])["total"] == 64 and F.particle_lines(F.GROUPINGS["U8"])["lines_between"] == 0
    assert [math.comb(8, k) for k in range(9)] == [sum(b["dim"] for b in F.sheet_blocks(k)) for k in range(9)]


def test_fields_and_yukawa_triangles():
    T = F.sm_table()
    assert T["Q"]["Y"] == Fr(1, 6) and T["u^c"]["Y"] == Fr(-2, 3) and T["d^c"]["Y"] == Fr(1, 3) and T["L"]["Y"] == Fr(-1, 2)
    assert T["e^c"]["Y"] == 1 and T["nu^c"]["Y"] == 0 and T["H_u"]["Y"] == Fr(1, 2) and T["H_d"]["Y"] == Fr(-1, 2)
    assert T["Q"]["colour"] == "3" and T["u^c"]["colour"] == "3bar" and T["Q"]["triality"] == 1 and T["d^c"]["triality"] == 2
    tri = F.yukawa_triangles()
    assert set(tri) == {"Q u^c H_u", "Q d^c H_d", "L e^c H_d", "L nu^c H_u"} and all(t["closes"] and t["sum_Y"] == 0 for t in tri.values())
    assert F.FIELDS["H_u"] == ("L", "A2") and F.FIELDS["H_d"] == ("L", "B")                       # the lines that close the triangles (DS-01)


def test_menus_one_block_per_conjugate_pair():
    assert [len(F.menu(p)) for p in ("U8", "SO16", "E8")] == [10, 22, 46]
    assert [sum(b["dim"] for b in F.menu(p)) for p in ("U8", "SO16", "E8")] == [24, 52, 116]      # half of 48, + 28, + 64
    spin = [b for b in F.menu("E8") if b["kind"].startswith("sheet")]
    assert len(spin) == 24 and sum(b["dim"] for b in spin) == 64
    assert all(q.denominator == (1 if st == "L" else 2) for b in spin for q, st in zip(b["q"], F.STACKS))   # halves on the odd-rank stacks
    e = F.eisenstein()
    assert e["directions"] == [30.0 * k for k in range(12)] and e["g2_roots"] == 12 and e["g2_weyl_order"] == 12 and e["sm_carriers"] == 12
    assert abs(e["root_over_weight"] - math.sqrt(3)) < 1e-12


# ------------------------------------------------------------------------------------------------------------------------------ 2.-3. the sealed rules
def test_touch_scan_and_seal_pr_s3():
    u8 = F.touch_scan("U8")
    assert u8["verdict"] == "FAIL" and [r["block"] for r in u8["singlets"]] == ["A1<-A2"]
    assert u8["singlets"][0]["reaches"] == ["L", "u^c", "e^c", "nu^c", "H_u"]                        # DS-1: not Q, d^c, H_d
    so = F.touch_scan("SO16")
    assert so["verdict"] == "PASS" and {r["block"] for r in so["singlets"] if r["down_sector"]} == {"B+A2", "A1+B", "L+L"}
    e8 = F.touch_scan("E8")
    assert e8["verdict"] == "PASS" and sum(r["down_sector"] for r in e8["singlets"]) == 9
    assert P.evaluate("PR-S3", "U8")["result"]["verdict"] == "FAIL" and P.evaluate("PR-S3", "E8")["result"]["verdict"] == "PASS"


def test_z6_and_seal_pr_s4():
    t = F.z6_test()
    assert t["verdict"] == "PASS" and all(r["trivial"] for r in t["rows"])
    for f, r in F.sm_table().items():
        ct = F.charge_triality(r)
        assert ct["charges_mod_1"] == [ct["two_thirds_triality"]], f                              # Q = (2/3) t mod 1 on every field
    bad = F.z6_test([F.line("c", "A1")], {**F.Y_STACK, "A1": Fr(-1, 3)})                           # a (3, 1)_{1/2} block: Z6 broken
    assert bad["verdict"] == "FAIL"
    assert P.evaluate("PR-S5")["result"]["verdict"] == "PASS"
    with pytest.raises(RuntimeError):
        P.evaluate("PR-S4")                                                                              # superseded by CC-64


# ------------------------------------------------------------------------------------------------------------------------------ 4. hypercharge
def test_hypercharge_two_routes():
    sp_ = F.anomaly_free_space()
    assert sp_["dimension"] == 2 and len(sp_["solutions"]) == 1 and set(sp_["free_symbols"]) == {"y_c", "y_B"}
    a = F.hypercharge_route_a()
    assert a["neutral_higgs_solutions"] == [-1, 1] and a["z6_passes_for_a_in"] == [-5, -3, -1, 1, 3, 5] and a["standard"] and a["swap_is_relabelling"]
    assert a["higgs_Y"] == "a/2"
    b = F.hypercharge_route_b()
    assert b["equals_standard"] and b["two_thirds"]
    assert F.hypercharge_e2()["agree"]


# ------------------------------------------------------------------------------------------------------------------------------ 5. the torus
def test_torus_geometry():
    ew = F.electroweak_torus()
    assert ew["pythagoras"] and 28.5 < ew["theta_W_deg"] < 29.0 and ew["pati_salam_equal_couplings"]["sin2"] == Fr(3, 8)
    ps = F.pati_salam_torus(1.1236e6)
    assert ps["pythagoras_Y"] and 3.6 < ps["R"] < 3.75 and all(v.real == 0 and v.imag > 0 for v in ps["tau"].values())
    sh = F.coupling_shapes()
    assert sh["two_point"] == sorted({1 / 9, 3 / 7, 1.0, 7 / 3, 9.0}) and sh["single_point"] == [0.2, 1.0, 5.0] and not sh["reaches_data"]


# ------------------------------------------------------------------------------------------------------------------------------ 6. the columns
def test_u8_column_recovers_pv02():
    s = F.exact_sm_solutions("U8", box=4)
    assert s["consistent"] and s["rank"] == 3 and s["equations"] == 5 and s["count_in_box"] == 18 and s["min_pairs"] == 6
    assert F.pv02_conditions_recovered(box=4)["conditions_hold"]
    best = s["solutions"][0]["phi"]
    assert best["c"] - best["L"] == 3 and best["B"] - best["c"] == 3 and best["A1"] + best["A2"] - 2 * best["c"] == 3


def test_so16_and_e8_columns_closed_by_parity():
    for parent in ("SO16", "E8"):
        s = F.exact_sm_solutions(parent, box=3)
        assert not s["consistent"] and s["count_in_box"] == 0
    par = F.so16_parity()
    assert par["even"] and par["Q_row_no_shift"] == ["2", "0", "0", "0", "0"] and len(par["exotic_doublets_with_shift"]) == 1
    assert par["lepton_rows_cancel"] and par["e8"]["exotic_row_is_minus_Q_row"]                       # the review's certificate; E8's three Q come with three exotics
    e8rows = F.chirality_system("E8")["rows"]
    phi6 = (6, 0, 0, 0, 0)
    assert sum(c * p for c, p in zip(e8rows[("3", 2, Fr(1, 6))], phi6)) == 3 and sum(c * p for c, p in zip(e8rows[("3", 2, Fr(-5, 6))], phi6)) == -3
    rows = F.chirality_system("SO16")["rows"]
    assert rows[("1", 2, Fr(-1, 2))] == rows[("1", 1, Fr(-1))]                                      # L and e^c counts coincide in SO(16)
    for parent in ("SO16", "E8"):
        sc = F.shift_scan(parent, denominator=48, bound=2)
        assert sc["inconsistent"] == 193 and sc["consistent_at"] == []
    e8 = F.exact_sm_solutions("E8", box=1)
    assert any("-5, 6" in r for r in e8["exotic_rows"])                                            # the exotic (3, 2)_{-5/6} sheet


def test_parent_columns_summary_and_nonsusy():
    c = F.parent_columns(box=3)
    assert c["U8"]["consistent"] and not c["SO16"]["consistent"] and not c["E8"]["consistent"]
    assert c["U8"]["min_pairs"] == 6 and c["SO16"]["best"] is None
    n = F.nonsusy_column()
    assert n["box3"] == 60 and n["box4"] == 156 and not n["supersymmetric"] and len(n["exemplars"]) == 3


# ------------------------------------------------------------------------------------------------------------------------------ 7. the clock
def test_modular_clock_and_centre():
    ck = F.modular_clock()
    assert ck["SL2Z_abelianisation"] == [12] and ck["PSL2Z_abelianisation"] == [6] and ck["subgroup_orders_of_Z12"] == [1, 2, 3, 4, 6, 12]
    assert ck["orders"]["S"] == 4 and ck["orders"]["ST"] == 6 and ck["orders"]["S^2 = -1"] and ck["orders"]["(ST)^3 = -1"]
    assert ck["level_143"] == {"g_minus_1": 12, "am_11_13": 12, "chi": -24, "spin_twist": (12, 0)}
    fm = F.five_means()
    assert (fm["HM"], fm["GM^2"], fm["AM"], fm["QM^2"], fm["CHM"]) == (Fr(143, 12), 143, 12, 145, Fr(145, 12)) and fm["HM*AM = GM^2"] and fm["AM*CHM = QM^2"]
    assert F.sm_centre()["z6_trivial_on_all_fields"] and "open" in F.sm_centre()["identification"]


def test_systoles_and_the_two_loops():
    rows = {r["N"]: r for r in F.systole_table()}
    assert [N for N, r in rows.items() if r["twelve"]] == [118, 141, 142, 143] and [N for N, r in rows.items() if r["golden"]] == [145, 155]
    assert rows[143]["disc"] == 12 and abs(rows[143]["length"] - 2 * math.log(2 + math.sqrt(3))) < 1e-12 and not rows[143]["golden_closes"]
    assert all(not r["golden_closes"] for r in rows.values() if r["even_level"])                    # even levels have even traces only
    lg = F.legendre_reading()
    assert (lg["(5/11)"], lg["(5/13)"], lg["(3/11)"], lg["(3/13)"]) == (1, -1, 1, 1)
    from mtft.matter import one_loop_adjoints as OL
    assert OL.systole_X0_143()["trace"] == 4 and OL.systole_X0_143()["trace"] == F.systole_trace(143)
    g, t = F.loop(F.GOLDEN), F.loop(F.TWELVE)
    assert g["trace"] == 3 and g["disc"] == 5 and abs(g["length"] - 4 * math.log((1 + 5 ** 0.5) / 2)) < 1e-12 and g["unit_circle"]["is_i"]
    assert g["crossings"]["Re tau = 0 (theta = 0)"]["form"] == (1, 0, 1) and g["crossings"]["Re tau = 1/2 (theta = pi)"] == \
        {"tau": "1/2 + i sqrt(5/4)", "form": (2, -2, 3), "disc": -20, "class": "theta_pi", "in_tile": True}
    assert t["trace"] == 4 and t["disc"] == 12 and t["unit_circle"]["is_rho"]
    assert {k: v["disc"] for k, v in t["crossings"].items()} == {"Re tau = 0 (theta = 0)": -8, "Re tau = 1/2 (theta = pi)": -11, "Re tau = -1/2 (theta = pi)": -3}
    assert t["crossings"]["Re tau = 0 (theta = 0)"]["class"] == "theta_0" and t["crossings"]["Re tau = 1/2 (theta = pi)"]["class"] == "theta_pi"
    arcs = F.loop_arcs()
    assert arcs["golden_check"] and arcs["twelve_check"] and abs(arcs["golden_arc"] - math.log((1 + 5 ** 0.5) / 2)) < 1e-12


def test_reversibility_icosahedron_polygons():
    r = F.reversibility()
    assert r["golden"]["S_reverses"] and r["golden"]["unit_norm"] == -1
    assert r["twelve"]["conjugators_up_to_8"] == [] and not r["twelve"]["minus_one_norm_mod_3"] and r["twelve"]["unit_norm"] == 1
    ic = F.icosahedron()
    assert (ic["index"], ic["cusps"], ic["order_2"], ic["order_3"], ic["genus"], ic["euler"]) == (60, 12, 30, 20, 0, 2)
    u = F.polygon_units()
    assert u["pentagon"]["is_phi"] and u["octagon"]["is_silver"] and str(u["dodecagon"]["ratios"][4]) == "sqrt(3) + 2" and str(u["dodecagon"]["tan75"]) == "sqrt(3) + 2"
    assert not u["gemini_claim_sin15_eq_(phi-1)/sqrt2"] and not u["phi_in_Q_sqrt2_sqrt3"]


# ------------------------------------------------------------------------------------------------------------------------------ the record
def test_trial_x17_logged_after_looking():
    t = {e["id"]: e for e in P.trials()}
    assert t["X17"]["selects"] == [118, 141, 142, 143] and t["X17"]["sealed"] is None and t["X17"]["includes_143"] and not t["X17"]["singles_out_143"]


def test_command_line():
    def run(*argv):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = matter_main(["forces", *argv])
        return code, buf.getvalue()
    code, out = run("--box", "2"); assert code == 0 and "touch (PR-S3) U8   : FAIL" in out and "2 phi_c = 3" in out and "SL(2,Z)^ab = Z/12" in out
    code, out = run("--touch", "SO16"); d = json.loads(out); assert code == 0 and d["verdict"] == "PASS"
    code, out = run("--column", "U8", "--box", "2"); d = json.loads(out); assert code == 0 and d["consistent"] and d["rank"] == 3
    code, out = run("--json", "--box", "2"); d = json.loads(out); assert code == 0 and d["columns"]["SO16"]["consistent"] is False


def test_the_record():
    gf = [e for e in L.E if e["wave"] == "GF-01"]
    assert [e["id"] for e in gf] == [f"GF-{i}" for i in range(1, 14)] and all(e["module"] == "mtft.matter.forces" for e in gf)
    assert all(e["test"] == "test_gf01_forces.py" for e in gf)
    E = {e["id"]: e for e in L.E}
    assert all(E[i]["status"] == "superseded" and "CC-64" in E[i]["note"] for i in ("GF-2", "GF-3", "GF-7", "GF-8"))
    assert all(E[i]["status"] == "active" for i in ("GF-10", "GF-11", "GF-12", "GF-13")) and "PR-S5" in E["CC-64"]["statement"]
    a = Z.assignment()
    assert a["GF-1"] == "V" and a["GF-2"] == "IX" and a["GF-13"] == "IX" and a["GF-5"] == "VI" and a["GF-8"] == "II" and a["GF-12"] == "II" and a["GF-9"] == "IV"
    subs = Z.subs()
    assert subs["IX.10"][1].status == "answered" and subs["V.7"][1].status == "answered" and subs["V.8"][1].status == "answered" and subs["I.1"][1].status == "open"
    assert "GF-11" in subs["I.1"][1].claims and "GF01_FORCES_FIRST.md" in subs["I.1"][1].docs
    assert {s["id"] for s in P.SEALS} >= {"PR-S3", "PR-S4", "PR-S5"}
    arcs = F.loop_arcs()
    assert "twice" in arcs["twelve_points"][1] and "twice" in arcs["twelve_points"][3]
