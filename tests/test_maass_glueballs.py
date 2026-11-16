"""mtft.maass.glueballs: H24 against the lattice.  Every expected number below is the lab's recorded one (records MP15, MP18,
MP20, 12 to 14 November 2026: needed_shifts.json and mp20_numbers.json) or is built in the test and says so."""
import math

import pytest

from mtft.falsify import Prediction
from mtft.maass import glueballs as G
from mtft.maass import spectrum as SP

# the lab's records, state by state: (lattice ratio, its error, curve ratio, needed shift of lambda, its error, sigma)
RECORD = {"0-+": (1.5494860499, 0.0163089494, 1.5590607204, -0.0023936676, 0.0040772374, 0.5870807643),
          "0++*": (1.7195301028, 0.0160454095, 1.7920300867, -0.0181249960, 0.0040113524, 4.5184252700),
          "0-+*": (2.1409691630, 0.0403980032, 2.2345494633, -0.0233950751, 0.0100995008, 2.3164585593),
          "A1+-": (2.7371512482, 0.0839468618, 1.9219314122, 0.2038049590, 0.0209867154, -9.7111412979),
          "A1--": (2.8722466960, 0.1362518857, 1.8565188644, 0.2539319579, 0.0340629714, -7.4547799961)}


def test_the_lattice_table_and_its_ratios():
    assert {k: v[:2] for k, v in G.LATTICE.items()} == {"0++": (3.405, 0.021), "0++*": (5.855, 0.041), "0-+": (5.276, 0.045),
                                                         "0-+*": (7.29, 0.13), "A1+-": (9.32, 0.28), "A1--": (9.78, 0.46)}      # MP15's targets
    assert "2007.06422" in G.SOURCE and "not re-read against the paper" in G.TRANSCRIPTION
    for state, rec in RECORD.items():
        r, e = G.lattice_ratio(state)
        assert r == pytest.approx(rec[0], abs=1e-9) and e == pytest.approx(rec[1], abs=1e-9)
    C = G.covariance(["0-+", "0++*", "0-+*"])                                     # the three share the error of the 0++
    cor = lambda i, j: C[i][j] / math.sqrt(C[i][i] * C[j][j])
    assert (round(cor(0, 1), 2), round(cor(0, 2), 2), round(cor(1, 2), 2)) == (0.39, 0.19, 0.22)      # MP20's figures
    assert G.covariance(["0-+", "0++*"], correlated=False)[0][1] == 0.0
    assert set(G.DICTIONARY) == {"D1", "D2", "D3", "D4", "D4'", "D5a", "D5b"}


def test_needed_shifts_reproduce_the_record():
    need = G.needed_shifts()
    assert list(need) == list(RECORD)
    for state, rec in RECORD.items():
        d = need[state]
        assert d["curve_ratio"] == pytest.approx(rec[2], abs=1e-8) and d["shift"] == pytest.approx(rec[3], abs=1e-9)
        assert d["shift_error"] == pytest.approx(rec[4], abs=1e-9) and d["sigma"] == pytest.approx(rec[5], abs=1e-6)
        assert d["curve_ratio"] == pytest.approx(1 + 4 * SP.form(d["form"])["R"] ** 2, abs=1e-12)       # D1: lambda over the edge
    assert [round(need[k]["percent"], 1) for k in need] == [-0.6, -4.0, -4.2, 42.4, 54.7]               # MP20's table


def test_reading_b_misses_all_five():
    rb = G.reading_b()
    assert rb["forms"] == {"0++": "u.e.mm.1", "0++*": "u.e.pp.1", "0-+": "u.o.mm.1", "0-+*": "u.o.pp.1"}
    assert rb["N0"]["holds"] is False and rb["N0"]["value"] == ("u.o.mm.1", 1)    # the lowest cusp form is odd
    assert [round(rb[k]["value"], 3) for k in ("N1", "N2", "N3", "N4")] == [0.870, 2.236, 1.247, 2.152]      # MP15's result
    assert [rb[k]["window"] for k in ("N1", "N2", "N3")] == [(1.472, 1.627), (1.634, 1.806), (2.020, 2.262)]
    assert not any(rb[k]["holds"] for k in ("N0", "N1", "N2", "N3", "N4")) and rb["all_miss"] is True
    assert rb["N4"]["form"] == "u.o.mp.1" and SP.form("u.o.mp.1")["R"] == pytest.approx(0.8451232343, abs=1e-9)
    assert [f["id"] for f in G.channel(1, 1)][:2] == ["u.o.mm.1", "u.o.pp.1"] and G.channel(0, -1)[0]["eps11"] * G.channel(0, -1)[0]["eps13"] == -1


def test_reading_a_and_the_eta_twist():
    ra = G.reading_a_checks()
    assert ra["S2"]["holds"] and ra["S2"]["R"] == pytest.approx(0.3738518157, abs=1e-9) and round(ra["S2"]["sigma"], 1) == 0.6
    assert ra["MP18-P1"]["holds"] is False and ra["MP18-P1"]["twisted_forms_below"] == ["t.o.mm.1", "t.e.pp.1", "t.o.pm.1", "t.e.pm.1"]
    assert [round(G.curve_ratio(i), 4) for i in ra["MP18-P1"]["twisted_forms_below"]] == [1.8565, 1.9219, 1.9954, 2.3695]      # MP18's table
    assert ra["MP19-P4"]["holds"] and ra["MP19-P4"]["even_forms_with_fricke_plus"] == []
    rows = G.comparison()
    assert all(isinstance(p, Prediction) for p in rows) and len(rows) == 5
    assert [p.status for p in rows] == ["PASS", "FAIL", "TENSION", "FAIL", "FAIL"]
    assert [round(p.sigma, 2) for p in rows] == [0.59, 4.52, 2.32, 9.71, 7.45]
    assert all("registered miss" in p.multiplicity for p in rows[3:]) and rows[0].observed == pytest.approx(1.5495, abs=5e-5)
    s = G.summary()
    assert s["reading_b_all_miss"] and s["S2_holds"] and not s["eta_twist_bound_holds"] and "fails" in s["status"] and s["seals"] == {"GB-S1": True}


def test_the_fits_reproduce_the_record():
    p1 = G.self_interaction_fit()                                                 # MP20's P1
    assert p1["kappa"] == pytest.approx(-0.005881013130578237, abs=1e-9) and p1["chi2"] == pytest.approx(10.203614993506449, abs=1e-5)
    assert p1["p"] == pytest.approx(0.006085736669066765, abs=1e-8) and p1["passes"] is False and p1["dof"] == 2
    assert p1["chi2_without"] == pytest.approx(24.362908864937044, abs=1e-5) and p1["needed_ratio"] == pytest.approx(0.132064, abs=1e-5)
    assert p1["actual_ratio"] == pytest.approx(0.811807, abs=1e-6) and [round(x, 3) for x in p1["fitted"]] == [-0.011, -0.014, -0.013]
    one = G.one_coupling_fit()                                                    # "one coupling for all five": independent errors
    assert one["kappa"] == pytest.approx(-0.003221308117902022, abs=1e-9) and one["chi2"] == pytest.approx(169.60315797918201, abs=1e-4)
    assert one["chi2_without"] == pytest.approx(176.0, abs=0.05) and one["dof"] == 4 and one["signs_needed"] == [-1, -1, -1, 1, 1]
    full = G.one_coupling_fit(correlated=True)                                    # found here: with the 0++ error shared by all five
    assert 170 < full["chi2"] < 185 and full["chi2"] > one["chi2"] and full["p"] < 1e-30
    tw = G.twist_energy()                                                         # "a constant twist energy"
    assert tw["mu"] == pytest.approx(0.21759746061293436, abs=1e-9) and tw["error"] == pytest.approx(0.017867681390085623, abs=1e-10)
    assert tw["chi2"] == pytest.approx(1.5697330428463627, abs=1e-6) and tw["p"] == pytest.approx(0.21024578659183224, abs=1e-7)
    assert [round(v, 2) for v in tw["states"].values()] == [2.73, 2.79, 2.87, 3.24] and round(tw["state_error"], 2) == 0.07
    assert round(tw["mu_over_edge"], 2) == 0.87
    ms = G.minus_sector_fit()                                                     # MP20's P2
    assert ms["kappa"] == pytest.approx(0.09297735694987283, abs=1e-8) and ms["kappa_error"] == pytest.approx(0.009574297613182117, abs=1e-9)
    assert ms["lowest_odd"] == "t.o.mm.1" and ms["landed"]["t.o.mm.1"] == pytest.approx(2.5644385128710825, abs=1e-7)
    assert ms["sigma"] == pytest.approx(-1.9919356158357262, abs=1e-6) and ms["sigma_lattice_error_alone"] == pytest.approx(-2.259111363244011, abs=1e-6)
    assert ms["passes_by_the_letter"] is True and ms["as_a_fit"]["chi2"] == pytest.approx(3.9678074976348534, abs=1e-6)
    assert ms["as_a_fit"]["p"] == pytest.approx(0.04637812256973021, abs=1e-7) and ms["as_a_fit"]["p"] < 0.05      # P1's rule would fail it
    assert G.chi2_sf(0.0, 3) == 1.0 and G.chi2_sf(3.841458820694124, 1) == pytest.approx(0.05, abs=1e-12)      # the 5 % point of one degree


def test_participation_ratios_by_the_package_own_quadrature():
    mine = G.participation(recompute=True)                                        # about 12 s: mtft.maass.ford on the Ford domain
    assert list(mine) == list(G.PARTICIPATION)
    assert max(abs(mine[k] - G.PARTICIPATION[k]) for k in mine) < 1e-8            # the lab's record, to its tenth digit
    assert all(1.5 <= v <= 4.5 for v in mine.values())                            # MP20's F1
    p1 = G.self_interaction_fit(mine)
    assert p1["chi2"] == pytest.approx(10.2036, abs=1e-4) and not p1["passes"]
    assert G.one_coupling_fit(mine)["chi2"] == pytest.approx(169.603, abs=1e-3)


def test_the_two_open_forecasts_are_sealed(monkeypatch):
    seal = G.SEALS[0]
    assert G.verify_seals() == {"GB-S1": True} and len(seal["sha256"]) == 64 and seal["outcome"] == "pending" and seal["evaluated"] is None
    assert seal["sealed"] == "2026-11-16" and seal["registered"].startswith("2026-11-14 08:21") and "before any lookup" in G.TWIST_RULE
    tw, p = G.twist_energy(), G.TWIST_PARAMS                                      # the sealed numbers are the fit's, to four decimals
    assert abs(p["second_odd"] - tw["second_odd"]) < 5e-5 and abs(p["second_even"] - tw["second_even"]) < 5e-5
    assert abs(p["forecast_error"] - tw["state_error"]) < 5e-5 and (p["unit"], p["unit_error"]) == G.LATTICE["0++"][:2]
    blank = G.evaluate_open()
    assert blank["seal"] == "GB-S1" and blank["result"]["overall"] == "NO DATA"
    # the rule's arithmetic, on inputs built here from the forecast itself (not lattice values)
    unit, err = p["unit"], 0.30
    def at(n_sigma, fc):
        r = fc
        for _ in range(60):                                                       # solve r so that the rule returns n_sigma
            e = r * math.hypot(err / (r * unit), p["unit_error"] / unit)
            r = fc + n_sigma * math.hypot(p["forecast_error"], e)
        return (r * unit, err)
    assert G.twist_test(at(0.0, p["second_odd"]), at(1.5, p["second_even"]))["overall"] == "PASS"
    assert G.twist_test(at(2.5, p["second_odd"]), at(0.5, p["second_even"]))["overall"] == "TENSION"
    assert G.twist_test(at(1.0, p["second_odd"]), at(3.5, p["second_even"]))["overall"] == "FAIL"
    one = G.twist_test(at(1.0, p["second_odd"]), None)
    assert one["overall"] == "PASS (one channel)" and one["A1-- first excited (second odd)"]["sigma"] == pytest.approx(1.0, abs=1e-9)
    # a seal that no longer matches is refused: changed numbers, or a changed hash
    assert G.seal_hash(dict(seal, params=dict(seal["params"], thresholds=[2, 4]))) != seal["sha256"]
    monkeypatch.setitem(G.TWIST_PARAMS, "second_odd", 2.95)
    assert G.verify_seals() == {"GB-S1": False}
    with pytest.raises(RuntimeError):
        G.evaluate_open()


def test_a_new_dictionary_is_tested_forecast_first():
    rec = G.register_forecast("a made-up dictionary", "for this test only", {"0-+": (1.60, 0.02), "0++*": 1.70}, "2026-11-16 07:10 America/New_York")
    assert len(rec["sha256"]) == 64 and rec["predictions"] == {"0-+": [1.60, 0.02], "0++*": [1.70, 0.0]}
    rows = G.compare_forecast(rec)
    assert [type(p) for p in rows] == [Prediction, Prediction] and rows[0].status == "PASS" and rows[1].status == "PASS"
    assert rows[0].sigma == pytest.approx(abs(1.60 - G.lattice_ratio("0-+")[0]) / math.hypot(G.lattice_ratio("0-+")[1], 0.02))
    changed = dict(rec, predictions={"0-+": [1.55, 0.02], "0++*": [1.70, 0.0]})   # moved after sealing
    for bad in (changed, {"name": "x", "predictions": {"0-+": [1.5, 0.0]}}, None):
        with pytest.raises(ValueError):
            G.compare_forecast(bad)
    for kwargs in ({"predictions": {"0++": 1.0}}, {"predictions": {"2++": 1.4}}, {"predictions": {}}, {"registered": " "}):
        args = {"name": "n", "statement": "s", "predictions": {"0-+": 1.5}, "registered": "2026-11-16 07:10", **kwargs}
        with pytest.raises(ValueError):
            G.register_forecast(**args)


def test_the_register_of_forecasts_and_outcomes():
    reg = G.register()
    ids = [r["id"] for r in reg]
    assert len(ids) == len(set(ids)) == 23 and all(set(r) == {"id", "record", "registered", "whose", "probability", "statement", "outcome", "result"} for r in reg)
    by = {r["id"]: r for r in reg}
    rb, ra = G.reading_b(), G.reading_a_checks()
    for k in ("N0", "N1", "N2", "N3", "N4"):                                      # the register agrees with what is recomputed
        assert (by[f"MP15-{k}"]["outcome"] == "missed") == (not rb[k]["holds"])
    assert (by["MP15-S2"]["outcome"] == "passed") == ra["S2"]["holds"] and (by["MP18-P1"]["outcome"] == "missed") == (not ra["MP18-P1"]["holds"])
    assert (by["MP19-P4"]["outcome"] == "held") == ra["MP19-P4"]["holds"] and (by["MP20-P1"]["outcome"] == "failed") == (not G.self_interaction_fit()["passes"])
    assert by["MP20-F1"]["outcome"] == "held" and all(1.5 <= v <= 4.5 for v in G.PARTICIPATION.values())
    assert by["GB-S1"]["outcome"] == "open" and by["MP15-F4"]["outcome"] == "missed" and by["MP15-F4"]["probability"] == 0.9
    reg[0]["outcome"] = "changed"                                                 # a copy
    assert G.register()[0]["outcome"] == "missed"
    assert sorted({r["record"] for r in reg}) == ["MP15", "MP18", "MP19", "MP20"]
