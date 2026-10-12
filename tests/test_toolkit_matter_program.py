"""Toolkit gates: every preset reproduces its wave's numbers; the ledger is consistent; the command line works."""
import contextlib, io
from mtft.matter.background import PRESETS, ps_symmetric_gamma_tradeoff
from mtft.matter import ledger as L


def test_presets_reproduce_the_waves():
    m = PRESETS["matter04"].analyse()
    assert m["indices"] == {"(4,2,1)": 3, "(4bar,1,2)": 3, "(1,2,2)": 6} and m["zero_modes"]["total"] == 288 and m["d_flat"]["d_flat"]
    assert m["yukawa"].startswith("tree-level Yukawas vanish") and m["wilson_line_stability"].startswith("minimum") and m["ckm"]["admissible"]
    assert PRESETS["three-cusps"].analyse()["zero_modes"]["total"] == 96 and not PRESETS["three-cusps"].analyse()["ckm"]["admissible"]
    assert PRESETS["two-cusps-dparity"].analyse()["zero_modes"]["total"] == 216 and not PRESETS["two-cusps-dparity"].analyse()["ckm"]["admissible"]
    u = PRESETS["matter-untwisted"].analyse()
    assert u["zero_modes"]["total"] == 896 and u["yukawa"].startswith("tree-level Yukawas possible") and u["wilson_line_stability"].startswith("maximum")
    assert abs(PRESETS["matter04"].messenger_scale() - 3.537) < 1e-3
    assert PRESETS["ps-cheap"].analyse()["higgs"].startswith("global") and "gamma-localised" in PRESETS["ps-nu-optimal"].analyse()["higgs"]
    assert all(PRESETS[k].analyse()["indices"] == {"(4,2,1)": 3, "(4bar,1,2)": 3, "(1,2,2)": 6} for k in PRESETS)
    assert all(PRESETS[k].analyse()["d_flat"]["d_flat"] for k in PRESETS)


def test_gamma_flux_tradeoff():
    t = ps_symmetric_gamma_tradeoff()
    assert t["min_abs_nu_c_charge_with_global_higgs"] == 2 and not t["nu_c_charge_one_possible_with_global_higgs"]


def test_ledger():
    c = L.check(); assert all(c.values()), c
    assert {e["id"] for e in L.entries(status="retracted")} == {"AL1-4", "CP-5", "XQ-4", "TF-3", "CD-4", "LEP-4"}
    assert L.to_markdown().count("\n") == len(L.E) + 4


def test_command_line():
    from mtft.matter.__main__ import main
    for argv in (["presets"], ["background", "matter04"], ["background", "three-cusps", "--json"], ["ledger"], ["ledger", "--markdown"], ["epsilon"]):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf): assert main(argv) == 0
        assert buf.getvalue().strip()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf): assert main(["background", "no-such-preset"]) == 2
