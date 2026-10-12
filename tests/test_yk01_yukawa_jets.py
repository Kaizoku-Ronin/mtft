"""YK-01: the jet–jet Yukawa theorem, the twist classification and the price of tree-level Yukawas."""
import sympy as sp
from mtft.matter import yukawa_jets as YJ, cusp_families as CF


def test_jet_representatives_wedge_to_zero():
    x, xb, y, yb = sp.symbols("x xbar y ybar")
    chi = sp.Function("chi")(x * xb + y * yb); h1, h2 = sp.Function("h1")(x, y), sp.Function("h2")(x, y)
    a = (sp.diff(chi, xb) * h1, sp.diff(chi, yb) * h1); b = (sp.diff(chi, xb) * h2, sp.diff(chi, yb) * h2)
    assert sp.simplify(a[0] * b[1] - a[1] * b[0]) == 0                         # (dbar chi h1) ^ (dbar chi h2) = 0


def test_only_the_matter_untwisting_wilson_line_allows_tree_yukawas():
    c = YJ.classify()
    assert c["solutions"] == [{"4": 0, "L": 1, "R": 1}, {"4": 1, "L": 0, "R": 0}] and c["matter_blocks_share_fermion_twist_when_higgs_untwisted"]
    assert not YJ.tree_yukawa_possible({"4": 0, "L": 0, "R": 0})               # the jet-only backgrounds since Matter-04


def test_standard_model_pattern():
    assert YJ.sm_pattern() == {"quark_lepton_fermions": "untwisted", "squarks_sleptons": "jets", "higgs_bosons": "untwisted",
                               "higgsinos": "jets", "gauge_bosons": "untwisted", "gauginos": "jets"}


def test_price_on_the_three_cusp_background():
    e = YJ.untwisted_matter_entry()
    assert e["Hom(L->4)"] == (0, 13, 13) and e["Hom(4->L)"] == (0, 15, 12) and (e["net"], e["pairs"]) == (3, 25)
    assert CF.chi((1, 0, 0, -1, -1, -1)) == -3 == e["Hom(4->L)"][0] - e["Hom(4->L)"][1] + e["Hom(4->L)"][2]
    assert CF.chi((-1, 0, 0, 1, 1, 1)) == 0 == e["Hom(L->4)"][0] - e["Hom(L->4)"][1] + e["Hom(L->4)"][2]
    s = YJ.spectrum_comparison(); assert s["jet_only_total"] == 96 and s["matter_untwisted_total"] == 896
