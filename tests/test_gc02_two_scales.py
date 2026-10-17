"""GC-02 (2026-10-17): the flux-split couplings matched at the curve's compactification scale M_C above the Pati–Salam scale, with
Pati–Salam running between — R rises, and the designs R = 5 and 9 are met with M_PS at the b–tau crossing."""
import math

from mtft.matter import two_scales as G
from mtft.matter import gauge_couplings as G1
from mtft.matter import ledger as L
from mtft.matter import puzzle as Z

BT = 1.1236e6


def test_ps_coefficients_and_left_right_symmetry():
    b = G.ps_b()
    assert abs(b["4"] + 32 / 3) < 1e-12 and abs(b["L"] - 1.0) < 1e-12 and b["L"] == b["R"]
    for nH in (0, 1, 2, 3):
        b = G.ps_b(n_higgs=nH); assert b["L"] == b["R"]
    b = G.ps_b(n_delta=1); assert b["R"] > b["L"]                                              # only a (4bar,1,2) scalar breaks L-R symmetry


def test_continuity_and_rise():
    a = G.ps_inverse_couplings(BT, BT)
    assert abs(a["4"] - G1._inv("3", BT)) < 1e-12 and abs(a["L"] - G1._inv("2", BT)) < 1e-12 and abs(G.two_scale_ratio(BT, BT) - G1.ratio(BT)) < 1e-12
    rs = [G.two_scale_ratio(BT * 10 ** k, BT) for k in (0, 0.5, 1, 1.5, 2)]
    assert all(x < y for x, y in zip(rs, rs[1:])) and abs(rs[0] - 3.678) < 1e-3


def test_compactification_scales():
    r = G.two_scale_report(BT)
    for R in (1 / 9, 3 / 7, 1.0, 7 / 3, 1 / 5):
        assert r["designs"][R] is None                                                          # below the SM value: only M_C = M_PS (GC-01)
    d5, d9 = r["designs"][5.0], r["designs"][9.0]
    assert 1.6e7 < d5["M_C"] < 1.8e7 and 2.6e8 < d9["M_C"] < 2.9e8
    assert 2.4e8 < min(d5["M_KK"]) and max(d9["M_KK"]) < 4.5e9 and abs(d5["y_b_over_y_tau_at_M_PS"] - 1) < 0.02
    assert 1.0e13 < G.compactification_scale(5.0, 8.9e11) < 1.2e13                              # GC-01's lowest matching scale as M_PS
    assert 6e7 < G.compactification_scale(5.0, BT, n_bidoublet_fermions=0) < 8e7                 # content sensitivity


def test_record():
    e = {x["id"]: x for x in L.E}
    assert e["GC-6"]["class"] == "EXACT" and e["GC-7"]["class"] == "NUMERICAL" and e["GC-6"]["wave"] == "GC-02"
    for cc in ("GC-4", "CC-48", "SC-2"):
        assert e[cc]["status"] == "caveat" and "GC-02" in e[cc]["note"]
    sub = Z.subs()["VI.4"][1]
    assert sub.status == "open" and "GC-7" in sub.claims and "GC02_TWO_SCALES.md" in sub.docs
    assert Z.WAVE_DATES["GC-02"] == "2026-10-17" and all(Z.check().values()) and all(L.check().values())
