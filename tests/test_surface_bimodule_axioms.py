"""v0.27.1: full real-spectral-triple gate set in surface.bimodule (reality gate was missing before)."""
import numpy as np

from mtft.surface import bimodule as BM


def test_two_state_cyclic_control_reproduces_addendum_table():
    r = BM.two_state_cyclic_control(3)
    assert r["[[0,Wh],[Wh^T,0]]"]["status"] == "FAIL" and r["[[0,Wh],[Wh^T,0]]"]["reality_JD_equals_DJ"] > 0.1
    assert r["[[0,Wh],[Wh,0]]"]["status"] == "FAIL" and r["[[0,Wh],[Wh,0]]"]["selfadjoint"] > 0.1
    assert r["K_h=Wh+hW^-1"]["status"] == "PASS" and r["K_h=Wh+hW^-1"]["rank_d"] == 4
    assert r["K=W+W^-1"]["status"] == "PASS" and abs(r["K=W+W^-1"]["max_one_form_size"] - 2.0) < 1e-9
    assert r["pairing"]["unit_in_radical"] and not r["pairing"]["nondegenerate_on_supplied_basis"]


def test_two_factor_involution_has_zero_calculus():
    r = BM.two_state_cyclic_control(2)
    assert r["K_h=Wh+hW^-1"]["status"] == "PASS" and r["K_h=Wh+hW^-1"]["rank_d"] == 0


def test_af09_frozen_alphabet_passes_full_gates_with_symmetric_dirac():
    from mtft.surface import frozen as FR
    d = FR.x0143()
    Jint = d["intersection_cycles"].astype(float)
    G = Jint @ d["J_true"]; G = (G + G.T) / 2
    G = -G if np.linalg.eigvalsh(G)[0] < 0 else G
    T = BM.RealTriple.build({k: d[k].astype(float) for k in ("T2", "T3", "U13")}, G, d["W13"].astype(float))
    K = BM.symmetric_dirac_block(T.twist, np.eye(26))
    Z = np.zeros((26, 26))
    g = T.axiom_gates(np.block([[Z, K], [K.T, Z]]))
    assert g["status"] == "PASS"                       # AF-09 with the W13 twist is a genuine real even triple ...
    assert g["max_one_form_size"] < 1e-9              # ... with identically vanishing calculus (involutive twist)
