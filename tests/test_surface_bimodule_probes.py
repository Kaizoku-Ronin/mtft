"""v0.27.3: adversarial probes from Astra's release audit must fail on the right gate."""
import numpy as np

from mtft.surface import bimodule as BM


def test_scalar_twist_probe_fails_representation_gate():
    T = BM.RealTriple.build({"I": np.eye(4)}, None, 2 * np.eye(4))
    Z = np.zeros((4, 4))
    g = T.axiom_gates(np.block([[Z, Z], [Z, Z]]))
    assert g["status"] == "FAIL" and g["representation"] == "FAIL"
    assert T.representation_gates()["twist_orthogonal"] > 0.1      # W = 2I is not orthogonal
    assert abs(T.left(np.eye(4))[4:, 4:] - np.eye(4)).max() < 1e-12  # π is now unital even so


def test_af09_alphabet_representation_gates_pass():
    from mtft.surface import frozen as FR
    d = FR.x0143()
    Jint = d["intersection_cycles"].astype(float)
    G = Jint @ d["J_true"]; G = (G + G.T) / 2
    G = -G if np.linalg.eigvalsh(G)[0] < 0 else G
    T = BM.RealTriple.build({k: d[k].astype(float) for k in ("T2", "T3")}, G, d["W13"].astype(float))
    assert T.representation_gates()["status"] == "PASS"
    r = BM.two_state_cyclic_control(3)
    assert r["K_h=Wh+hW^-1"]["status"] == "PASS" and r["K_h=Wh+hW^-1"]["representation"] == "PASS"
