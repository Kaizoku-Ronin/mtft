"""TA-01: the two-point design's mixing test."""
from mtft.matter import two_point_audit as T


def test_fibre_geometry():
    g = T.fibre_geometry(201); assert 0.62 < g["deep_hole"] < 0.63 and g["torsion"][3][-1] == 0.5983 and g["torsion"][2] == [0.4943, 0.5629]


def test_constraints_and_mixing():
    b = T.blowups(0.55); assert b["ample_on_fibre"] and b["d_flat"] and abs(b["r2"] - 0.4184) < 1e-3 and abs(b["r1"] - 0.2676) < 1e-3
    m = T.mixing(0.55); assert 0.33 < m["top_tail_at_jets"] < 0.35 and 0.03 < m["jet_tail_at_top"] < 0.05
    m3 = T.mixing(1 / 3); assert 0.15 < m3["top_tail_at_jets"] < 0.17


def test_required_separation_exceeds_fibre():
    assert T.required_separation(1 / 3) > 0.9 and T.required_separation(0.55) > 1.2 > T.fibre_geometry(101)["deep_hole"]
