"""v0.27.0: cycle frame <-> canonical homology frame, and frame-independent Hodge cross-checks."""
import numpy as np

from mtft.surface import frozen as FR, intertwiner as IT


def test_frozen_gates_include_intertwiner():
    d = FR.x0143()
    assert all(d["gates"].values())
    assert all(k in d["gates"] for k in ("intertwiner_Pi_unimodular", "intertwiner_poincare_duality"))
    assert all(IT.canonical_frame_gates().values())


def test_cross_frame_j_invariants():
    r = IT.cross_frame_hodge_check(40)
    for label, ref in (("(+,+)", IT.J_143A1), ("(-,-)", IT.J_11A3)):
        x = r[label]
        assert x["rank"] == (2, 2) and x["cross_frame_rel"] < 1e-10
        assert abs(x["j_cycle"].real - ref) / abs(ref) < 1e-10
    assert abs(r["(-,-)"]["j_cycle"].real + 4096 / 11) < 1e-8        # the (-,-) sector is 11a3
