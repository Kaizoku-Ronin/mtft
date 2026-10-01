"""Independent review regressions for degenerate loci and omitted representation data.

These checks do not change the catalogue's physics assumptions or finite scans.
"""
from fractions import Fraction
import pytest
import sympy as sp
from mtft import parent as P
from mtft.research import anomaly8d as A8
from mtft.surface import blowup as B


def test_identical_flat_stacks_have_all_positive_shape_loci():
    stacks = (("c", 3, (0, 0)), ("L", 2, (0, 0)), ("1", 1, (0, 0)))
    candidate = P.ParentCandidate("FLAT", 8, "U(6)", (("adjoint", 1),), stacks=stacks, higgs_blocks={"u": (0, 0)})
    gates = P.gates(candidate)
    assert gates["4_colour_slope"]["pass"] is True
    assert gates["4_colour_slope"]["witness"]["r_colour"] == "all"
    assert gates["7_light_higgs"]["pass"] is True
    assert gates["6_tachyon_free"]["witness"]["polystable_r"] == "all"


def test_identically_zero_higgs_slope_can_share_a_nonzero_block_locus():
    candidate = P.ParentCandidate("MIXED", 8, "U(6)", (("adjoint", 1),), higgs_blocks={"u": (0, 0), "d": (2, -4)})
    gates = P.gates(candidate)
    assert gates["7_light_higgs"]["pass"] is True
    r = gates["7_light_higgs"]["witness"]["slope_free_loci"]["d"]
    assert r > 0 and all(a + b * r == 0 for a, b in candidate.higgs_blocks.values())


def test_unknown_anomaly_content_is_open_rather_than_silently_omitted():
    candidate = P.ParentCandidate("CHARGE", 8, "U(1)", (("u1(1)", 1),))
    assert A8.I10([(A8.u1(1), 1)]) != 0
    verdict = P.gates(candidate)["9_anomaly"]
    assert verdict["pass"] is None
    assert verdict["status"].startswith("OPEN")
    assert verdict["witness"]["unsupported_content"] == [("u1(1)", 1)]


def test_identical_blowup_stacks_fill_the_ample_cone():
    locus = B.polystable_locus((("c", 3, (0, 0, 0)), ("L", 2, (0, 0, 0))))
    assert len(locus) == 1
    r, t = locus[0]["r"], locus[0]["eps_over_A_E"]
    assert r.free_symbols and t.free_symbols
    for rr, tt in ((Fraction(1, 2), Fraction(1, 3)), (3, Fraction(1, 2)), (Fraction(1, 4), Fraction(1, 3))):
        assert bool(locus[0]["ample"].subs({r: rr, t: tt})) == B.ample(rr, 1, tt)


def test_collinear_scaling_changes_chirality_but_not_d_flat_shape():
    for a in (1, 2, -1):
        family = B.collinear_family(a)
        t = family["d_flat_locus"]
        assert t == sp.Rational(1, 3)
        assert family["ample_iff"] == "r > 1/3"
        assert B.ample(sp.Rational(1, 2), 1, t)
        assert all(B.slope(d, sp.Rational(1, 2), 1, t) == 0 for d in family["stacks"].values())
        assert family["net_per_unit_charge"] == -3 * a
    with pytest.raises(ValueError, match="nonzero"):
        B.collinear_family(0)


def test_hypercharge_normalisation_preserves_weak_triplet_exotic():
    key = ("1", "3w", (1,))
    assert A8.sm_net_content({key: 3}, {key: Fraction(1)}) == {("1", "3w", Fraction(1)): 3}


def test_nonzero_colour_product_can_have_exact_sm_net_content():
    # Outside the released scan boxes: a counterexample to necessity of its
    # imposed colour-zero restriction, not to its finite counts or trilemma.
    stacks = (("c", 3, (1, 1)), ("L", 2, (0, 2)), ("s1", 1, (0, 0)),
              ("s2", 1, (0, -12)), ("s3", 1, (-3, -1)))
    y = (Fraction(-1, 3), Fraction(1, 2), Fraction(1), Fraction(0), Fraction(0))
    blocks = A8.b1_content_blocks(stacks, 8)
    hypercharges = {key: sum(q * yi for q, yi in zip(key[2], y)) for key in blocks}
    net = A8.sm_net_content(blocks, hypercharges)
    assert stacks[0][2][0] * stacks[0][2][1] == 1
    assert net[("1", "1", Fraction(0))] == 63
    assert {key: value for key, value in net.items() if key != ("1", "1", Fraction(0))} == A8.SM_NET(3)
