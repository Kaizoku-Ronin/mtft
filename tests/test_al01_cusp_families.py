"""AL-01: Atkin–Lehner-symmetric backgrounds and the one-family-per-cusp model."""
from mtft.matter import cusp_families as CF


def test_full_atkin_lehner_symmetry_forbids_three_families():
    ob = CF.al_obstruction(box=4)
    assert ob["three_family_configurations"] == []
    assert len(ob["per_homomorphism_hits"]) == 16


def test_single_point_reproduces_matter04():
    a = CF.analyse("one point")
    assert a["fermion_total"] == 288 and a["yukawa_tree_level"].startswith("single point")


def test_two_cusp_left_right_parity_kills_the_yukawa():
    a = CF.analyse("two cusps, D-parity")
    assert [a["blocks"][r]["net"] for r in ("(4,2,1)", "(4bar,1,2)", "(1,2,2)")] == [3, 3, 6]
    assert a["fermion_total"] == 216 and a["yukawa_tree_level"].startswith("zero")


def test_one_family_per_cusp():
    a = CF.analyse("three cusps, one family each")
    b = a["blocks"]
    assert [b[r]["net"] for r in ("(4,2,1)", "(4bar,1,2)", "(1,2,2)")] == [3, 3, 6]
    assert b["(4,2,1)"]["pairs_per_entry"] == 0 == b["(4bar,1,2)"]["pairs_per_entry"]          # no vector-like quarks or leptons
    assert b["(1,2,2)"]["pairs_per_entry"] == 3 and a["fermion_total"] == 96
    assert a["yukawa_tree_level"] == "diagonal in the cusp basis"
    amp = CF.three_cusp_ampleness()
    assert all(v == 0 for v in amp["slopes"].values())
    import sympy as sp
    AX, AE = sp.symbols("A_X A_E", positive=True)
    assert sp.simplify(amp["C_x_O"] - (AX - AE)) == 0 and sp.simplify(amp["fibre"] - 2 * AE / 3) == 0
    assert CF.dot(CF.K(3), CF.K(3)) == -3 and CF.net((1, 0, 0, -1, -1, -1)) == -3
