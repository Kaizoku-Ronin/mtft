"""HP-01: the light Higgs of the ratio-3 background — localisation, type, the missing coupling, the gamma-neutral alternative."""
from fractions import Fraction
from mtft.matter import higgs_profile as H


def test_localisation_and_graph():
    assert H.localisation(1) == {"modes": 4, "torsion_points": 1, "at_cusps": 4, "elsewhere": 0}
    assert H.localisation(3)["modes"] == 36 and H.localisation(3)["elsewhere"] == 32
    assert H.phi_taylor(7) == [Fraction(1), Fraction(0), Fraction(-1, 3), Fraction(-1, 2), Fraction(-1, 5), Fraction(0), Fraction(-2, 7)]


def test_a_type_only():
    for c in (1, 3):
        for e in (6, -6, 0):
            assert H.chi_entry(c, e) == -4 * c * c - Fraction(e * (e - 1), 2) < 0


def test_no_tree_yukawa_in_ratio_three():
    assert not H.superpotential_allowed(("A", "A", "A")) and H.superpotential_allowed(("A", "A", "Phi"))
    assert H.coupling_higgs_charge(-3) == 6
