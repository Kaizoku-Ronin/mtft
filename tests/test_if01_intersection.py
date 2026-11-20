"""IF-01: the intersection form of X₀(143) on the Hecke-block coordinates, deg φ = 4, the homology lattice and the flux-free
lattice ker φ_* (record FF-19)."""
import pytest
import sympy as sp

from mtft.lines import intersection as I


@pytest.fixture(scope="module")
def form():
    return I.intersection_form(24)


def test_the_form_reproduces_every_crossing(form):
    J = form["J"]
    assert form["geodesics"] == 716 and form["pairs checked"] == 716 ** 2
    assert form["pairs that disagree"] == 0 and form["largest entry off the Hecke blocks"] == 0
    assert J == -J.T and J.rank() == 26
    assert J[0:2, 0:2] == sp.Matrix([[0, sp.Rational(1, 4)], [sp.Rational(-1, 4), 0]])
    assert I.modular_degree(J) == 4


def test_the_lattice_and_the_flux_free_classes(form):
    assert I.homology_lattice()["index"] == 1327104 == 2 ** 14 * 3 ** 4
    ff = I.flux_free_lattice(form["J"])
    assert ff["form on Λ integral"] and ff["det on Λ"] == 1
    assert ff["Smith invariants of φ_*"] == [1, 1]
    assert (ff["rank"], ff["det"], ff["rank over Q"], ff["rank mod 2"]) == (24, 16, 24, 22)
    assert ff["largest isotropic rank"] == 12               # < 13 meridians: no handlebody of the curve carries the γ-flux (FF-19)


def test_cup_form_against_the_crossings_mod_2():
    r = I.cup_against_crossings(24)
    assert r["geodesics"] == 716 and r["pairs with odd crossings"] > 100000
    assert r["pairs that disagree mod 2"] == 0


def test_pullbacks_vanish_on_flux_free_classes():
    assert I.pullbacks_on_flux_free() == {"kernel rank in Z^29": 27, "pullback values that are not 0": 0}


def test_too_few_geodesics_raise():
    with pytest.raises(ArithmeticError):
        I.intersection_form(5)
