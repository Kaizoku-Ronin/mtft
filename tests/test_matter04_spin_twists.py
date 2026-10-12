"""Matter-04: spin structures on 143a1 and the Scherk–Schwarz twist — adjoint and charged fermion zero modes, flux minimum, running."""
from fractions import Fraction
import sympy as sp
from mtft.matter import spin_twists as ST
from mtft.matter.collinear import cohomology


def test_two_torsion_of_143a1_is_an_S3_cubic_over_Q_sqrt_minus_11():
    d = ST.two_division_polynomial()
    assert d["polynomial"] == 4 * sp.Symbol("x") ** 3 - 4 * sp.Symbol("x") ** 2 - 4 * sp.Symbol("x") - 7
    assert d["Delta"] == -1859 and d["discriminant"] == -29744 and d["disc_equals_16_Delta"] and d["irreducible_over_Q"]
    assert d["galois_group"] == "S3" and d["quadratic_resolvent"] == "Q(sqrt(-11))"
    a1, a2, a3, a4, a6 = ST.A_INVARIANTS
    for x, y in d["two_torsion_points"]:                                  # each root is a genuine 2-torsion point (y = -1/2)
        assert abs(y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)) < 1e-9 and abs(2 * y + a1 * x + a3) < 1e-12
    assert [s["parity"] for s in ST.elliptic_spin_structures()] == ["odd", "even", "even", "even"]


def test_adjoint_fermion_zero_modes_by_twist():
    assert ST.adjoint_fermion_zero_modes(True, True)["h"] == (1, 14, 13)                 # BHV: gaugino + 27 adjoint chirals
    assert (cohomology(0).h0, cohomology(0).h1, cohomology(0).phi_dim) == (1, 14, 13)   # agrees with Matter-02's q = 0 row
    assert ST.adjoint_fermion_zero_modes(False, True)["h"] == (0, 12, 12)                # a curve twist keeps g - 1
    for zeta in (True, False):
        assert ST.adjoint_fermion_zero_modes(zeta, False)["total"] == 0                   # any even spin structure on 143a1 removes all
    assert ST.bosonic_adjoint_zero_modes() == {"vector": 1, "wilson_line_scalars": 14, "phi_scalars": 13, "complex_scalars": 27}


def test_charged_fermions_collapse_to_exceptional_jets():
    for q in range(-5, 6):
        jets = ST.charged_fermion_block(q, eta_trivial=False)
        if q: assert jets["chiral"] == cohomology(q).h1 - cohomology(q).h2 and jets["phi"] == 0
        else: assert jets["chiral"] == 0
        assert ST.charged_fermion_block(q, True)["chiral"] - ST.charged_fermion_block(-q, True)["chiral"] == ST.charged_fermion_block(q, False)["chiral"] - ST.charged_fermion_block(-q, False)["chiral"]   # net chirality unchanged
    assert ST.fermion_census(True)["totals"] == {"Qminus": 279, "Qplus": 279}
    assert ST.fermion_census(False)["totals"] == {"Qminus": 153, "Qplus": 153}
    assert ST.family_chirality((-1, 0, 3)) == 3


def test_flux_minimum_under_scherk_schwarz():
    m = ST.minimal_flux_assignment()
    assert m["minimum"] == 288 and m["matter02"] == 1152
    best = [r for r in m["scan"] if r["fermion_zero_modes"] == 288]
    assert {r["phi_A1"] for r in best} == {0, 3} and all(r["exact_SM"] and r["d_flat"] and r["coinciding_classes"] for r in best)


def test_one_loop_sign_and_running():
    VB0, VF0, V0 = ST.ss_wilson_potential(0.0); VBpi, VFpi, Vpi = ST.ss_wilson_potential(3.141592653589793)
    assert V0 < Vpi and ST.ss_mass_squared_at_symmetric_point() > 0
    assert abs(ST.ss_wilson_potential(0.7)[2] - (-2 * sum(__import__("math").cos(n * 0.7) / n ** 5 for n in range(1, 401, 2)))) < 1e-12   # only odd windings
    assert ST.qcd_b0(6, 27, True) == -78 and ST.qcd_b0(6, 0, True) == 3 and ST.qcd_b0(6) == Fraction(7)
