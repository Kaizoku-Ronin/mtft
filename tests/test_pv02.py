"""PV-02 (2026-09-29): mtft.surface.blowup (B3) and the exact claims of the arithmetic-supersymmetry note — twin-prime uniqueness of
X0(143), the SUSY-twist chirality e - 24 b, the D-term identity E_split - E_HYM = slope variance / r, and the D-flat Standard-Model search."""
import sympy as sp
from fractions import Fraction
from mtft.surface import blowup as B
from mtft.research import product_surface as PS, vacuum_energy as VE, anomaly8d as A8


def test_lattice_and_canonical_class():
    assert B.dot((1, 0), (0, 1)) == 1 and B.dot((1, 0), (1, 0)) == 0 and B.dot((0, 0, 1, 0), (0, 0, 1, 0)) == -8 and B.dot((0, 0, 0, 1), (0, 0, 0, 1)) == -1
    assert B.dot((a := 3, 1), B.K_TILDE) == 24 * 1 and B.dot((0, 0, 0, 1), B.K_TILDE) == -1
    assert not B.is_spin() and B.dot(B.K_TILDE, B.K_TILDE) == -1 and B.chi((0, 0, 0)) == 0


def test_susy_twist_chirality_is_euler_characteristic_times_torus_flux():
    box = range(-6, 7)
    assert all(B.net_chirality((x, y, 0)) == -24 * y == PS.adjoint_net_chirality((x, y), (0, 0))["net_chirality"] for x in box for y in box)
    assert all(PS.adjoint_net_chirality((x, y))["net_chirality"] == 0 for x in box for y in box)                     # spin twist on the product
    assert B.net_chirality((-1, 0, 3)) == 3 and B.slope((-1, 0, 3), 5, 3, 1) == 0 and B.ample(5, 3, 1)               # the D-flat three-family block


def test_chi_matches_product_riemann_roch():
    for a in range(-5, 6):
        for b in range(-3, 4):
            assert B.chi((a, b, 0)) == (a - 12) * b


def test_twin_prime_uniqueness_of_143():
    def genus_sqfree(ps):
        mu = 1; n2 = 1; n3 = 1
        for p in ps:
            mu *= p + 1; n2 *= 1 + sp.jacobi_symbol(-4 % p, p); n3 *= 1 + sp.jacobi_symbol(-3 % p, p)
        return 1 + sp.Rational(mu, 12) - sp.Rational(n2, 4) - sp.Rational(n3, 3) - sp.Rational(2 ** len(ps), 2)
    hits = []
    for k in range(1, 200):
        p, q = 6 * k - 1, 6 * k + 1
        if sp.isprime(p) and sp.isprime(q):
            g = genus_sqfree([p, q]); assert g == k * (3 * k + 1) - 1
            if g - 1 == (p + q) // 2: hits.append((p, q, g))
    assert hits == [(11, 13, 13)] and sp.factor(3 * sp.Symbol("k") ** 2 - 5 * sp.Symbol("k") - 2) == (sp.Symbol("k") - 2) * (3 * sp.Symbol("k") + 1)


def test_d_term_identity_split_minus_hym_is_slope_variance():
    r = VE.R
    for m in ("S1", "S2"):
        st = VE.STACK_MODELS[m]; N = sum(n for _, n, _ in st)
        s = [d[0] + d[1] * r for _, _, d in st]; sbar = sum(n * si for (_, n, _), si in zip(st, s)) / N
        var = sum(n * (si - sbar) ** 2 for (_, n, _), si in zip(st, s)) / r
        assert sp.simplify(VE.flux_energy(model=m)["E_over_4pi2"] - VE.hym_floor(model=m)["E_HYM_over_4pi2"] - var) == 0


def test_five_means_of_11_and_13():
    a, b = 11, 13
    assert (Fraction(2 * a * b, a + b), a * b, Fraction(a + b, 2), Fraction(a * a + b * b, 2), Fraction(a * a + b * b, a + b)) == (Fraction(143, 12), 143, 12, 145, Fraction(145, 12))
    assert 143 + 145 == 2 * 144 and PS.SPIN_TWIST == (12, 0)


def test_sm_phi_conditions_and_solution_space():
    sol = B.sm_solution_space(); pc, pa = sol["free"]
    for c_, a_ in ((0, 0), (0, 3), (2, -5), (-7, 11)):
        phis = {"c": c_, "L": int(sol["phi_L"].subs(pc, c_)), "B": int(sol["phi_B"].subs(pc, c_)), "A1": a_, "A2": int(sol["phi_A2"].subs({pc: c_, pa: a_}))}
        rec = B.sm_phi_conditions(phis); assert rec["exact_SM"] and rec["vector_like_pairs"]["doublets"] >= 6
    assert not B.sm_phi_conditions({"c": 0, "L": -3, "A1": 3, "A2": 3, "B": 3})["exact_SM"]                       # phi_A1 + phi_A2 must be 3


def test_collinear_family_and_search():
    fam = B.collinear_family()
    assert fam["exact_SM"] and fam["traceless"] and fam["d_flat_locus"] == sp.Rational(1, 3) and fam["net_per_unit_charge"] == -3
    assert fam["vector_like_pairs"] == {"u^c": 6, "e^c": 6, "doublets": 9, "sterile_singlets": 15}
    s = B.d_flat_sm_search(n_box=6)
    assert s["n_solutions"] > 0 and all(B.sm_phi_conditions({**sol["phi"]})["exact_SM"] for sol in s["solutions"])
    for sol in s["solutions"]:
        AE, AX, eps = sol["kahler_normal_(A_E,A_X,eps)"]
        assert B.ample(AX, AE, eps) and all(B.slope(d, AX, AE, eps) == 0 for d in sol["stacks"].values())
    assert any(sol["collinear"] for sol in s["solutions"]) and min(sol["extras"] for sol in s["solutions"]) <= 36


def test_polystable_locus_of_collinear_family():
    fam = B.collinear_family(); stacks = tuple((n, {"c": 3, "L": 2, "B": 1, "A1": 1, "A2": 1}[n], d) for n, d in fam["stacks"].items())
    loc = B.polystable_locus(stacks)
    assert len(loc) == 1 and loc[0]["eps_over_A_E"] == sp.Rational(1, 3)
