import sympy as sp, pytest
from mtft.research import product_surface as PS, theta_torus as TT

def test_full_bidegree_rule_and_two_solutions():
    assert not PS.yukawa_bidegree_rule((3, -1), (-1, -3))["allowed"]            # the R2C-08 "solution" sums to (2,2): CC-32
    assert PS.yukawa_bidegree_rule((3, 1), (1, -3))["higgs_bidegree"] == (1, 0)
    sols = PS.surface_solutions(); keys = {(r["Q"], r["u"], r["H"], r["higgs_modes"], r["slope_free_locus"] / PS.A_E) for r in sols}
    assert keys == {((3, 1), (1, -3), (-4, 2), 32, 2), ((1, -3), (3, 1), (-4, 2), 32, 2), ((-3, 1), (1, 3), (2, -4), 4, sp.Rational(1, 2)), ((1, 3), (-3, 1), (2, -4), 4, sp.Rational(1, 2))}
    assert all(r["slope_free_locus"] is None or r["higgs_modes"] == 0 for r in PS.scan_family_pairs_full() if {abs(r["Q"][0]), abs(r["u"][0])} == {3})   # curve x curve never works
    assert TT.up_mass_rank_bound()["rank_bound"] == 2

@pytest.mark.slow
def test_theta_torus_factors():
    import mpmath as mp
    tau = TT.tau_143a1(); assert abs(mp.re(tau) - 0.5) < 1e-12 and abs(mp.im(tau) - 1.0232745927) < 1e-8
    s1 = TT.torus_factor(1, 2, tau, N=16); assert s1["gram_offdiag"] < 1e-12 and s1["closure_residual"] < 1e-10 and s1["rank"] == 2
    assert abs(s1["singular_values"][0] - 1.103) < 0.02 and abs(s1["singular_values"][1] - 0.859) < 0.02
    s2 = TT.torus_factor(1, 3, tau, N=16); assert s2["rank"] == 3 and abs(s2["singular_values"][0] - 1.124) < 0.02
