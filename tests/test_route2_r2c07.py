import sympy as sp
from mtft.research import product_surface as PS

def test_index_and_slope_decouple_on_the_product_surface():
    fam = PS.product_block(3, -1); assert fam["fermion_index"] == -3 and fam["fermion_zero_modes"] == (0, 3) and fam["c1_squared"] == -6
    H = PS.triangle((3, -1), (3, -1)); assert H["L_H"] == (-6, 2) and H["block"]["fermion_index"] == -12 and H["block"]["vector_zero_modes_H1"] == 36
    assert H["slope_zero_locus"] == [3 * PS.A_E] and sp.solve(sp.Eq(fam["slope"], 0), PS.A_X) == [3 * PS.A_E]                          # same locus
    r = PS.m4_record(); d = sp.Symbol("delta"); assert sp.expand(r["off_locus"]["higgs_slope"] + 2 * d) == 0 and sp.expand(r["off_locus"]["family_vector_slope"] - d) == 0   # opposite signs
    # independent checks: Riemann–Roch on the curve factor (h^1(O(-2 sum P)) = 18 as in SM-08) and the torus (deg 2 -> 2 sections)
    assert PS.curve_cohomology(-6, False) == (0, 18) and PS.torus_cohomology(2) == (2, 0) and PS.curve_cohomology(3, True) == (3, 0)
    # on the curve alone the same triangle has slope = degree = -6 (R2C-06); here index -12 but slope 2(A_X - 3A_E)
    assert sp.expand(H["block"]["slope"] - (2 * PS.A_X - 6 * PS.A_E)) == 0
