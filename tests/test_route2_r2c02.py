import numpy as np, sympy as sp, pytest
from mtft.research import vector_higgs as VH

def test_vector_higgs_tachyon_exact_identity():
    r = VH.tachyon_mass2(-6); R = r["R"]; assert r["cancellation"] and sp.simplify(r["m2_LLL"] + sp.Rational(1, 4) / R ** 2) == 0
    for d in (-3, -6, -12, 3): assert sp.simplify(VH.tachyon_mass2(d)["m2_LLL"] + abs(d) / (24 * VH.tachyon_mass2(d)["R"] ** 2)) == 0      # -|d|/((2g-2) R^2)
    lam = np.array([1.25] * 18 + [1.8]); m2 = VH.mass_operator_spectrum(lam, -6, 48 * np.pi); assert np.allclose(m2[:18], -0.25) and m2[18] > 0

@pytest.mark.slow
def test_vector_higgs_spectrum_fem_and_yukawa_gauge_ratio():
    from mtft.surface import magnetic as MG
    s = VH.m1_vector_higgs_spectrum(0.3, 5, nev=22); assert abs(s["tachyon_mean"] / s["tachyon_exact"] - 1) < 0.06 and s["first_massive"] > 0 and s["multiplicity"] == 18
    r = MG.m1_fem_yukawa(0.3, 5, n_kk=1); q = VH.yukawa_gauge_ratio(r["Y"], r["A_c"], n=1500); assert 0.4 < q["s_max_percentiles"][1] < 0.9      # DIAGNOSTIC band
