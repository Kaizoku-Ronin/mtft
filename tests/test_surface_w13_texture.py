import numpy as np
from mtft.surface import hym as HY, rrspace as RR

def test_degeneracy_theorem_in_corrected_normalisation():
    d = HY.load_m2_tensors(); Yres = RR.up_yukawa_M1(); G = HY.w13_graded_normalisation(Yres, d["N_fam_untwisted"], d["N_H_untwisted"])
    assert G["gram_offgrade_family"] < 1e-2 and G["gram_offgrade_higgs"] < 1e-2           # W13-invariant metric (mesh-level block diagonality)
    rng = np.random.default_rng(1)
    for _ in range(5):
        v = rng.standard_normal(G["Ye"].shape[2]) + 1j * rng.standard_normal(G["Ye"].shape[2]); s = np.sort(np.linalg.svd(np.einsum("ijk,k->ij", G["Ye"], v), compute_uv=False))[::-1]
        assert abs(s[1] / s[0] - 1) < 1e-6 and s[2] / s[0] < 1e-6                           # even VEV: (1, 1, 0); tolerance 1e-6 over a measured 1e-9
    sc = HY.w13_epsilon_scan(G["Ye"], G["Yo"], eps_values=(1e-3, 1e-1), n=100)
    assert 3e-4 < sc[1e-3]["median_m1_m3"] < 3e-3 and sc[1e-3]["median_m2_m3"] > 0.99      # m1/m3 ~ eps, heavy pair unsplit
