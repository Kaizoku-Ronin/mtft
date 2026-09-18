import numpy as np, pytest
from mtft.surface import magnetic as MG

@pytest.mark.slow
def test_landau_levels_realise_h0():
    mm = MG.MagneticMesh(0.25, 6); rng = np.random.default_rng(0)
    r = mm.divisor_spectrum(30, [(ti, 1) for ti in rng.choice(mm.nT, 30, replace=False)], 24); e = r["eigenvalues"]     # h^0 = 18 for every degree-30 bundle
    assert abs(e[:18].mean() / r["landau"] - 1) < 0.03 and e[:18].std() < 0.05 and e[18] - e[17] > 0.2
    r = mm.divisor_spectrum(15, [(ti, 1) for ti in rng.choice(mm.nT, 15, replace=False)], 8); e = r["eigenvalues"]      # h^0 = 3 for every degree-15 bundle (Riemann–Roch)
    assert abs(e[:3].mean() / r["landau"] - 1) < 0.03 and e[:3].std() < 0.02


@pytest.mark.slow
def test_fem_yukawa_reproduces_corrected_normalisation():
    """Cross-pipeline gate (KK-TOWER-02): FEM eigenmodes (no Cholesky) vs the CC-26-corrected modular-form pipeline."""
    from mtft.surface import hym as HY
    r = MG.m1_fem_yukawa(0.25, 6, n_kk=1); R = MG.ratio_distribution(r["Y"], n=1500)
    d = HY.load_m2_tensors(); Yc = HY.normalise_yukawa(d["Y_M1_up_raw"], d["N_fam_untwisted"], d["N_fam_untwisted"], d["N_H_untwisted"]); Rc = MG.ratio_distribution(Yc, n=1500)
    for col in (0, 1):
        assert abs(np.log(np.median(R[:, col]) / np.median(Rc[:, col]))) < 0.15        # medians agree within 15% (measured 1–3% at h = 0.2)
    assert np.median(R[:, 0]) > 0.05                                                    # and both are O(1): the retracted 1e-8 is excluded by 6 orders
