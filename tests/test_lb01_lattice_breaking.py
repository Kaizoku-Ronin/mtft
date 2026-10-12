"""LB-01: lattice breaking at the blow-up scale and the bound on jet Yukawas."""
from mtft.matter import lattice_breaking as LB


def test_lattice_data():
    d = LB.lattice_data()
    assert abs(d["j"] + 262144 / 1859) < 1e-8 and abs(d["G4"] - 1.40162) < 1e-4 and abs(d["G6"] - 3.85955) < 1e-4
    assert d["shell_angles"] == [-116.0, -64.0, 0.0, 64.0, 116.0, 180.0] and abs(d["sixfold_sum"] - 5.6612) < 1e-3


def test_parity_and_breaking():
    assert LB.parity_channels() == [(3, 0, 2, 6), (4, 0, 1, 6), (4, 1, 2, 7)]
    d = LB.lattice_data()
    assert abs(LB.breaking(6, data=d) - 7.684e-4) < 2e-6 and abs(LB.breaking(6, across=True, data=d) - 0.04918) < 1e-4


def test_jet_top_bound():
    top = LB.top_yukawa_at_crossing(); b = LB.jet_yukawa_bound(0.14)
    assert 0.68 < top["y_t"] < 0.71 and b < 0.0100 and top["y_t"] / b > 70 and abs(top["y_b"] - 0.0104) < 3e-4
