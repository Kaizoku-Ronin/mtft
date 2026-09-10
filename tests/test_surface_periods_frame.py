"""v0.27.2: the periods-frame integer map; three frames are one lattice."""
import numpy as np

from mtft.surface import frozen as FR, intertwiner as IT


def test_periods_frame_gates():
    d = FR.x0143()
    assert all(d["gates"].values())
    assert all(k in d["gates"] for k in ("periods_X_unimodular", "periods_J_true_equals_minus_periods_J",
                                          "intertwiner_Pi_unimodular"))
    X = IT.periods_frame_map()
    assert abs(round(np.linalg.det(X.astype(float)))) == 1


def test_full_hodge_structure_agrees_across_frames():
    import mtft.periods as P
    d = FR.x0143()
    X = d["X_periods"].astype(float)
    Jm = P.hodge_complex_structure(40)
    Jp = np.array([[float(Jm[i, j]) for j in range(26)] for i in range(26)])
    Jt = X @ d["J_true"] @ np.linalg.inv(X)
    assert np.linalg.norm(Jt + Jp) / np.linalg.norm(Jp) < 1e-12
