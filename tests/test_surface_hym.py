"""SM-08: the compact uniformising metric of X0(143) — Gauss–Bonnet gate; Green's function flux consistency."""
import numpy as np

from mtft.surface import hym as HY, rrspace as RR


def test_gauss_bonnet_and_greens_flux():
    g = HY.gauss_bonnet_gate(143, Y0=2.0, h=0.25, nx=6)
    assert g["newton_iterations"] <= 12 and g["residual"] < 1e-8
    assert abs(g["compact_area"] / (48 * np.pi) - 1) < 0.01                 # 48 pi within the mesh error (56 pi would fail)
    assert abs(g["hyperbolic_area"] / (56 * np.pi - 2) - 1) < 0.01
    cm = RR.cm_classes(); P = cm["P"]
    gr, Ac = HY.greens_functions(g["cx"], g["faces"], g["K"], g["M"], g["b"], g["u"], {"P1": P[1]}, {"cusp0": 143}, Y0=2.0)
    e2u = np.exp(2 * g["u"]); Me = g["M"] @ e2u
    assert abs(float(Me @ gr["P1"])) < 1e-6 * abs(Me).sum() and abs(float(Me @ gr["cusp0"])) < 1e-6 * abs(Me).sum()   # zero-mean normalisation
    assert gr["P1"].min() < -0.3 and gr["cusp0"].min() < -0.3                  # logarithmic wells at the sources


def test_neumann_solve_with_disconnected_slivers():
    # a fine mesh where isolated sliver triangles appear (two 3-node components at h=0.15): the pinned solve must still work
    import scipy.sparse as sps
    K = sps.csc_matrix(np.array([[1, -1, 0, 0], [-1, 1, 0, 0], [0, 0, 1, -1], [0, 0, -1, 1]], float)); Me = np.array([1.0, 1.0, 1e-6, 1e-6])
    G = HY.solve_neumann(K, np.array([1.0, -1.0, 0.0, 0.0]), Me)
    assert abs((G[0] - G[1]) - 1.0) < 1e-12 and abs(float(Me @ G)) < 1e-12


def test_vectorised_character_reduction_matches_pointwise():
    cm = RR.cm_classes(); pts = np.concatenate([cm["Z"], cm["P"]]); d13 = RR.load_chi13(); F4 = d13["F4"].astype(float)
    w, jt, dm = HY.reduce_gamma0_vec(pts)
    q = np.exp(2j * np.pi * w); qp = q[:, None] ** np.arange(F4.shape[1])[None, :]
    chi = np.array([1.0 if pow(int(d), 6, 13) == 1 else -1.0 for d in dm])
    Vvec = np.array([chi * jt ** (-4) * (qp @ F4[i]) for i in range(40)]); Vref, _ = RR.chi13_values(pts, 4, F4)
    assert np.abs(Vvec - Vref).max() / np.abs(Vref).max() < 1e-10
