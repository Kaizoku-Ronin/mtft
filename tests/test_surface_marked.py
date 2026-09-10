"""v0.28.0: marked geometry instrument and hardened oldsector, pinned to Astra's exact tables."""
import sympy as sp
import pytest

from mtft.surface import marked as MK, oldsector as OS


def test_marked_observables_separate_isospectral_inputs():
    S = OS.two_prime_sector(13, 17, 4, -2)
    Gp, Gq = S["local_p"]["G"], S["local_q"]["G"]
    Lf, Li = MK._orthonormal(Gp, Gq)
    Z, I2 = sp.diag(1, -1), sp.eye(2)
    hsyn = Li.T * (60 * sp.kronecker_product(Z, I2) + 48 * sp.kronecker_product(I2, Z)) * Lf.T
    Hc = Li.T * OS.connected_part(S["H_pq"], Gp, Gq) * Lf.T
    assert sorted(hsyn.eigenvals()) == sorted(Hc.eigenvals()) == [-108, -12, 12, 108]
    da, ds = MK.decomposition(Hc, Gp, Gq), MK.decomposition(hsyn, Gp, Gq)
    assert (da["centered_norm2"], da["connected_norm2"], da["eta"]) == (5904, 5904, 1)
    assert (ds["local_p_norm2"], ds["local_q_norm2"], ds["connected_norm2"], ds["eta"]) == (3600, 2304, 0, 0)


def test_maps_readout_and_two_objectives():
    S = OS.two_prime_sector(13, 17, 4, -2)
    maps = MK.degeneracy_maps(S)
    assert all(maps["gates"].values())
    a = sp.Symbol("alpha")
    H = S["H_prod"] + a * S["skew_correction"]
    r = MK.cross_map_readout(S, H, maps)
    A, R = S["local_p"]["A"], S["local_p"]["R"]
    assert (r["K"][0][0] + 34 * A).is_zero_matrix and sp.simplify(r["K"][0][1] - 146 * A - 160 * a * R).is_zero_matrix
    obj = MK.two_objectives(S)
    assert obj["least_leakage_trace"] == 7200 * a**2 + sp.Rational(409600, 81)   # minimum at alpha = 0
    assert obj["herm_hecke_distance2"] == sp.Rational(57600, 41) * (a - 1)**2       # zero at alpha = 1
    assert obj["readout_alpha"] == a


def test_oldsector_hardening():
    S = OS.two_prime_sector(13, 17, 4, -2)
    sel = OS.hermitian_hecke_selection(S)
    assert sel["status"] == "unique" and sel["alpha_selected"] == 1
    c = OS.good_prime_replacement_control(S, -2)
    assert c["skew_correction_zero"] and c["connected_term_zero"]
    with pytest.raises(ValueError):
        OS.two_prime_sector(13, 13, 4, 4)
    with pytest.raises(ValueError):
        OS.local_block(12, 1)
    with pytest.raises(ValueError):
        OS.local_block(13, 9)                     # Hasse bound
