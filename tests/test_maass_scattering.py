"""mtft.maass.scattering: the scattering matrix of the cusps of Gamma_0(143) and its Eisenstein series."""
import mpmath as mp
import pytest

from mtft.maass import scattering as S


def _mx(M):
    return max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))


def test_scattering_matrix_identities():
    mp.mp.dps = 40
    assert S.cusps(143) == ["inf", "1/11", "1/13", "0"] and S.cusps(11) == ["inf", "0"] and S.cusps(1) == ["inf"]
    fe = max(_mx(S.scattering_matrix(s) * S.scattering_matrix(1 - s) - mp.eye(4)) for s in (mp.mpf("0.31"), mp.mpc("0.7", "1.3"), mp.mpc("2.2", "-0.4")))
    un = max(_mx(S.scattering_matrix(mp.mpc(0.5, t)) * S.scattering_matrix(mp.mpc(0.5, t)).H - mp.eye(4)) for t in (0.3, 2.0, 7.5))
    assert 0 < fe < mp.mpf(10) ** -30 and 0 < un < mp.mpf(10) ** -30          # the record: 3.5e-41, 4.6e-41
    P = S.scattering_matrix(mp.mpf("1.3"))
    assert _mx(P - P.T) == 0 and abs(P[0, 0] - P[3, 3]) == 0                  # symmetric; the flip of both types
    assert _mx(P - mp.eye(4)) > 0.1                                           # and not the identity
    mp.mp.dps = 15


def test_phi_at_one_half_is_minus_the_identity():
    mp.mp.dps = 40
    h = mp.mpf(10) ** -25
    assert _mx(S.scattering_matrix(mp.mpf(1) / 2 + h) + mp.eye(4)) < mp.mpf(10) ** -20      # the limit
    assert _mx(S.scattering_matrix(mp.mpf(1) / 2) + mp.eye(4)) == 0                         # the value returned
    assert abs(S.phi(mp.mpf(1) / 2 + h) + 1) < mp.mpf(10) ** -20
    assert _mx(S.local_matrix(11, mp.mpf(1) / 2) - mp.eye(2)) < mp.mpf(10) ** -35
    mp.mp.dps = 15


def test_residue_at_one_is_one_over_the_area():
    with mp.workdps(80):
        h = mp.mpf(10) ** -30
        res = S.scattering_matrix(1 + h) * h
        dev = max(abs(res[i, j] - 3 / (168 * mp.pi)) for i in range(4) for j in range(4))
        assert dev < mp.mpf(10) ** -20                                        # the record: 5.2e-32
        assert abs((S.scattering_matrix(1 + h, 11) * h)[0, 1] - 3 / (12 * mp.pi)) < mp.mpf(10) ** -20


def test_eisenstein_series_are_invariant_and_vanish_at_the_edge():
    mp.mp.dps = 20
    z = mp.mpc("-0.31", "0.05")
    w = (12 * z + 1) / (143 * z + 12)                                         # the Pell unit of Gamma_0(143)
    a, b = S.eisenstein(z, 2), S.eisenstein(w, 2)
    assert max(abs(a[k] - b[k]) / abs(a[k]) for k in a) < mp.mpf(10) ** -15
    assert abs(a["inf"] - a["0"]) / abs(a["inf"]) > 0.1                       # four different series
    u = -1 / (143 * z)                                                        # the Fricke involution swaps inf and 0
    c = S.eisenstein(u, 2)
    assert abs(c["0"] - a["inf"]) / abs(a["inf"]) < mp.mpf(10) ** -15 and abs(c["1/11"] - a["1/13"]) / abs(a["1/13"]) < mp.mpf(10) ** -15
    ratios = []
    for h in (mp.mpf(10) ** -6, mp.mpf(10) ** -9):
        E = S.eisenstein(mp.mpc("0.123", "0.6"), mp.mpf(1) / 2 + h)
        ratios.append(max(abs(v) for v in E.values()) / h)
    assert ratios[1] < 10 and abs(ratios[0] - ratios[1]) < 1e-3               # E_a(z, 1/2 + h) = O(h): all four vanish at s = 1/2
    assert abs(ratios[1] - mp.mpf("4.573")) < 1e-2                            # the record's value at this point
    mp.mp.dps = 15


def test_error_paths():
    for fun in (S.cusps, lambda N: S.scattering_matrix(2, N), lambda N: S.eisenstein(1j, 2, N)):
        with pytest.raises(ValueError):
            fun(12)
