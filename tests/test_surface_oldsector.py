"""v0.27.4: exact two-prime old sector (11a1 at 13, 17) and the Herm(Hecke) selection rule."""
import sympy as sp

from mtft.surface import oldsector as OS

x = sp.Symbol("x")


def test_composite_correspondence_numbers():
    S = OS.two_prime_sector(13, 17, 4, -2)
    assert S["cross_commutators_zero"]
    assert sp.factor((x * sp.eye(4) - S["H_pq"]).det()) == sp.factor((x**2 - 92 * x - 1584) * (x**2 + 100 * x - 1584))
    C = OS.connected_part(S["H_pq"], S["local_p"]["G"], S["local_q"]["G"])
    assert sorted(C.eigenvals().keys()) == [-108, -12, 12, 108] and OS.schmidt_rank(C) == 2
    Cp = OS.connected_part(S["H_prod"], S["local_p"]["G"], S["local_q"]["G"])
    assert set(Cp.eigenvals().keys()) == {-48, 48} and OS.schmidt_rank(Cp) == 1
    assert OS.connected_part(S["H_add"], S["local_p"]["G"], S["local_q"]["G"]).is_zero_matrix


def test_hermitian_hecke_selects_alpha_one():
    S = OS.two_prime_sector(13, 17, 4, -2)
    r = OS.hermitian_hecke_selection(S)
    assert r["dim_Herm_Hecke"] == 4 and r["alpha_selected"] == 1
    # good-prime control: skew part vanishes, so no interaction term exists
    L = OS.local_block(13, 4)
    assert not (L["R"]).is_zero_matrix           # bad-prime U_13 is non-normal
    assert (L["R"].T * L["G"] + L["G"] * L["R"]).is_zero_matrix   # R is G-skew
