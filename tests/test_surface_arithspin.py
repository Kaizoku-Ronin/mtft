"""ARITH-SPIN-01 / GAL-01 gates: cuspidal group, K ~ 6D, theta characteristics (5, 2), unit constants, purity."""
import numpy as np

from mtft.surface import arithspin as AS


def test_cuspidal_group_and_theta_characteristics():
    lat = AS.eta_quotient_lattice(143)
    assert list(lat["exponents"].tolist()).count([-1, -1, 1, 1]) == 1          # the unit u is found
    H = lat["hermite"]
    g = AS.cuspidal_group(H); assert g["order"] == 4200 and sorted(g["invariant_factors"]) == [1, 10, 420] or g["invariant_factors"] == [420, 10, 1]
    assert (AS.cusp_difference_order(H, 0, 3), AS.cusp_difference_order(H, 1, 3), AS.cusp_difference_order(H, 2, 3)) == (420, 60, 70)
    th = AS.x0143_theta_characteristics(H)
    assert len(th["two_torsion"]) == 3 and all(all(c["al_invariant"].values()) for c in th["characteristics"])
    assert AS.h0_cuspidal_class([3, 3, 3, 3]) == 5 and AS.h0_cuspidal_class([6, 6, 0, 0]) == 2 and AS.h0_cuspidal_class([0, 0, 0, 0]) == 13
    assert AS.effective_representative(H, [-24, 6, 30, 0]) == [3, 3, 3, 3]


def test_unit_constants_and_purity():
    k = AS.functional_constants()
    assert abs(k["C_W13"] + 1 / 13) < 1e-8 and abs(k["c_W11"] + 1) < 1e-8
    pu = AS.three_family_purity()
    assert pu["u_squared_max_dev_from_-1/13"] < 1e-8 and pu["signs_present"] == [-1, 1] and pu["h0_L3"] == 3
