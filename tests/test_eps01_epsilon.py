"""EPS-01: what sets epsilon — exact selection rules at the blown-up point and the gamma-flux estimate."""
import numpy as np
from mtft.matter import epsilon as E


def test_selection_rules():
    assert E.allowed_entries([0]) == {(0, 0)}                                          # torus-constant D-term: rank one
    assert E.d_charges_from_condensate([0, 1, 2]) == [0, 1, 2]
    assert E.allowed_entries([0, 1, 2]) == {(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0)}
    assert E.plane_wave_density_is_constant() < 1e-12                                 # single-band gamma mode: torus-constant density
    assert all(E.outer_product_rank(s) == 1 for s in range(5))                        # rotating the families keeps rank one


def test_texture():
    p = E.texture_properties(0.1)
    s = p["singular_values"]; assert abs(np.prod(s) - abs(p["det"])) < 1e-12 and 0.005 < s[1] / s[0] < 0.03 and 5e-5 < s[2] / s[0] < 5e-4
    Y = E.texture(0.1); assert Y[1, 2] == 0 and Y[2, 1] == 0 and Y[2, 2] == 0


def test_up_down_and_size():
    assert abs(E.eps_ratio_needed() - 2.305) < 2e-3
    rho = E.arakelov_cusp_density(); assert abs(rho - 0.02434) < 2e-4
    lo, hi = E.r_for_eps(0.14, rho=rho); assert 200 < lo < 270 and 2000 < hi < 2600
    a, b = E.eps_gamma(176, rho=rho); assert 0.14 < a < 0.18 and 0.45 < b < 0.55
