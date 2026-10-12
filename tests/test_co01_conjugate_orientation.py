"""CO-01: no conjugate condensate; what distinguishes right from left."""
from mtft.matter import conjugate_orientation as C


def test_closed_form_and_flat_family():
    m = C.d_term_minimum(); assert abs(m["v2_over_xi"] - 8 / 13) < 1e-12 and abs(m["residual_U1"] + 5 / 13) < 1e-12
    v = [C.d_term_potential(x, 1.0) for x in (0.6, 8 / 13, 0.63)]; assert v[1] < v[0] and v[1] < v[2]
    assert abs(C.d_term_potential(1.0, 1.0, 0.4) - C.d_term_potential(0.6, 1.0, 0.0)) < 1e-12
    lift = C.lifted_minimum(n=150); assert lift["vbar2"] == 0.0 and lift["v2"] > 0.5


def test_asymmetries():
    s = C.scalar_masses(); assert abs(s["e_c_partner"] - 10 / 13) < 1e-12 and abs(s["coloured_partner"] - 2 / 13) < 1e-12
    assert 0.88 < C.gauge_loop_asymmetry() < 0.92
