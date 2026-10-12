"""GS-01: relaxation, the slepton direction, the vector-like symmetry."""
from mtft.matter import green_schwarz as G


def test_relaxation_and_fixed_moduli():
    m = G.minimum_dynamical_xi(n=120); assert m["v2"] == 0.0 and m["xi"] == 0.0 and abs(m["V"]) < 1e-15
    assert abs(G.minimum_fixed_xi() - 8 / 13) < 1e-12


def test_slepton_direction_and_u1v():
    assert abs(G.su4_shift("lepton") + 3 / 8) < 1e-12 and abs(G.su4_shift("colour") - 1 / 8) < 1e-12
    r = G.report(); assert r["cubic_PhiAA"] and not r["dirac_AA"] and not r["kk_Phi_A"]
