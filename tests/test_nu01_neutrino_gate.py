"""NU-01: the neutrino gate — inventory, localisation, Atkin–Lehner at the cusps, seesaw requirements."""
from fractions import Fraction
from mtft.matter import neutrino_gate as NG


def test_inventory_and_the_cheap_gamma_flux():
    inv = NG.inventory()["rows"]
    lep, nuc, h1, h2 = inv["nu_L (A1 <- L)"], inv["nu^c (A2 <- A1)"], inv["bidoublet (L <- B)"], inv["bidoublet (L <- A2)"]
    assert (lep["modes"], lep["conjugate_modes"], abs(lep["net"]), lep["vector_like_pairs"]) == (19, 22, 3, 19) and lep["Y"] == Fraction(-1, 2)
    assert (nuc["modes"], nuc["conjugate_modes"], abs(nuc["net"]), nuc["vector_like_pairs"]) == (19, 22, 3, 19) and nuc["B_minus_L"] == 1 and nuc["Y"] == 0
    assert (h1["modes"], h1["conjugate_modes"], h1["gamma_charge"]) == (21, 15, 0) and h1["B_minus_L"] == 0 == h2["B_minus_L"]
    assert NG.inventory()["BL_neutral_singlets"] == 0
    assert sum(NG.localised_extras(NG.GAMMA_CHEAP).values()) == 256 < sum(NG.localised_extras(NG.GAMMA_PS01).values()) == 384
    v = NG.fi_check()["vacua"]; assert len(v) == 1 and v[0]["rep"] == "(4,1,2)/(4bar,1,2)"          # the cheap flux still gives the PS condensate alone
    assert sum({"c": 3, "A1": 1, "L": 2, "B": 1, "A2": 1}[s] * NG.GAMMA_CHEAP[s] for s in NG.GAMMA_CHEAP) == 0


def test_localisation_and_atkin_lehner():
    assert NG.localisation_sets(1)["bosons_are_the_cusps"] and NG.localisation_sets(2)["bosons"] == 16 and NG.localisation_sets(2)["disjoint"]
    al = NG.atkin_lehner_on_cusps()
    assert al["simply_transitive"] and al["character_multiplicities"] == {(1, 1): 1, (1, -1): 1, (-1, 1): 1, (-1, -1): 1}
    assert al["action"][143][143] == 1 and al["action"][11][143] == 13


def test_seesaw_requirements():
    y = NG.up_yukawas_at()
    assert 0.66 < y["y_t"] < 0.73 and 115 < y["m_D_GeV"]["3"] < 127
    req = NG.seesaw_requirements()
    assert 2e14 < req["3"]["typeI_M_R_GeV"] < 4e14 and 0.003 < req["3"]["inverse_mu_GeV"] < 0.006 and 3e-9 < req["3"]["mu_over_M"] < 5e-9
    assert 40 < req["2"]["inverse_mu_GeV"] < 90 and req["2"]["mu_over_M"] < 1e-4
