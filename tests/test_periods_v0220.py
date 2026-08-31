"""v0.22: involutions, oldspace torus, Hamiltonian layer, Bergman channels."""
import math
import mpmath as mp

import mtft
import mtft.periods as P
from mtft.periods import gates as G


def test_version_and_data_ship():
    assert mtft.__version__ == "0.24.0"
    assert P.data_path("X0_143_atkin_lehner_v022.json").exists()


def test_involutions_exact():
    r = G.gate_involutions()
    assert r["census"] == (1, 6, 5, 1)
    assert r["eps"][11] == {"ell": "+", "q4": "-", "q6": "+"}
    assert r["eps"][13] == {"ell": "+", "q4": "+", "q6": "-"}
    assert r["route2"] == {"ell": 2, "ghost": 2, "q4": 0, "q6": 0}


def test_oldtorus_chain():
    r = G.gate_oldtorus()
    assert r["l9_index"] == 9 and r["mod3_rank"] == 2
    assert r["polarization"]["type"] == (2, 18)
    assert abs(r["entropy"] - 2 * math.log(2 + math.sqrt(5))) < 1e-14


def test_hamiltonian_layer():
    r = G.gate_hamiltonian(50)
    assert abs(r["width"]["rho"] - 0.1234286299) < 1e-8
    assert abs(r["distance"]["rho"] - 0.4248813827) < 1e-8
    assert r["width"]["new_new"] > 0.85 and r["distance"]["new_new"] > 0.85
    assert r["width"]["intra_block"] > 0.7 > 0.5 > r["distance"]["intra_block"]


def test_bergman_channels():
    r = G.gate_bergman_channels(40)
    assert abs(mp.mpf(r["crossover_ratio"])
               - mp.mpf("2.302140221833918907")) < mp.mpf("1e-15")


def test_star_orbit_symmetry():
    orb = P.star_charge_orbit()
    assert len(orb) == 4 and orb[0][2] == 1


def test_degree_null_control_exact():
    rep = P.channel_report("degree", 50)
    assert rep["hamiltonian_antilinear_fraction"] < 1e-13


def test_star_odd_residual_api():
    r = P.cp_channel_report("degree", 50)
    assert math.isnan(r["commutator_star_odd_rel"])
    assert r["commutator_star_odd_abs"] < 1e-13
    w = P.cp_channel_report("width", 50)
    assert w["commutator_star_odd_rel"] < 1e-10
