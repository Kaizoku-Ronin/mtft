"""NU-02: the modular parametrisation of 143a1, the fermionic localisation fibres, the Atkin–Lehner action on them, and the suppression estimates."""
import math
import numpy as np
import pytest
from mtft.matter import neutrino_mu as NM


@pytest.fixture(scope="module")
def MP():
    return NM.ModularParametrisation()


@pytest.fixture(scope="module")
def fib(MP):
    return NM.fibres(MP)


def test_newform_and_lattice(MP):
    assert MP.a[1:14].astype(int).tolist() == [1, 0, -1, -2, -1, 0, -2, 0, -2, 0, -1, 2, -1]
    assert abs(MP.Omega - 1.96992316457) < 1e-9
    lc = NM.lattice_check(MP)
    for row in lc["periods"]: assert np.allclose(row, np.round(row), atol=1e-8) and np.any(np.abs(np.round(row)) > 0)
    assert np.allclose(lc["W143_shift"], 0, atol=1e-9)


def test_fibres_are_free_atkin_lehner_orbits(fib):
    for name, v in fib.items():
        assert v["distinct"] and v["AL_free"] and len(v["points"]) == 4
        assert max(v["stretch"]) - min(v["stretch"]) < 1e-6                         # Atkin–Lehner isometries: same stretch on the fibre
    assert abs(fib["real 1/2"]["rho"] - 3.5891) < 1e-3 and abs(fib["tau/2"]["rho"] - 3.3973) < 1e-3
    assert abs(fib["real 1/2"]["stretch"][0] - 0.1965) < 1e-3 and abs(fib["tau/2"]["stretch"][0] - 0.2049) < 1e-3
    assert abs(fib["real 1/2"]["sigma"] - 0.7675) < 1e-3 and abs(fib["tau/2"]["sigma"] - 0.9956) < 1e-3


def test_suppression_requirements(fib):
    req = NM.requirements(fib)
    lo, hi = req["real 1/2"]["third family"]; assert 600 < lo < 700 and 1000 < hi < 1150
    lo, hi = req["tau/2"]["second family"]; assert 90 < lo < 100 and 170 < hi < 185
    a, b = req["tau/2"]["at_r_56pi"]; assert 1e-6 < a < 3e-6 and 5e-5 < b < 7e-5


def test_flux_options():
    fo = NM.flux_options()
    assert fo["nu-optimal (1,1,-2,0,0)"]["nu_c_gamma_charge"] == -1 and fo["nu-optimal (1,1,-2,0,0)"]["d_flat_reps"] == ["(4,1,2)/(4bar,1,2)"]
    assert fo["cheap (1,1,-1,-1,-1)"]["nu_c_gamma_charge"] == -2 and fo["cheap (1,1,-1,-1,-1)"]["extras"] == 256
