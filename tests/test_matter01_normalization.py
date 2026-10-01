"""Basis independence, physical matching, and the null-wave distinction."""
import math
import numpy as np
import pytest
from mtft.matter.normalization import (
    canonical_u1_vertex, four_dimensional_coupling, electroweak_matching,
    breathing_gauge_kinetics, electromagnetic_invariants,
)


def test_u1_vertex_is_universal_after_nonorthogonal_basis_changes():
    rng = np.random.default_rng(143)
    k = np.array([[3,1+0.4j],[1-0.4j,2]],dtype=complex)
    for _ in range(12):
        a = rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))+3*np.eye(2)
        transformed = a.conj().T@k@a
        result = canonical_u1_vertex(transformed,-1,0.30282212)
        assert result["kinetic_residual"] < 2e-12
        assert np.allclose(result["vertex"],-0.30282212*np.eye(2),atol=2e-13,rtol=2e-12)


def test_bad_kinetic_forms_do_not_produce_physical_modes():
    for k in ([[1,2],[2,1]], [[1,1],[0,1]], [[0]], [[float('nan')]]):
        with pytest.raises(ValueError):
            canonical_u1_vertex(k,-1,0.3)


def test_dimensional_matching_scaling():
    assert four_dimensional_coupling(2,16) == 0.5
    assert four_dimensional_coupling(2,64) == 0.25
    for x in (0,-1,float('nan'),float('inf'),True):
        with pytest.raises(ValueError):
            four_dimensional_coupling(x,1)


def test_electroweak_matching_preserves_neutral_higgs_generator():
    # These coefficients are matching-scale inputs. This is not alpha(0).
    g2, gy = 1, math.sqrt(3/5)
    r = electroweak_matching(g2,gy)
    assert r["sin2_theta_w"] == pytest.approx(3/8)
    assert 1/r["e"]**2 == pytest.approx(1+5/3)
    neutral_mass_matrix = np.array([[g2*g2,-g2*gy],[-g2*gy,gy*gy]])
    photon = np.array([gy,g2])/math.hypot(g2,gy)
    assert np.linalg.norm(neutral_mass_matrix@photon) < 1e-15
    for t3,y in ((-0.5,0.5),(-0.5,-0.5),(0,-1),(0.5,1/6)):
        assert g2*photon[0]*t3+gy*photon[1]*y == pytest.approx(r["e"]*(t3+y),abs=1e-15)


def test_ads02_breathing_gauge_coefficient_and_8d_distinction():
    r = breathing_gauge_kinetics(0.2,2)
    assert r["Z"] == pytest.approx(math.exp(0.2))
    assert r["g_eff_over_g0"] == pytest.approx(math.exp(-0.1))
    assert r["mp_dZ_dchi"] == pytest.approx(r["Z"])
    assert breathing_gauge_kinetics(0.2,4)["canonical_exponent"] == pytest.approx(math.sqrt(4/3))
    assert r["implemented_evolution"] is False


def test_null_light_wave_has_energy_without_direct_scalar_f_squared_source():
    plane = electromagnetic_invariants([1,0,0],[0,1,0])
    assert plane["F_squared"] == 0
    assert plane["energy_density"] == 1
    assert plane["poynting"] == [0,0,1]
    standing_antinode = electromagnetic_invariants([2,0,0],[0,0,0])
    assert standing_antinode["F_squared"] == -8
    assert standing_antinode["energy_density"] == 2
