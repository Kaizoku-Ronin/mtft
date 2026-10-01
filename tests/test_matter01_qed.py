"""Matter-01 checks physical identities, not a second copy of the amplitude code."""

import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from mtft.matter import (
    ALPHA_REFERENCE, ELECTRON_MASS_GEV, GEV2_TO_NB, Kinematics,
    boost, centre_of_mass, diagram_amplitudes, differential_cross_section,
    dirac_adjoint, gamma_matrices, integrate_cross_section, minkowski_dot,
    slash, transverse_polarizations, u_spinor,
    unpolarized_amplitude_squared, v_spinor,
)
from mtft.matter.report import benchmark, numerical_gates


M = ELECTRON_MASS_GEV


def _point(beta, c, phi=0):
    return centre_of_mass(2 * M / math.sqrt(1 - beta * beta), c, phi=phi)


def test_clifford_algebra_and_gamma_hermiticity():
    gamma = gamma_matrices()
    ident = np.eye(4)
    for i in range(4):
        for j in range(4):
            target = 2 * (1 if i == 0 else -1) * ident if i == j else np.zeros((4, 4))
            np.testing.assert_allclose(gamma[i] @ gamma[j] + gamma[j] @ gamma[i], target, atol=0)
        expected = gamma[i] if i == 0 else -gamma[i]
        np.testing.assert_allclose(gamma[i].conjugate().T, expected, atol=0)
    # API returns copies, so callers cannot mutate the engine's convention.
    gamma[0][0, 0] = 7
    assert gamma_matrices()[0][0, 0] == 1


@pytest.mark.parametrize("velocity", [[0, 0, 0], [.23, -.31, .17], [.2, .1, .65]])
def test_spinor_dirac_equations_and_completeness(velocity):
    kin = _point(.73, .42).boosted(velocity)
    ident = np.eye(4)
    for p in (kin.p_minus, kin.p_plus):
        for spinor_function, sign in ((u_spinor, 1), (v_spinor, -1)):
            summed = np.zeros((4, 4), dtype=np.clongdouble)
            for spin in (0, 1):
                psi = spinor_function(p, spin)
                np.testing.assert_allclose((slash(p) - sign * M * ident) @ psi, 0, atol=2e-17)
                np.testing.assert_allclose(dirac_adjoint(psi) @ psi, sign * 2 * M, atol=2e-17)
                np.testing.assert_allclose(psi.conjugate() @ psi, 2 * p[0], atol=2e-17)
                summed += np.outer(psi, dirac_adjoint(psi))
            np.testing.assert_allclose(summed, slash(p) + sign * M * ident, rtol=2e-13, atol=2e-17)


@pytest.mark.parametrize("c", [-1, -.37, 0, .8, 1])
def test_on_shell_conserving_cm_and_transverse_basis(c):
    kin = _point(.67, c, phi=.29)
    np.testing.assert_allclose(kin.p_minus + kin.p_plus, kin.k1 + kin.k2, atol=2e-18)
    for p in (kin.p_minus, kin.p_plus):
        assert minkowski_dot(p, p) == pytest.approx(M * M, rel=1e-12)
    for k in (kin.k1, kin.k2):
        assert abs(minkowski_dot(k, k)) < 2e-20
        assert k[0] > 0
        basis = transverse_polarizations(k)
        for i, e in enumerate(basis):
            assert abs(minkowski_dot(e, k)) < 2e-17
            for j, other in enumerate(basis):
                assert minkowski_dot(e, other) == pytest.approx(-1 if i == j else 0, abs=2e-15)


@pytest.mark.parametrize("beta,c", [(.01, .17), (.6, -.8), (.999, .42)])
def test_two_diagrams_ward_cancellation_and_bose_exchange(beta, c):
    kin = _point(beta, c, phi=.41)
    e1s, e2s = transverse_polarizations(kin.k1), transverse_polarizations(kin.k2)
    scale = math.sqrt(unpolarized_amplitude_squared(kin))
    single_diagram = 0
    for sm in (0, 1):
        for sp in (0, 1):
            for e1 in e1s:
                for e2 in e2s:
                    args = (kin.p_minus, kin.p_plus, kin.k1, kin.k2)
                    normal = diagram_amplitudes(*args, e1, e2, sm, sp)
                    for a, b in ((kin.k1 / kin.k1[0], e2), (e1, kin.k2 / kin.k2[0])):
                        ward = diagram_amplitudes(*args, a, b, sm, sp)
                        assert abs(ward["total"]) / scale < 1e-10
                        single_diagram = max(single_diagram, abs(ward["t"]), abs(ward["u"]))
                    swapped = diagram_amplitudes(kin.p_minus, kin.p_plus, kin.k2, kin.k1,
                                                  e2, e1, sm, sp)
                    assert abs(swapped["total"] - normal["total"]) / scale < 1e-12
                    assert abs(swapped["t"] - normal["u"]) / scale < 1e-12
                    assert abs(swapped["u"] - normal["t"]) / scale < 1e-12
    # A Ward identity with each exchange diagram separately would be wrong.
    assert single_diagram > 1e-4


def test_gauge_shifts_and_transverse_basis_rotation():
    kin = _point(.9, .51, phi=.8)
    a, b = transverse_polarizations(kin.k1)
    c, d = transverse_polarizations(kin.k2)
    angle1, angle2 = .31, .59
    rotated1 = (math.cos(angle1) * a + math.sin(angle1) * b,
                -math.sin(angle1) * a + math.cos(angle1) * b)
    rotated2 = (math.cos(angle2) * c + math.sin(angle2) * d,
                -math.sin(angle2) * c + math.cos(angle2) * d)
    shifted1 = (rotated1[0] + 1.4 * kin.k1 / kin.k1[0],
                rotated1[1] - .9 * kin.k1 / kin.k1[0])
    shifted2 = (rotated2[0] - 1.2 * kin.k2 / kin.k2[0],
                rotated2[1] + .5 * kin.k2 / kin.k2[0])
    physical = unpolarized_amplitude_squared(kin)
    shifted = unpolarized_amplitude_squared(kin, polarizations1=shifted1, polarizations2=shifted2)
    assert shifted == pytest.approx(physical, rel=2e-12)


@pytest.mark.parametrize("velocity", [[.3, .17, -.1], [-.5, .2, .4], [.05, .1, .85]])
def test_unpolarized_probability_is_lorentz_covariant(velocity):
    kin = _point(.98, -.77, phi=.14)
    bases1, bases2 = transverse_polarizations(kin.k1), transverse_polarizations(kin.k2)
    original = unpolarized_amplitude_squared(kin)
    transformed = kin.boosted(velocity)
    covariant = unpolarized_amplitude_squared(
        transformed, polarizations1=tuple(boost(e, velocity) for e in bases1),
        polarizations2=tuple(boost(e, velocity) for e in bases2))
    fresh = unpolarized_amplitude_squared(transformed)
    assert covariant == pytest.approx(original, rel=2e-11)
    assert fresh == pytest.approx(original, rel=2e-11)
    np.testing.assert_allclose(boost(boost(kin.k1, velocity), -np.array(velocity)),
                               kin.k1, rtol=2e-12, atol=2e-17)


def test_named_amplitudes_spin_average_and_identical_particle_normalization():
    report = benchmark(.01, .37)
    data = report["amplitude_sample"]["amplitudes"]
    assert len(data) == 16
    summed = 0
    for row in data:
        amp = complex(row["total"]["real"], row["total"]["imag"])
        t = complex(row["t"]["real"], row["t"]["imag"])
        u = complex(row["u"]["real"], row["u"]["imag"])
        assert amp == pytest.approx(t + u, rel=1e-12, abs=1e-15)
        summed += abs(amp) ** 2
    cross = report["differential"]
    assert summed / 4 == pytest.approx(cross["msq"], rel=1e-12)
    assert cross["event_dsigma_domega_gev2"] == pytest.approx(
        cross["msq"] / (128 * math.pi ** 2 * .01 ** 2 * cross["beta"]), rel=1e-13)
    assert cross["event_dsigma_dcostheta_gev2"] == pytest.approx(
        2 * math.pi * cross["event_dsigma_domega_gev2"], rel=1e-13)
    assert cross["event_dsigma_dcostheta_nb"] == pytest.approx(
        GEV2_TO_NB * cross["event_dsigma_dcostheta_gev2"], rel=1e-13)


def test_threshold_sigma_velocity_and_high_energy_central_limit():
    beta = .001
    energy = 2 * M / math.sqrt(1 - beta * beta)
    total = integrate_cross_section(energy, order=64)
    # Fixed-order massive QED: lim(beta->0) sigma*2beta = pi*alpha^2/m^2.
    target = math.pi * ALPHA_REFERENCE ** 2 / M ** 2
    assert total["sigma_times_2beta_gev2"] == pytest.approx(target, rel=3e-6)
    c = .37
    central = differential_cross_section(1.0, c)["msq"]
    e2 = 4 * math.pi * ALPHA_REFERENCE
    massless_target = 4 * e2 ** 2 * (1 + c * c) / (1 - c * c)
    assert central == pytest.approx(massless_target, rel=3e-6)


def test_massive_forward_peak_retained_and_rapidity_integration_converges():
    beta = .9999
    energy = 2 * M / math.sqrt(1 - beta * beta)
    centre = differential_cross_section(energy, 0)["msq"]
    forward = differential_cross_section(energy, 1)["msq"]
    backward = differential_cross_section(energy, -1)["msq"]
    assert math.isfinite(forward) and forward > 1000 * centre
    assert forward == pytest.approx(backward, rel=1e-13)
    order48 = integrate_cross_section(energy, order=48)
    order96 = integrate_cross_section(energy, order=96)
    assert order48["event_sigma_gev2"] == pytest.approx(order96["event_sigma_gev2"], rel=1e-10)
    cut = integrate_cross_section(energy, order=64, cos_theta_max=.9)
    assert cut["event_sigma_gev2"] < order96["event_sigma_gev2"]


def test_coupling_mass_and_azimuth_scaling():
    one = differential_cross_section(.01, .37, phi=.0)
    rotated = differential_cross_section(.01, .37, phi=.9)
    doubled_charge_squared = differential_cross_section(.01, .37, alpha=2 * ALPHA_REFERENCE)
    scaled_mass_and_energy = differential_cross_section(.02, .37, mass_gev=2 * M)
    assert rotated["msq"] == pytest.approx(one["msq"], rel=1e-12)
    assert doubled_charge_squared["msq"] == pytest.approx(4 * one["msq"], rel=1e-12)
    assert scaled_mass_and_energy["msq"] == pytest.approx(one["msq"], rel=1e-12)
    assert scaled_mass_and_energy["event_dsigma_dcostheta_gev2"] == pytest.approx(
        one["event_dsigma_dcostheta_gev2"] / 4, rel=1e-12)


@pytest.mark.parametrize("energy,c,mass", [
    (2 * M, 0, M), (M, 0, M), (-.01, 0, M), (float("nan"), 0, M),
    (float("inf"), 0, M), (.01, 1.01, M), (.01, float("nan"), M),
    (.01, 0, 0), (.01, 0, -M), (.01, 0, float("inf")),
])
def test_invalid_kinematics_rejected(energy, c, mass):
    with pytest.raises(ValueError):
        centre_of_mass(energy, c, mass_gev=mass)


@pytest.mark.parametrize("alpha", [0, -1, float("nan"), float("inf"), True, ".01", .1 + .01j])
def test_invalid_coupling_rejected(alpha):
    with pytest.raises(ValueError):
        differential_cross_section(.01, .37, alpha=alpha)


def test_generalized_input_contracts_reject_nonphysical_data():
    kin = _point(.6, .37)
    with pytest.raises(ValueError):
        Kinematics(kin.p_minus, kin.p_plus, kin.k1 * 1.01, kin.k2)
    with pytest.raises(ValueError):
        Kinematics(kin.p_minus * 1.1, kin.p_plus, kin.k1, kin.k2)
    with pytest.raises(ValueError):
        boost(kin.p_minus, [1, 0, 0])
    with pytest.raises(ValueError):
        boost(kin.p_minus, [0, float("nan"), 0])
    with pytest.raises(ValueError):
        boost(kin.p_minus, [0, .2 + .1j, 0])
    with pytest.raises(ValueError):
        slash([1 + 0j, 0, 0, 0])
    with pytest.raises(ValueError):
        slash(["1", "0", "0", "0"])
    with pytest.raises(ValueError):
        centre_of_mass(.01, ".37")
    with pytest.raises(ValueError):
        centre_of_mass(.01, .37, phi=True)
    with pytest.raises(ValueError):
        u_spinor(kin.p_minus, 0.0)
    with pytest.raises(ValueError):
        v_spinor(kin.p_plus, True)
    bad_basis = (np.array([0, 0, 0, 0]), np.array([0, 1, 0, 0]))
    with pytest.raises(ValueError):
        unpolarized_amplitude_squared(kin, polarizations1=bad_basis)
    for order in (0, 1, 7, 8.5, True):
        with pytest.raises(ValueError):
            integrate_cross_section(.01, order=order)
    for cut in (0, -1, 1.1, float("nan")):
        with pytest.raises(ValueError):
            integrate_cross_section(.01, cos_theta_max=cut)
    # The code must not quietly round a massive process to beta=1.
    with pytest.raises(ValueError):
        centre_of_mass(1e12, 1)


def test_terminal_benchmark_json_is_deterministic(tmp_path):
    repository = Path(__file__).resolve().parents[1]
    environment = dict(os.environ)
    existing = environment.get("PYTHONPATH", "")
    environment["PYTHONPATH"] = str(repository / "src") + (os.pathsep + existing if existing else "")
    output = tmp_path / "benchmark.json"
    cmd = [sys.executable, "-m", "mtft.matter", "benchmark", "--sqrt-s-gev", ".01",
           "--cos-theta", ".37", "--output", str(output)]
    subprocess.run(cmd, check=True, cwd=repository, env=environment, capture_output=True, text=True)
    first = output.read_bytes()
    subprocess.run(cmd, check=True, cwd=repository, env=environment, capture_output=True, text=True)
    assert output.read_bytes() == first
    result = json.loads(first)
    assert result["schema"] == "mtft.matter01.qed.benchmark.v1"
    assert result["differential"]["msq"] > 0


def test_numerical_gate_report_diagnostic_not_certification():
    report = numerical_gates()
    assert report["status"] == "DIAGNOSTIC"
    assert report["passed"]
    assert report["count"] == 12
    assert all(gate["passed"] for gate in report["gates"].values())


def test_custom_parameter_provenance_is_supplied_not_pdg():
    result = benchmark(.01, .37, mass_gev=.001, alpha=.01, include_amplitudes=False)
    inputs = result["inputs"]
    assert inputs["electron_mass_gev"] == .001
    assert inputs["alpha"] == .01
    for name in ("electron_mass_gev", "alpha"):
        assert inputs["parameter_provenance"][name]["origin"] == "USER OVERRIDE"
        assert inputs["parameter_provenance"][name]["source_url"] is None
    assert inputs["reference_defaults"]["electron_mass_gev"] == M
    assert "source_url" not in inputs
    # Overriding one parameter must not detach provenance from the other.
    mixed = benchmark(.01, .37, alpha=.01, include_amplitudes=False)["inputs"]
    assert mixed["parameter_provenance"]["electron_mass_gev"]["origin"] == "PDG 2026 DEFAULT"
    assert mixed["parameter_provenance"]["electron_mass_gev"]["source_url"].startswith("https://pdg.")
    assert mixed["parameter_provenance"]["alpha"]["origin"] == "USER OVERRIDE"
