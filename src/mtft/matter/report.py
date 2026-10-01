"""Deterministic numerical diagnostics and a browser-ready QED reference report."""

import math
from typing import Dict

import numpy as np

from .qed import (
    ALPHA_REFERENCE, ELECTRON_MASS_GEV, GEV2_TO_NB, REFERENCE_INPUTS, Kinematics,
    amplitude_sample, boost, centre_of_mass, diagram_amplitudes,
    differential_cross_section, dirac_adjoint, gamma_matrices,
    integrate_cross_section, slash, transverse_polarizations,
    u_spinor, unpolarized_amplitude_squared, v_spinor,
)


def _relative(actual, expected):
    return float(abs(actual - expected) / max(abs(expected), 1e-30))


def numerical_gates() -> Dict[str, object]:
    """Check the explicit-amplitude engine's structural identities.

    These are numerical diagnostics with stated tolerances.  Independent
    analytic cross-section validation is supplied separately in Matter-01's
    oracle study; it is not manufactured from the amplitude implementation.
    """
    gamma = gamma_matrices()
    identity = np.eye(4)
    clifford_error = max(float(np.max(np.abs(gamma[a] @ gamma[b] + gamma[b] @ gamma[a]
                       - (2 * (1 if a == 0 else -1) * identity if a == b else 0))))
                       for a in range(4) for b in range(4))
    kin = centre_of_mass(0.01, 0.37, phi=0.29)
    m = ELECTRON_MASS_GEV
    spinor_error, completeness_error = 0.0, 0.0
    for p in (kin.p_minus, kin.p_plus):
        for kind, sign in ((u_spinor, 1), (v_spinor, -1)):
            accum = np.zeros((4, 4), dtype=np.clongdouble)
            for spin in (0, 1):
                psi = kind(p, spin)
                residual = (slash(p) - sign * m * identity) @ psi
                scale = float(p[0]) * float(np.max(np.abs(psi)))
                spinor_error = max(spinor_error, float(np.max(np.abs(residual))) / scale)
                accum += np.outer(psi, dirac_adjoint(psi))
            expected = slash(p) + sign * m * identity
            completeness_error = max(completeness_error,
                                     float(np.max(np.abs(accum - expected))) / float(p[0]))
    bases = (transverse_polarizations(kin.k1), transverse_polarizations(kin.k2))
    amplitude_scale = math.sqrt(unpolarized_amplitude_squared(kin))
    ward, gauge_shift, bose = 0.0, 0.0, 0.0
    individual_ward = 0.0
    for sm in (0, 1):
        for sp in (0, 1):
            for i in (0, 1):
                for j in (0, 1):
                    eps1, eps2 = bases[0][i], bases[1][j]
                    args = (kin.p_minus, kin.p_plus, kin.k1, kin.k2)
                    original = diagram_amplitudes(*args, eps1, eps2, sm, sp)
                    first = diagram_amplitudes(*args, kin.k1 / kin.k1[0], eps2, sm, sp)
                    second = diagram_amplitudes(*args, eps1, kin.k2 / kin.k2[0], sm, sp)
                    ward = max(ward, abs(first["total"]) / amplitude_scale,
                               abs(second["total"]) / amplitude_scale)
                    individual_ward = max(individual_ward, abs(first["t"]), abs(first["u"]))
                    shifted = diagram_amplitudes(*args, eps1 + 0.83 * kin.k1 / kin.k1[0],
                                                 eps2 - 0.47 * kin.k2 / kin.k2[0], sm, sp)
                    gauge_shift = max(gauge_shift, abs(shifted["total"] - original["total"])
                                      / amplitude_scale)
                    swapped = diagram_amplitudes(kin.p_minus, kin.p_plus, kin.k2, kin.k1,
                                                  eps2, eps1, sm, sp)
                    bose = max(bose, abs(swapped["total"] - original["total"])
                               / amplitude_scale)
    initial_msq = unpolarized_amplitude_squared(kin)
    velocity = np.array([0.23, -0.31, 0.17])
    moved = kin.boosted(velocity)
    covariant_msq = unpolarized_amplitude_squared(
        moved, polarizations1=tuple(boost(e, velocity) for e in bases[0]),
        polarizations2=tuple(boost(e, velocity) for e in bases[1]))
    covariance = _relative(covariant_msq, initial_msq)
    fresh_basis = _relative(unpolarized_amplitude_squared(moved), initial_msq)
    coupling = _relative(unpolarized_amplitude_squared(kin, alpha=2 * ALPHA_REFERENCE),
                         4 * initial_msq)
    order48 = integrate_cross_section(0.01, order=48)["event_sigma_gev2"]
    order96 = integrate_cross_section(0.01, order=96)["event_sigma_gev2"]
    convergence = _relative(order48, order96)
    unit_conversion = _relative(
        differential_cross_section(0.01, 0.37)["event_dsigma_dcostheta_nb"],
        differential_cross_section(0.01, 0.37)["event_dsigma_dcostheta_gev2"] * GEV2_TO_NB)
    values = {
        "clifford_algebra": (clifford_error, 1e-15),
        "dirac_equations": (spinor_error, 1e-12),
        "spinor_completeness": (completeness_error, 1e-12),
        "ward_full_sum": (ward, 1e-11),
        "gauge_shift_invariance": (gauge_shift, 1e-11),
        "bose_exchange": (bose, 1e-12),
        "boost_covariance": (covariance, 1e-10),
        "boost_fresh_transverse_basis": (fresh_basis, 1e-10),
        "alpha_squared_scaling": (coupling, 1e-12),
        "rapidity_quadrature_convergence": (convergence, 2e-10),
        "nb_conversion": (unit_conversion, 1e-14),
    }
    gates = {name: {"measured_error": error, "tolerance": tolerance,
                    "passed": bool(error <= tolerance)}
             for name, (error, tolerance) in values.items()}
    gates["ward_requires_two_diagrams"] = {
        "individual_diagram_amplitude": float(individual_ward),
        "required_minimum": 1e-5,
        "passed": bool(individual_ward > 1e-5 and ward < 1e-11),
    }
    return {"status": "DIAGNOSTIC", "passed": all(g["passed"] for g in gates.values()),
            "count": len(gates), "gates": gates}


def _input_ledger(mass_gev, alpha):
    """Record supplied values separately from the frozen default reference.

    The baseline citation is not evidence for a user-supplied override.  Each
    active parameter carries its own origin and citation, if applicable.
    """
    reference_url = REFERENCE_INPUTS["source_url"]
    mass_is_default = float(mass_gev) == ELECTRON_MASS_GEV
    alpha_is_default = float(alpha) == ALPHA_REFERENCE
    provenance = {}
    for name, value, is_default in (
            ("electron_mass_gev", float(mass_gev), mass_is_default),
            ("alpha", float(alpha), alpha_is_default)):
        provenance[name] = {
            "value": value, "status": "INPUT",
            "origin": "PDG 2026 DEFAULT" if is_default else "USER OVERRIDE",
            "source_url": reference_url if is_default else None,
        }
    provenance["gev_inverse_squared_to_nb"] = {
        "value": GEV2_TO_NB, "status": "REFERENCE",
        "origin": "PDG 2026 ROUNDED TABLE CONVERSION", "source_url": reference_url,
    }
    return {
        "electron_mass_gev": float(mass_gev), "alpha": float(alpha),
        "alpha_inverse": float(1 / alpha), "gev_inverse_squared_to_nb": GEV2_TO_NB,
        "status": "INPUT: see per-parameter provenance; no values are MTFT predictions",
        "parameter_provenance": provenance,
        "reference_defaults": {
            "description": "PDG baseline applies only to parameters labelled PDG 2026 DEFAULT",
            "electron_mass_gev": ELECTRON_MASS_GEV, "alpha": ALPHA_REFERENCE,
            "alpha_inverse": REFERENCE_INPUTS["alpha_inverse"],
            "source_url": reference_url,
        },
    }


def benchmark(sqrt_s_gev=0.01, cos_theta=0.37, mass_gev=ELECTRON_MASS_GEV,
              alpha=ALPHA_REFERENCE, phi=0.0, include_amplitudes=True) -> Dict[str, object]:
    """One named reference process, with inputs, external momenta and amplitudes."""
    kin = centre_of_mass(sqrt_s_gev, cos_theta, phi=phi, mass_gev=mass_gev)
    result = {
        "schema": "mtft.matter01.qed.benchmark.v1",
        "process": "electron + positron -> photon + photon",
        "status": "REFERENCE: tree-level QED; physical parameters are INPUT",
        "inputs": _input_ledger(mass_gev, alpha),
        "kinematics": kin.as_dict(),
        "differential": differential_cross_section(sqrt_s_gev, cos_theta, phi=phi,
                                                   mass_gev=mass_gev, alpha=alpha),
        "conventions": {"metric": "+---", "units": "natural GeV",
                        "spin_average": "1/4", "polarization_sum": "2 x 2 transverse",
                        "identical_photon_factor": "1/2 for full-sphere event cross sections"},
    }
    if include_amplitudes:
        result["amplitude_sample"] = amplitude_sample(kin, alpha=alpha)
    return result


def benchmark_report() -> Dict[str, object]:
    """Deterministic low-energy reference grid suitable for offline visualization."""
    energies = [0.0011, 0.002, 0.01, 0.1]
    curves = []
    for energy in energies:
        curves.append({"sqrt_s_gev": energy, "points": [
            differential_cross_section(energy, float(c))
            for c in np.linspace(-0.999, 0.999, 161)]})
    return {
        "schema": "mtft.matter01.qed.report.v1",
        "process": "electron + positron -> photon + photon",
        "model": "massive tree-level QED reference engine",
        "status": "REFERENCE; numerical checks are DIAGNOSTIC; not MTFT-derived physics",
        "inputs": _input_ledger(ELECTRON_MASS_GEV, ALPHA_REFERENCE),
        "conventions": {
            "metric": "+---", "units": "natural GeV", "spinor_normalization": "2m",
            "initial_spin_average": "1/4", "final_polarization_sum": "two transverse per photon",
            "cross_section": "events; identical-photon 1/2; integrate full cos(theta) [-1,1]",
            "photon_number_histogram": "twice the event angular cross section",
            "diagram_phase": "common physically irrelevant amplitude phase omitted",
        },
        "verification": numerical_gates(),
        "benchmarks": [benchmark(energy, 0.37) for energy in energies],
        "angular_curves": curves,
        "integrated": [integrate_cross_section(energy, order=96) for energy in energies],
        "limitations": [
            "No MTFT vacuum, spectrum, coupling, or mass is derived by this QED kernel.",
            "No radiative corrections, running alpha, weak interactions, or hadrons.",
            "No Fock-state time evolution, real-time particle creation, or many-body matter simulation.",
            "No dynamical metric or quantum-gravity calculation.",
            "Fixed-order threshold limit excludes Coulomb resummation and positronium binding.",
            "No infrared-safe higher-order observable; only the tree-level two-photon final state.",
            "Finite-precision diagnostics; unsupported massive energy ratios are rejected.",
        ],
    }
