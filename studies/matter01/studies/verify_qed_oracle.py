#!/usr/bin/env python3
"""Compare independent scalar formulas with the explicit-spinor QED engine.

Predeclared finite-grid numerical diagnostics. No engine spinor, gamma, boost,
polarization, Mandelstam, or kinematics-construction helper is used here.
The engine's Kinematics container and public scattering APIs are the subject
under test. Independent external momenta/boosts and analytic scalar formulas
are constructed locally. This is not an interval certificate or data fit.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

STUDY_ROOT = Path(__file__).resolve().parents[1]
REPO = STUDY_ROOT / "repo_candidate" if (STUDY_ROOT / "repo_candidate").exists() else STUDY_ROOT.parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np
from scipy.integrate import quad
from mtft.matter import oracle

PHI_GRID = (0.0, 0.61, 1.73)
BOOST = (0.12, -0.07, 0.2)


def relative_error(actual, expected):
    return abs(float(actual) - float(expected)) / max(abs(float(expected)), 1e-300)


def independently_constructed_momenta(beta, cosine, phi=0.0, mass=oracle.MASS_GEV):
    energy = mass / math.sqrt(1.0 - beta*beta)
    momentum = energy * beta
    transverse = math.sqrt(1.0-cosine*cosine)
    k = energy * np.asarray((transverse*math.cos(phi), transverse*math.sin(phi), cosine))
    return (
        np.asarray((energy, 0.0, 0.0, momentum)),
        np.asarray((energy, 0.0, 0.0, -momentum)),
        np.concatenate(([energy], k)),
        np.concatenate(([energy], -k)),
    )


def independently_boosted(vector, velocity=BOOST):
    velocity = np.asarray(velocity, dtype=float)
    speed2 = float(velocity @ velocity)
    gamma = 1.0 / math.sqrt(1.0-speed2)
    projection = float(velocity @ vector[1:])
    output = np.empty(4)
    output[0] = gamma*(vector[0]+projection)
    output[1:] = vector[1:] + velocity*((gamma-1.0)*projection/speed2+gamma*vector[0])
    return output


def run(include_engine=True):
    if include_engine:
        from mtft.matter import qed
    records, gates = [], []

    def gate(name, measured, limit, samples):
        gates.append({"name": name, "maximum_relative_error": measured,
                      "limit": limit, "samples": samples, "passed": measured <= limit})

    maxima = {name: 0.0 for name in ("cm_vs_invariants", "cm_vs_heitler", "quad_vs_total",
              "engine_msq", "engine_boosted_msq", "engine_domega", "engine_dcos",
              "engine_units", "engine_total", "engine_integration_convergence")}
    for beta in oracle.BETA_GRID:
        s = 4.0*oracle.MASS_GEV**2/(1.0-beta**2)
        expected_total = oracle.sigma_total(beta)
        heitler = oracle.sigma_total_heitler(beta)
        integrated, quadrature_error = quad(lambda c: oracle.dsigma_dcosine(beta, c),
                                            -1.0, 1.0, epsabs=1e-9, epsrel=3e-13,
                                            limit=300)
        maxima["cm_vs_heitler"] = max(maxima["cm_vs_heitler"], relative_error(heitler, expected_total))
        maxima["quad_vs_total"] = max(maxima["quad_vs_total"], relative_error(integrated, expected_total))
        total_record = {"beta": beta, "sqrt_s_gev": math.sqrt(s),
                        "event_sigma_gev2": expected_total,
                        "event_sigma_nb": expected_total*oracle.GEV2_TO_NB,
                        "heitler_sigma_gev2": heitler,
                        "independent_quad_sigma_gev2": integrated,
                        "quad_reported_absolute_error": quadrature_error,
                        "sigma_times_2beta_gev2": oracle.sigma_v_moller(beta)}
        if include_engine:
            numerical64 = qed.integrate_cross_section(math.sqrt(s), mass_gev=oracle.MASS_GEV,
                                                       alpha=oracle.ALPHA, order=64)
            numerical128 = qed.integrate_cross_section(math.sqrt(s), mass_gev=oracle.MASS_GEV,
                                                        alpha=oracle.ALPHA, order=128)
            maxima["engine_total"] = max(maxima["engine_total"],
                relative_error(numerical128["event_sigma_gev2"], expected_total))
            maxima["engine_integration_convergence"] = max(maxima["engine_integration_convergence"],
                relative_error(numerical64["event_sigma_gev2"], numerical128["event_sigma_gev2"]))
            total_record["engine64"] = numerical64
            total_record["engine128"] = numerical128
        records.append(total_record)

        for cosine in oracle.COSINE_GRID:
            stu = oracle.mandelstam(beta, cosine)
            expected_msq = oracle.msq_cm(beta, cosine)
            invariant_msq = oracle.msq_mandelstam(*stu)
            maxima["cm_vs_invariants"] = max(maxima["cm_vs_invariants"], relative_error(invariant_msq, expected_msq))
            if not include_engine:
                continue
            for phi in PHI_GRID:
                momenta = independently_constructed_momenta(beta, cosine, phi)
                kin = qed.Kinematics(*momenta, mass_gev=oracle.MASS_GEV)
                actual_msq = qed.unpolarized_amplitude_squared(kin, alpha=oracle.ALPHA)
                maxima["engine_msq"] = max(maxima["engine_msq"], relative_error(actual_msq, expected_msq))
                differential = qed.differential_cross_section(math.sqrt(s), cosine, phi=phi,
                                          mass_gev=oracle.MASS_GEV, alpha=oracle.ALPHA)
                maxima["engine_domega"] = max(maxima["engine_domega"], relative_error(
                    differential["event_dsigma_domega_gev2"], oracle.dsigma_domega(beta, cosine)))
                maxima["engine_dcos"] = max(maxima["engine_dcos"], relative_error(
                    differential["event_dsigma_dcostheta_gev2"], oracle.dsigma_dcosine(beta, cosine)))
                maxima["engine_units"] = max(maxima["engine_units"], relative_error(
                    differential["event_dsigma_dcostheta_nb"], oracle.dsigma_dcosine(beta, cosine)*oracle.GEV2_TO_NB))
            boosted = qed.Kinematics(*(independently_boosted(p) for p in momenta), mass_gev=oracle.MASS_GEV)
            actual_boosted_msq = qed.unpolarized_amplitude_squared(boosted, alpha=oracle.ALPHA)
            maxima["engine_boosted_msq"] = max(maxima["engine_boosted_msq"],
                                               relative_error(actual_boosted_msq, expected_msq))

    gate("CM scalar versus invariant scalar", maxima["cm_vs_invariants"], 3e-12, 42)
    gate("CM total versus independent Heitler total", maxima["cm_vs_heitler"], 3e-12, 6)
    gate("adaptive direct-angle integration versus analytic total", maxima["quad_vs_total"], 3e-12, 6)
    if include_engine:
        gate("explicit spinor Msq versus independent oracle", maxima["engine_msq"], 2e-10, 126)
        gate("independently boosted spinor Msq versus invariant oracle", maxima["engine_boosted_msq"], 2e-10, 42)
        gate("engine event differential per solid angle", maxima["engine_domega"], 2e-10, 126)
        gate("engine event differential per cosine", maxima["engine_dcos"], 2e-10, 126)
        gate("engine nb conversion", maxima["engine_units"], 2e-10, 126)
        gate("engine integrated event total versus oracle", maxima["engine_total"], 2e-10, 6)
        gate("engine integration orders 64 versus 128", maxima["engine_integration_convergence"], 2e-10, 6)

    threshold_reference = math.pi*oracle.ALPHA**2/oracle.MASS_GEV**2
    threshold = [{"beta": beta, "ratio_sigma_times_2beta_to_threshold":
                  oracle.sigma_v_moller(beta)/threshold_reference} for beta in (1e-3, 1e-5, 1e-7)]
    gate("formal Born sigma*v threshold limit", abs(threshold[-1]["ratio_sigma_times_2beta_to_threshold"]-1), 1e-12, 3)
    high_energy = [{"beta": beta, "ratio_exact_to_asymptotic":
                   oracle.sigma_total(beta)/oracle.sigma_high_energy(beta)} for beta in (0.999, 0.9999, 0.99999)]
    gate("high-energy asymptotic total", abs(high_energy[-1]["ratio_exact_to_asymptotic"]-1), 1e-4, 3)
    scaling_error = max(relative_error(oracle.sigma_total(beta, mass=2*oracle.MASS_GEV)*4,
                                      oracle.sigma_total(beta)) for beta in oracle.BETA_GRID)
    gate("independent mass dimension scaling", scaling_error, 1e-13, 6)
    alpha_scaling_error = max(relative_error(oracle.msq_cm(beta, c, 2*oracle.ALPHA)/4,
                                             oracle.msq_cm(beta, c)) for beta in oracle.BETA_GRID for c in oracle.COSINE_GRID)
    gate("independent alpha squared scaling", alpha_scaling_error, 1e-13, 42)

    return {"study": "Matter-01 independent massive QED oracle",
            "status": "FINITE-GRID NUMERICAL DIAGNOSTIC; ordinary Born QED reference",
            "engine_included": include_engine,
            "protocol": {"beta_grid": oracle.BETA_GRID, "cosine_grid": oracle.COSINE_GRID,
                         "phi_grid": PHI_GRID, "independent_boost": BOOST,
                         "event_identical_photon_factor": 0.5,
                         "engine_quadrature_orders": [64, 128]},
            "inputs": {"electron_mass_gev": oracle.MASS_GEV, "alpha_inverse": 1/oracle.ALPHA,
                       "gev2_to_nb": oracle.GEV2_TO_NB, "provenance": "PDG 2026 supplied inputs"},
            "gates": gates, "all_passed": all(g["passed"] for g in gates),
            "total_grid": records, "formal_threshold_checks": threshold,
            "angular_grid": [{"beta": beta, "energy_ratio": 1/math.sqrt(1-beta**2),
                "cos_theta": cosine, "dsigma_dcos_nb": oracle.dsigma_dcosine(beta, cosine)*oracle.GEV2_TO_NB}
                for beta in (0.1, 0.4, 0.7, 0.9, 0.98) for cosine in (-1.0, -0.5, 0.0, 0.5, 1.0)],
            "high_energy_checks": high_energy,
            "scope": "No parameter derivation, positronium, Coulomb resummation, radiative corrections, material effects, or MTFT vacuum validation."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-engine", action="store_true")
    parser.add_argument("--output", type=Path, default=STUDY_ROOT / "results" / "qed_oracle_checks.json")
    args = parser.parse_args()
    report = run(include_engine=not args.skip_engine)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"all_passed": report["all_passed"], "gates": report["gates"],
                      "output": str(args.output)}, indent=2))
    raise SystemExit(0 if report["all_passed"] else 1)
