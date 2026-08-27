"""
MTFT — Modular Time Field Theory  (v0.15.0 — certification wave)
============================================================================

39 modules covering the complete MTFT framework from arithmetic
weights through lattice gauge theory, materials science, LHC
confrontation, computation theory, and sonification — plus the
primon-gas spectral reconstruction toolkit (ledger, chain, expansion,
ep; Integration Plan v0.1 stages 1–3, Add. AZ–BL).  coupled.py
(stage 4) and the studies re-pointing (stage 5) are pending; the
All five Integration Plan stages are landed; studies/ holds the re-pointed suite.

Install: pip install mtft
CLI:     python -m mtft verify | report | tower | screen | info
"""

__version__ = "0.22.0"

# ── Tier 0: Constants & Arithmetic ───────────────────────────
from mtft.constants import (
    SM, GAUGE, HIGGS, PDG, PhysicalConstants, CriticalDepths, MassRatios,
    FEIGENBAUM_DELTA, FEIGENBAUM_ALPHA, T_INF, LAMBERT_OMEGA, EULER_GAMMA,
    QUARKS, LEPTONS, ZETA_2, XI,
)
from mtft.arithmetic import (
    weight, weight_array, damped_weights,
    S1, stiffness_S, center_stiffness, mass_gap_stiffness,
    su3_hessian_eigenvalues, find_confinement_depth,
)

# ── Tier 0b: Combinatorial Ancestry (2024 lineage, v0.13.0) ──
from mtft import combinatorial  # noqa: F401
from mtft.combinatorial import (
    figurate_decomposition, sigma_parity_check, mean_tano_weight_exact,
)

# ── Tier 1: Modular Geometry & Forms ─────────────────────────
from mtft.modular import TauField, sl2z_transform
from mtft.forms import dedekind_eta, jacobi_theta3, eisenstein_E2k, verify_spectral_identity
from mtft.modular_curve import ModularCurve, X0

# ── Tier 2: Gauge-Higgs ──────────────────────────────────────
from mtft.hosotani import HosotaniPotential, HosotaniMTFT
from mtft.particles import Particle, StandardModel

# ── Tier 3: Phenomenology ────────────────────────────────────
from mtft.dark_sector import TauVortexHalo, rotation_curve, rotation_curve_kpc, tully_fisher
from mtft.info_geometry import fisher_rao_metric, ricci_scalar_logistic
from mtft.cosmology import FriedmannMTFT

# ── Tier 4: Lattice & Fermions ───────────────────────────────
from mtft.lattice import LatticeConfig, MTFTAction, metropolis_sweep, avg_plaquette, avg_polyakov
from mtft.x0_143 import (
    tano_mass_predictions, koide_angle_tano, generation_count, JACOBIAN,
    ORBIT_TRACES_VERIFIED, TRACE_TOTALS_50, ROOT_NUMBERS_LIST,
    rankin_selberg_Q, verify_complex_eigenvalue,
)
from mtft.koide import koide_ratio, koide_leptons, predict_tau_mass, koide_manifold_point
from mtft.burning_ship import burning_ship_iterate, ANISOTROPIC

# ── Tier 5: Dimensional Bridge & Decay ───────────────────────
from mtft.dimensional_bridge import electron_mass_from_eta, charge_from_feigenbaum
from mtft.decay import ModularDecay

# ── Tier 5b: Falsifiability Engine ───────────────────────────
from mtft.falsify import (
    prediction_table, coupling_shift, coupling_shift_table,
    falsification_test, holonomy_flux, desert_check,
)

# ── Tier 5c: Multi-N Tower ───────────────────────────────────
from mtft.tower import (
    tower_stiffness, even_n_universality, phase_transition_scaling,
    boundary_tracking, arithmetic_genome, arithmetic_periodic_table,
    character_orthogonality, tower_report,
)

# ── Tier 5d: Riemann Explicit Formula ────────────────────────
from mtft.riemann import (
    RIEMANN_ZEROS, explicit_formula, bakry_emery_curvature,
    rh_diagnostic, tower_rigidity, gamma_suppression_table,
    skeleton_stiffness, normalized_oscillation, envelope_slope,
    corrected_rh_diagnostic,
    ZETAPRIME_ZEROS, ZETAPRIME_CENSUS_HEIGHT,
    zetaprime_zero_count_berndt, zetaprime_negative_zero,
    zetaprime_refine, zetaprime_logcurvature, hadamard_zetaprime_check,
    dirichlet_curvature, von_mangoldt_curvature, divisor_log_weights,
    weighted_theta, filtered_moment_identity, weighted_theta_cusp_fit,
)
from mtft.marked_gas import (  # the marked primon gas (note v0.1.1)
    ALPHA_COLD, B_COLD, Certified,
    z1, z2, zD_certified_interval,
    spectrum, flow_phase, kms_check, bc_deformation,
    weights_sieve, psi_coefficients, cold_gas_report,
    correlator, spectral_function, edge_mass,
    gates as marked_gas_gates,  # aliased: bare "gates" is too generic flat
)

# ── Tier 5e: Materials Science ───────────────────────────────
from mtft.tano_metric import (
    tano_contrast, geometry_index, predict_Tc, materials_screening,
    seebeck_diagnostic, josephson_holonomy, ELEMENTS,
)

# ── Tier 6: Quantum Computing ────────────────────────────────
from mtft.quantum import (
    TopologicalQudit, HolonomyGate, ArithmeticCode,
    gell_mann_matrices, topological_spectrum_info, skyrmion_number,
)

# ── Tier 7: Cryptography ────────────────────────────────────
from mtft.monster_hash import MonsterHash
from mtft import crypto  # noqa: F401  (package: primitives + jacobian_order)

# ── Tier 8: LHC Confrontation ───────────────────────────────
# (optional — requires uproot)
try:
    from mtft.lhcb_analysis import LHCbNtuple
except ImportError:
    pass  # uproot not installed

# ── Tier 9: Computation, Jacobian Engine & Sonification ──────
from mtft.arithmetic_machine import (
    Primitive, PrimitiveLevel, classify_computation, computational_stiffness,
    search_space_compression, arithmetic_entropy, analyze_halting_surface,
)
from mtft.arithmetic_wick import (
    dirichlet_ensemble, laplace_ensemble, wick_rotate, wick_at_critical_depths,
)
from mtft.busy_beaver import (
    HeckeSign, hecke_sign, bb_genus, bb_fatou, bb_sample,
    faulhaber_decompose, verify_telescoping,
)
from mtft.jacobian import JacobianStiffness
from mtft.music import VacuumSonifier, ModularScale, MonsterComposer, spectral_fingerprint

# ── Tier 10: Critical Ensemble (Li coefficients) ─────────────
from mtft.critical_ensemble import (
    lambda_1_closed_form, li_lambda, li_lambda_batch,
    logxi_taylor, logxi_taylor_cauchy, li_lambda_cauchy,
    li_lambda_zero_sum, certify, li_criterion_report,
    BOMBIERI_LAGARIAS_CAVEAT, XI_ANALYTICITY_RADIUS,
    THREE_ENSEMBLE_TABLE,
)

# ── Tier 11: Certificates & Standards ────────────────────────
from mtft.jc_counterexample import (
    verify_all as jc_verify_all, JCCertificate, apply_F as jc_apply,
    COLLISION_FIBER, COLLISION_TARGET,
)
from mtft.estimator_standards import (
    binned_log_slope, stride_resonance_check, recommended_samples_per_decade,
)

# ── Tier 5f: Primon-Gas Spectral Reconstruction (v0.11.0, stages 1–5) ──
# The rung-4/rung-5 spectral toolkit: the T = DKD kernel chain,
# exceptional points, the remainder expansion, the coupled (Bloch /
# Kesten) model, and the certified ledger of constants (Integration
# Plan v0.1 §3; audited Add. BI/BK/BN, dispositioned Add. BL).
# Exposed AS MODULES, not flattened — the function names (gap,
# newton, census) are too generic for this namespace.
# expansion.richardson is deliberately NOT part of the public surface
# pending BI.F1 (the audited Neville fix landed in chain(1)/expansion(1)
# and is verified, Add. BN §4; the re-export decision is the author's).
# coupled.selftest now asserts through the _L guard (BN-F1 closed; the
# quadrature cache is single-eigh, BN-F2 closed).  Historical note: it
# previously asserted two literals pending its switch to the
# _L guard (BN-F1; the numbers are registered as tau_c_star / V_b_tree
# in the ledger, Add. BN §7).
from mtft import ledger, chain, expansion, ep, coupled  # noqa: F401
from mtft.ledger import verify as spectral_ledger_verify
from mtft.chain import selftest as spectral_chain_selftest
from mtft.expansion import selftest as spectral_expansion_selftest
from mtft.ep import selftest as spectral_ep_selftest
from mtft.coupled import selftest as spectral_coupled_selftest

# ── Tier 12: Promotion Wave (v0.14.0) ─────────────────────────
# moments (Tano weight closed forms), curvature (Brioschi machinery
# of the (beta, lambda) manifold), hecke (Manin/Merel engine for
# X0(143)), eisenstein (congruence moduli).  Exposed AS MODULES
# (tier-5f precedent): the names are generic, the engines exact.
from mtft import moments, curvature, hecke, eisenstein  # noqa: F401

# ── Tier 13: v0.15.0 wave ─────────────────────────────────────
# weil (Gabor-compressed Weil explicit-formula form, W1): prime-side
# vs zero-side identity of the compressed matrix G, E2-certified at
# 3.472e-09 by the independent Kimi implementation.  Module-level
# import, same tier-5f precedent.
from mtft import weil  # noqa: F401

# ── The Legend (lazy: keeps `python -m mtft.legend` warning-free) ──
def __getattr__(name):
    _LEGEND_EXPORTS = {"legend": "legend", "what": "what", "card": "card",
                       "trace": "trace", "legend_status": "status"}
    if name in _LEGEND_EXPORTS:
        import importlib
        return getattr(importlib.import_module("mtft.legend"),
                       _LEGEND_EXPORTS[name])
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

# ── Tier 6: Origami / dimers / insertion calculus (v0.20.0) ──
from mtft import origami  # noqa: F401
from mtft.origami import (
    DimerGraph, ensemble_conservation, fisher_metric, cubic_tensor,
    brioschi_curvature, cumulant_curvature, path_independence,
    Theta, solve_perfect_branches, orbit_structure,
    galashin_24, prism_36, PRISM_C, PRISM_LAMBDA0,
)
from mtft import hardy_ramanujan  # noqa: F401
from mtft.hardy_ramanujan import (
    psi_direct, psi_modular, saddle_partition, hardy_ramanujan_asymptotic,
)
