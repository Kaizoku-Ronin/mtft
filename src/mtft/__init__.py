"""
MTFT — Modular Time Field Theory  (v0.7.0)
============================================

33 modules covering the complete MTFT framework from arithmetic
weights through lattice gauge theory, materials science, LHC
confrontation, computation theory, and sonification.

Install: pip install mtft
CLI:     python -m mtft verify | report | tower | screen | info
"""

__version__ = "0.7.1"

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
