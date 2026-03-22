"""
MTFT — Modular Time Field Theory  (v0.3.0)
============================================

16 modules covering the complete MTFT framework from arithmetic
weights through lattice gauge theory and experimental verification.
"""

__version__ = "0.4.1"

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

# ── Tier 2: Gauge-Higgs ──────────────────────────────────────
from mtft.hosotani import HosotaniPotential
from mtft.particles import Particle, StandardModel

# ── Tier 3: Phenomenology ────────────────────────────────────
from mtft.dark_sector import TauVortexHalo, rotation_curve, tully_fisher
from mtft.info_geometry import fisher_rao_metric, ricci_scalar_logistic
from mtft.cosmology import FriedmannMTFT

# ── Tier 4: Lattice & Fermions ───────────────────────────────
from mtft.lattice import LatticeConfig, MTFTAction, metropolis_sweep, avg_plaquette, avg_polyakov
from mtft.x0_143 import tano_mass_predictions, koide_angle_tano, generation_count, JACOBIAN
from mtft.koide import koide_ratio, koide_leptons, predict_tau_mass, koide_manifold_point
from mtft.burning_ship import burning_ship_iterate, ANISOTROPIC

# ── Tier 5: Dimensional Bridge & Decay ───────────────────────
from mtft.dimensional_bridge import electron_mass_from_eta, charge_from_feigenbaum
from mtft.decay import ModularDecay

# ── Tier 6: Quantum Computing ────────────────────────────────
from mtft.quantum import (
    TopologicalQudit, HolonomyGate, ArithmeticCode,
    gell_mann_matrices, topological_spectrum_info, skyrmion_number,
)
