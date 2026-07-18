#!/usr/bin/env python3
"""
Arithmetic Wick Rotation: Two Ensembles on the Same Weights
=============================================================

MIT License — Copyright (c) 2026 Roger Tano
See LICENSE file for full terms.

The Arithmetica Generale's holonomy weights w_n = Σ_{d|n} (log d)/d
can be assembled into two different statistical ensembles:

  LAPLACE (Minkowski / stiffness):
    Z_L(y) = Σ w_n e^{-2πyn},       energy E_n = n
    This is the MTFT stiffness ensemble. The mass gap lives here.
    Parameter: modular depth y ∈ (0, ∞)

  DIRICHLET (Euclidean / zeta):
    Z_D(β) = Σ w_n n^{-β},          energy E_n = log n
    This is the Dirichlet series of w_n. The zeta zeros live here.
    Parameter: inverse temperature β ∈ (1, ∞)

The "arithmetical Wick rotation" is the map between these two pictures:

    Laplace ←→ Dirichlet
    e^{-2πyn} ←→ n^{-β}
    E = n      ←→ E = log n
    y          ←→ β

This is implemented by the Mellin transform, which is the bridge
between additive and multiplicative number theory.

Both ensembles produce Fisher-Rao metrics via CURVE (Primitive V):
    g_L(y)  = Var_y(n)       — curvature in the Laplace picture
    g_D(β)  = Var_β(log n)   — curvature in the Dirichlet picture

The critical points of these metrics — where curvature extremizes —
are the "fixed points of the pipeline" in each picture.

MTFT Constants:
    Level N = 143 = 11 × 13
    Genus g = 13
    Confinement depth y_c = 0.18174

Roger Tano — MTFT Research Program — April 2026
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Dict, Optional


# ═══════════════════════════════════════════════════════════════
#  MTFT STRUCTURAL CONSTANTS
# ═══════════════════════════════════════════════════════════════

LEVEL = 143
GENUS = 13
DIM_NEW = 11
from mtft.constants import CriticalDepths

Y_C = CriticalDepths.y_conf  # Confinement depth 0.18174 (canonical, v0.6.1 audit)
Y_S1 = 0.1236      # Skeleton zero 1
Y_S2 = 0.2106      # Skeleton zero 2

# Meissel-Mertens constant
M_MERTENS = 0.2614972128476427


# ═══════════════════════════════════════════════════════════════
#  §1. HOLONOMY WEIGHTS — THE SHARED DATA
# ═══════════════════════════════════════════════════════════════

def sieve_primes(limit: int) -> List[int]:
    """Eratosthenes sieve."""
    s = [True] * (limit + 1)
    s[0] = s[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            for j in range(i*i, limit + 1, i):
                s[j] = False
    return [i for i in range(2, limit + 1) if s[i]]


def compute_weights(n_max: int) -> List[float]:
    """
    Compute holonomy weights w_n = Σ_{d|n} (log d)/d for n = 0..n_max.

    w_n is the master arithmetic function f(n)/n, where
    f(n) = (σ * Λ)(n) = Σ_{d|n} (n/d) log d.

    These weights are the SHARED DATA between both ensembles.
    The same w_n, assembled differently, produce different physics.
    """
    w = [0.0] * (n_max + 1)
    for d in range(2, n_max + 1):
        v = math.log(d) / d
        for m in range(d, n_max + 1, d):
            w[m] += v
    return w


def compute_skeleton_weights(n_max: int) -> List[float]:
    """
    Skeleton weights: w_n^skel = log(p)/p at primes, 0 elsewhere.

    The prime skeleton is the EXTRACT of the full weights — the
    irreducible component that generates the rest via ASSEMBLE.
    """
    primes = set(sieve_primes(n_max))
    w = [0.0] * (n_max + 1)
    for p in primes:
        w[p] = math.log(p) / p
    return w


# ═══════════════════════════════════════════════════════════════
#  §2. DIRICHLET ENSEMBLE (EUCLIDEAN PICTURE)
# ═══════════════════════════════════════════════════════════════

@dataclass
class DirichletEnsemble:
    """
    The Dirichlet (Euclidean) statistical ensemble.

    Z_D(β) = Σ_{n≥1} w_n n^{-β}

    Energy: E_n = log n
    Degeneracy: w_n
    Inverse temperature: β = Re(s) for Dirichlet variable s
    Gibbs measure: P_β(n) = w_n n^{-β} / Z_D(β)

    This is the "Euclidean" picture because:
    - The zeta zeros live in this picture (poles of -ζ'(s)/ζ(s))
    - β plays the role of Euclidean time
    - The partition function IS the Dirichlet series of w_n
    """
    beta: float
    n_max: int
    partition_fn: float        # Z_D(β)
    free_energy: float         # Φ = log Z_D(β)
    mean_energy: float         # ⟨E⟩ = ⟨log n⟩_β
    energy_variance: float     # Var_β(log n) = Fisher-Rao metric
    specific_heat: float       # C = β² Var_β(log n)
    entropy: float             # S = β⟨E⟩ + log Z


def dirichlet_ensemble(beta: float, n_max: int = 1000,
                       weights: Optional[List[float]] = None) -> DirichletEnsemble:
    """
    Compute the Dirichlet ensemble at inverse temperature β.

    The pipeline:
        ITERATE: primes → divisors → f(n)
        DIVIDE:  n/d inside f(n)
        ASSEMBLE: Z_D(β) = Σ w_n n^{-β}
        EXTRACT:  Φ = log Z_D
        CURVE:    g(β) = d²Φ/dβ² = Var_β(log n)
    """
    if weights is None:
        weights = compute_weights(n_max)

    # ASSEMBLE: partition function
    Z = 0.0
    for n in range(2, n_max + 1):
        if weights[n] > 0:
            Z += weights[n] * n ** (-beta)

    if Z <= 0:
        return DirichletEnsemble(
            beta=beta, n_max=n_max, partition_fn=0,
            free_energy=float('-inf'), mean_energy=0,
            energy_variance=0, specific_heat=0, entropy=0,
        )

    # EXTRACT: free energy
    Phi = math.log(Z)

    # CURVE: compute moments of log n under Gibbs measure
    # ⟨log n⟩ = Σ P(n) log n
    # ⟨(log n)²⟩ = Σ P(n) (log n)²
    mean_E = 0.0
    mean_E2 = 0.0
    for n in range(2, n_max + 1):
        if weights[n] > 0:
            p_n = weights[n] * n ** (-beta) / Z
            log_n = math.log(n)
            mean_E += p_n * log_n
            mean_E2 += p_n * log_n ** 2

    # Fisher-Rao metric = variance of energy
    var_E = mean_E2 - mean_E ** 2
    var_E = max(var_E, 0.0)  # numerical safety

    # Specific heat: C = β² Var(E)
    specific_heat = beta ** 2 * var_E

    # Entropy: S = βE + log Z
    entropy = beta * mean_E + Phi

    return DirichletEnsemble(
        beta=beta, n_max=n_max,
        partition_fn=Z, free_energy=Phi,
        mean_energy=mean_E, energy_variance=var_E,
        specific_heat=specific_heat, entropy=entropy,
    )


# ═══════════════════════════════════════════════════════════════
#  §3. LAPLACE ENSEMBLE (MINKOWSKI / STIFFNESS PICTURE)
# ═══════════════════════════════════════════════════════════════

@dataclass
class LaplaceEnsemble:
    """
    The Laplace (Minkowski / stiffness) statistical ensemble.

    Z_L(y) = Σ_{n≥1} w_n e^{-2πyn}

    Energy: E_n = n (linear)
    Degeneracy: w_n
    Inverse temperature: β = 2πy
    Gibbs measure: P_y(n) = w_n e^{-2πyn} / Z_L(y)

    This is the "Minkowski" picture because:
    - The mass gap lives here (μ_N(y) > 0)
    - y is the modular depth (imaginary part of τ)
    - The stiffness function is a filtered version of Z_L
    """
    y: float
    n_max: int
    partition_fn: float        # Z_L(y)
    free_energy: float         # Φ = log Z_L(y)
    mean_energy: float         # ⟨E⟩ = ⟨n⟩_y
    energy_variance: float     # Var_y(n) = Fisher-Rao metric
    specific_heat: float       # C = (2πy)² Var_y(n)
    stiffness_N3: float        # μ_3(y) with gauge filter
    entropy: float             # S = 2πy⟨n⟩ + log Z_L


def laplace_ensemble(y: float, n_max: int = 1000, N: int = 3,
                     weights: Optional[List[float]] = None) -> LaplaceEnsemble:
    """
    Compute the Laplace ensemble at modular depth y.

    The pipeline (same five primitives, different ASSEMBLE):
        ITERATE: primes → divisors → f(n)
        DIVIDE:  n/d inside f(n)
        ASSEMBLE: Z_L(y) = Σ w_n e^{-2πyn}
        EXTRACT:  Φ = log Z_L
        CURVE:    g(y) = d²Φ/dy² = (2π)² Var_y(n)
    """
    if weights is None:
        weights = compute_weights(n_max)

    two_pi = 2 * math.pi

    # ASSEMBLE: partition function
    Z = 0.0
    for n in range(2, n_max + 1):
        if weights[n] > 0:
            Z += weights[n] * math.exp(-two_pi * y * n)

    if Z <= 0:
        return LaplaceEnsemble(
            y=y, n_max=n_max, partition_fn=0,
            free_energy=float('-inf'), mean_energy=0,
            energy_variance=0, specific_heat=0,
            stiffness_N3=0, entropy=0,
        )

    # EXTRACT: free energy
    Phi = math.log(Z)

    # CURVE: moments of n under Gibbs measure
    mean_n = 0.0
    mean_n2 = 0.0
    stiffness = 0.0
    for n in range(2, n_max + 1):
        if weights[n] > 0:
            boltz = math.exp(-two_pi * y * n)
            p_n = weights[n] * boltz / Z
            mean_n += p_n * n
            mean_n2 += p_n * n * n
            # Stiffness with gauge filter
            stiffness += weights[n] * n * n * boltz * (1 - math.cos(two_pi * n / N))

    var_n = mean_n2 - mean_n ** 2
    var_n = max(var_n, 0.0)

    beta_eff = two_pi * y
    specific_heat = beta_eff ** 2 * var_n
    entropy = beta_eff * mean_n + Phi

    return LaplaceEnsemble(
        y=y, n_max=n_max,
        partition_fn=Z, free_energy=Phi,
        mean_energy=mean_n, energy_variance=var_n,
        specific_heat=specific_heat,
        stiffness_N3=stiffness, entropy=entropy,
    )


# ═══════════════════════════════════════════════════════════════
#  §4. THE WICK ROTATION MAP
# ═══════════════════════════════════════════════════════════════

@dataclass
class WickRotation:
    """
    The arithmetical Wick rotation between Laplace and Dirichlet pictures.

    The rotation is the Mellin transform:
        Laplace → Dirichlet: replace e^{-2πyn} with n^{-β}
        Dirichlet → Laplace: replace n^{-β} with e^{-2πyn}

    The key correspondence:
        Laplace parameter y  ←→  Dirichlet parameter β
        Energy E = n         ←→  Energy E = log n
        Additive structure   ←→  Multiplicative structure
        Mass gap             ←→  Zeta zeros
    """
    y: float                    # Laplace parameter
    beta: float                 # Dirichlet parameter
    laplace: LaplaceEnsemble
    dirichlet: DirichletEnsemble
    curvature_ratio: float      # g_D(β) / g_L(y)
    entropy_difference: float   # S_D - S_L
    mean_energy_ratio: float    # ⟨E⟩_D / ⟨E⟩_L


def wick_rotate(y: float, beta: float, n_max: int = 1000,
                weights: Optional[List[float]] = None) -> WickRotation:
    """
    Compute both ensembles at matched parameters and compare.

    The Wick rotation is the conceptual move: same weights w_n,
    different assembly. The comparison reveals how the arithmetic
    structure looks different in the two pictures.
    """
    if weights is None:
        weights = compute_weights(n_max)

    lap = laplace_ensemble(y, n_max, weights=weights)
    dir_ = dirichlet_ensemble(beta, n_max, weights=weights)

    curv_ratio = (dir_.energy_variance / lap.energy_variance
                  if lap.energy_variance > 1e-15 else float('inf'))
    ent_diff = dir_.entropy - lap.entropy
    energy_ratio = (dir_.mean_energy / lap.mean_energy
                    if lap.mean_energy > 1e-15 else float('inf'))

    return WickRotation(
        y=y, beta=beta,
        laplace=lap, dirichlet=dir_,
        curvature_ratio=curv_ratio,
        entropy_difference=ent_diff,
        mean_energy_ratio=energy_ratio,
    )


# ═══════════════════════════════════════════════════════════════
#  §5. CRITICAL POINT ANALYSIS
# ═══════════════════════════════════════════════════════════════

@dataclass
class CriticalPoint:
    """A critical point of a Fisher-Rao metric."""
    parameter: float           # y or β at the critical point
    curvature: float           # g(y) or g(β) at the point
    curvature_deriv: float     # dg/dy or dg/dβ (≈ 0 at critical)
    ensemble_type: str         # "Laplace" or "Dirichlet"
    critical_type: str         # "maximum", "minimum", "inflection"


def find_critical_points_dirichlet(
        beta_min: float = 1.1, beta_max: float = 10.0,
        n_points: int = 500, n_max: int = 1000,
        weights: Optional[List[float]] = None) -> List[CriticalPoint]:
    """
    Find critical points of the Dirichlet Fisher-Rao metric g_D(β).

    These are the β values where the curvature extremizes —
    "fixed points of the pipeline" in the Euclidean picture.
    """
    if weights is None:
        weights = compute_weights(n_max)

    # Sample g_D(β) at many points
    betas = [beta_min + i * (beta_max - beta_min) / n_points
             for i in range(n_points + 1)]
    curvatures = []
    for beta in betas:
        ens = dirichlet_ensemble(beta, n_max, weights=weights)
        curvatures.append(ens.energy_variance)

    # Find sign changes in the discrete derivative
    criticals = []
    for i in range(1, len(curvatures) - 1):
        d_left = curvatures[i] - curvatures[i - 1]
        d_right = curvatures[i + 1] - curvatures[i]

        if d_left * d_right < 0:
            # Sign change in derivative → critical point
            if d_left > 0 and d_right < 0:
                ctype = "maximum"
            else:
                ctype = "minimum"

            criticals.append(CriticalPoint(
                parameter=betas[i],
                curvature=curvatures[i],
                curvature_deriv=(d_left + d_right) / 2,
                ensemble_type="Dirichlet",
                critical_type=ctype,
            ))

    return criticals


def find_critical_points_laplace(
        y_min: float = 0.01, y_max: float = 1.0,
        n_points: int = 500, n_max: int = 1000,
        weights: Optional[List[float]] = None) -> List[CriticalPoint]:
    """
    Find critical points of the Laplace Fisher-Rao metric g_L(y).

    These are the y values where Var_y(n) extremizes.
    The confinement depth y_c ≈ 0.18174 should appear as
    a critical point of the stiffness-related curvature.
    """
    if weights is None:
        weights = compute_weights(n_max)

    ys = [y_min + i * (y_max - y_min) / n_points
          for i in range(n_points + 1)]
    curvatures = []
    for y in ys:
        ens = laplace_ensemble(y, n_max, weights=weights)
        curvatures.append(ens.energy_variance)

    criticals = []
    for i in range(1, len(curvatures) - 1):
        d_left = curvatures[i] - curvatures[i - 1]
        d_right = curvatures[i + 1] - curvatures[i]

        if d_left * d_right < 0:
            if d_left > 0 and d_right < 0:
                ctype = "maximum"
            else:
                ctype = "minimum"

            criticals.append(CriticalPoint(
                parameter=ys[i],
                curvature=curvatures[i],
                curvature_deriv=(d_left + d_right) / 2,
                ensemble_type="Laplace",
                critical_type=ctype,
            ))

    return criticals


# ═══════════════════════════════════════════════════════════════
#  §6. CURVATURE PROFILES — NUMERICAL COMPUTATION
# ═══════════════════════════════════════════════════════════════

@dataclass
class CurvatureProfile:
    """Full curvature profile of an ensemble over a parameter range."""
    parameters: List[float]
    partition_fn: List[float]
    free_energy: List[float]
    mean_energy: List[float]
    curvature: List[float]       # Fisher-Rao metric = Var(E)
    specific_heat: List[float]
    entropy: List[float]
    ensemble_type: str


def dirichlet_profile(beta_min: float = 1.1, beta_max: float = 10.0,
                      n_points: int = 200, n_max: int = 1000,
                      weights: Optional[List[float]] = None) -> CurvatureProfile:
    """Compute the full Dirichlet curvature profile."""
    if weights is None:
        weights = compute_weights(n_max)

    betas = [beta_min + i * (beta_max - beta_min) / n_points
             for i in range(n_points + 1)]

    Z_vals, Phi_vals, E_vals, g_vals, C_vals, S_vals = [], [], [], [], [], []

    for beta in betas:
        ens = dirichlet_ensemble(beta, n_max, weights=weights)
        Z_vals.append(ens.partition_fn)
        Phi_vals.append(ens.free_energy)
        E_vals.append(ens.mean_energy)
        g_vals.append(ens.energy_variance)
        C_vals.append(ens.specific_heat)
        S_vals.append(ens.entropy)

    return CurvatureProfile(
        parameters=betas, partition_fn=Z_vals, free_energy=Phi_vals,
        mean_energy=E_vals, curvature=g_vals, specific_heat=C_vals,
        entropy=S_vals, ensemble_type="Dirichlet",
    )


def laplace_profile(y_min: float = 0.01, y_max: float = 1.0,
                    n_points: int = 200, n_max: int = 1000,
                    weights: Optional[List[float]] = None) -> CurvatureProfile:
    """Compute the full Laplace curvature profile."""
    if weights is None:
        weights = compute_weights(n_max)

    ys = [y_min + i * (y_max - y_min) / n_points
          for i in range(n_points + 1)]

    Z_vals, Phi_vals, E_vals, g_vals, C_vals, S_vals = [], [], [], [], [], []

    for y in ys:
        ens = laplace_ensemble(y, n_max, weights=weights)
        Z_vals.append(ens.partition_fn)
        Phi_vals.append(ens.free_energy)
        E_vals.append(ens.mean_energy)
        g_vals.append(ens.energy_variance)
        C_vals.append(ens.specific_heat)
        S_vals.append(ens.entropy)

    return CurvatureProfile(
        parameters=ys, partition_fn=Z_vals, free_energy=Phi_vals,
        mean_energy=E_vals, curvature=g_vals, specific_heat=C_vals,
        entropy=S_vals, ensemble_type="Laplace",
    )


# ═══════════════════════════════════════════════════════════════
#  §7. SKELETON VS BULK COMPARISON
# ═══════════════════════════════════════════════════════════════

def skeleton_vs_bulk_dirichlet(beta: float, n_max: int = 1000) -> Dict[str, float]:
    """
    Compare Dirichlet ensemble built from full weights vs skeleton weights.

    The skeleton (primes only) is the EXTRACT of the full assembly.
    The ratio measures how much of the Dirichlet curvature comes from
    the prime skeleton vs. the composite (Lambertization) contribution.
    """
    full_w = compute_weights(n_max)
    skel_w = compute_skeleton_weights(n_max)

    full_ens = dirichlet_ensemble(beta, n_max, weights=full_w)
    skel_ens = dirichlet_ensemble(beta, n_max, weights=skel_w)

    curv_ratio = (skel_ens.energy_variance / full_ens.energy_variance
                  if full_ens.energy_variance > 1e-15 else 0.0)

    return {
        "beta": beta,
        "g_full": full_ens.energy_variance,
        "g_skeleton": skel_ens.energy_variance,
        "skeleton_fraction": curv_ratio,
        "Z_full": full_ens.partition_fn,
        "Z_skeleton": skel_ens.partition_fn,
    }


# ═══════════════════════════════════════════════════════════════
#  §8. WICK ROTATION AT MTFT CRITICAL DEPTHS
# ═══════════════════════════════════════════════════════════════

def wick_at_critical_depths(n_max: int = 1000) -> Dict[str, WickRotation]:
    """
    Compute the Wick rotation at MTFT's three critical depths.

    y_s1 = 0.1236  (skeleton zero 1)
    y_c  = 0.18174 (confinement depth / bulk zero)
    y_s2 = 0.2106  (skeleton zero 2)

    For each y, we pair it with β = 2πy (the natural correspondence)
    and compare the two ensembles.
    """
    weights = compute_weights(n_max)
    two_pi = 2 * math.pi

    results = {}
    for label, y in [("y_s1", Y_S1), ("y_c", Y_C), ("y_s2", Y_S2)]:
        beta = two_pi * y  # natural pairing
        wr = wick_rotate(y, beta, n_max, weights=weights)
        results[label] = wr

    return results


# ═══════════════════════════════════════════════════════════════
#  §9. ANALYSIS & OUTPUT
# ═══════════════════════════════════════════════════════════════

def run_full_analysis(n_max: int = 1000, verbose: bool = True):
    """Run the complete arithmetical Wick rotation analysis."""
    print("""
╔══════════════════════════════════════════════════════════════╗
║    ARITHMETICAL WICK ROTATION                               ║
║    Two Ensembles on the Same Weights                        ║
║                                                              ║
║    Laplace (mass gap) ←— Mellin —→ Dirichlet (zeta zeros)   ║
║                                                              ║
║    X₀(143) · genus 13 · y_c = 0.18174                      ║
║    Roger Tano — MTFT Research — April 2026                  ║
╚══════════════════════════════════════════════════════════════╝
""")

    weights = compute_weights(n_max)

    # §1: Wick rotation at critical depths
    print("═" * 65)
    print("  WICK ROTATION AT MTFT CRITICAL DEPTHS")
    print("═" * 65)

    results = wick_at_critical_depths(n_max)
    for label, wr in results.items():
        print(f"\n  {label}: y = {wr.y:.4f}, β = 2πy = {wr.beta:.4f}")
        print(f"    LAPLACE  (Minkowski):")
        print(f"      Z_L      = {wr.laplace.partition_fn:.6e}")
        print(f"      ⟨E⟩_L    = {wr.laplace.mean_energy:.4f}  (⟨n⟩)")
        print(f"      g_L(y)   = {wr.laplace.energy_variance:.4f}  (Var(n))")
        print(f"      C_L      = {wr.laplace.specific_heat:.4f}")
        print(f"      S_L      = {wr.laplace.entropy:.4f}")
        print(f"      μ₃(y)    = {wr.laplace.stiffness_N3:.6e}")
        print(f"    DIRICHLET (Euclidean):")
        print(f"      Z_D      = {wr.dirichlet.partition_fn:.6e}")
        print(f"      ⟨E⟩_D    = {wr.dirichlet.mean_energy:.4f}  (⟨log n⟩)")
        print(f"      g_D(β)   = {wr.dirichlet.energy_variance:.4f}  (Var(log n))")
        print(f"      C_D      = {wr.dirichlet.specific_heat:.4f}")
        print(f"      S_D      = {wr.dirichlet.entropy:.4f}")
        print(f"    ROTATION:")
        print(f"      g_D/g_L  = {wr.curvature_ratio:.6f}")
        print(f"      ΔS       = {wr.entropy_difference:.4f}")

    # §2: Critical points
    print("\n" + "═" * 65)
    print("  CRITICAL POINTS OF FISHER-RAO METRICS")
    print("═" * 65)

    print("\n  Dirichlet picture g_D(β) = Var_β(log n):")
    d_crits = find_critical_points_dirichlet(
        beta_min=1.1, beta_max=8.0, n_points=1000, n_max=n_max, weights=weights)
    if d_crits:
        for cp in d_crits:
            print(f"    β* = {cp.parameter:.4f}: g_D = {cp.curvature:.6f} "
                  f"({cp.critical_type})")
    else:
        print("    g_D(β) is monotonically decreasing — no interior extrema.")
        print("    This is expected: Var_β(log n) → 0 as β → ∞ (concentration).")

    print("\n  Laplace picture g_L(y) = Var_y(n):")
    l_crits = find_critical_points_laplace(
        y_min=0.01, y_max=0.5, n_points=1000, n_max=n_max, weights=weights)
    if l_crits:
        for cp in l_crits:
            print(f"    y* = {cp.parameter:.4f}: g_L = {cp.curvature:.6f} "
                  f"({cp.critical_type})")
    else:
        print("    g_L(y) is monotonically decreasing — no interior extrema.")

    # §3: Skeleton vs bulk
    print("\n" + "═" * 65)
    print("  SKELETON vs BULK IN DIRICHLET PICTURE")
    print("═" * 65)

    for beta in [1.5, 2.0, 3.0, 5.0]:
        sv = skeleton_vs_bulk_dirichlet(beta, n_max)
        print(f"\n  β = {beta:.1f}:")
        print(f"    g_full     = {sv['g_full']:.6f}")
        print(f"    g_skeleton = {sv['g_skeleton']:.6f}")
        print(f"    skeleton % = {sv['skeleton_fraction']*100:.1f}%")

    # §4: The bridge
    print("\n" + "═" * 65)
    print("  THE WICK ROTATION AS MELLIN BRIDGE")
    print("═" * 65)
    print("""
    LAPLACE (Minkowski)              DIRICHLET (Euclidean)
    ───────────────────              ─────────────────────
    Z_L(y) = Σ w_n e^{-2πyn}        Z_D(β) = Σ w_n n^{-β}
    Energy E = n (linear)            Energy E = log n (logarithmic)
    β_eff = 2πy                      β = Re(s)
    g_L = Var_y(n)                   g_D = Var_β(log n)
    Mass gap lives here              Zeta zeros live here
    Additive structure               Multiplicative structure

                    ←— Mellin Transform —→
                    n^{-s} = ∫ e^{-2πyn} y^{s-1} dy

    The "arithmetical Wick rotation" is:
    Replace e^{-2πyn} with n^{-β}, or equivalently,
    replace linear energy E = n with logarithmic energy E = log n.

    Same weights w_n. Same five primitives. Different physics.
    The Mellin transform IS the Wick rotation of number theory.
""")


# ═══════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        cmd = sys.argv[1]
        if cmd == "dirichlet":
            beta = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
            ens = dirichlet_ensemble(beta)
            print(f"  Dirichlet ensemble at β = {beta}:")
            print(f"    Z_D = {ens.partition_fn:.6e}")
            print(f"    g_D = {ens.energy_variance:.6f}")
        elif cmd == "laplace":
            y = float(sys.argv[2]) if len(sys.argv) > 2 else Y_C
            ens = laplace_ensemble(y)
            print(f"  Laplace ensemble at y = {y}:")
            print(f"    Z_L = {ens.partition_fn:.6e}")
            print(f"    g_L = {ens.energy_variance:.6f}")
        elif cmd == "wick":
            y = float(sys.argv[2]) if len(sys.argv) > 2 else Y_C
            beta = float(sys.argv[3]) if len(sys.argv) > 3 else 2 * math.pi * y
            wr = wick_rotate(y, beta)
            print(f"  Wick rotation: y={y}, β={beta:.4f}")
            print(f"    g_D/g_L = {wr.curvature_ratio:.6f}")
        elif cmd == "critical":
            d_crits = find_critical_points_dirichlet()
            l_crits = find_critical_points_laplace()
            print("  Dirichlet critical points:")
            for cp in d_crits:
                print(f"    β* = {cp.parameter:.4f}: {cp.critical_type}")
            print("  Laplace critical points:")
            for cp in l_crits:
                print(f"    y* = {cp.parameter:.4f}: {cp.critical_type}")
        elif cmd == "full":
            n_max = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
            run_full_analysis(n_max)
        else:
            print("Usage: python arithmetic_wick.py "
                  "[dirichlet [β]|laplace [y]|wick [y] [β]|critical|full [n_max]]")
    else:
        run_full_analysis()
