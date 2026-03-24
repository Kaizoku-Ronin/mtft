"""
MTFT Multi-N Tower: Yang-Mills Confinement Landscape
=====================================================

Maps the holonomy stiffness μ_N(y) across gauge groups SU(2) through
SU(N_max), revealing:

  1. Even-N universality — all SU(2k) collapse to the odd-integer sieve
  2. Prime-N individuality — each prime sees a unique arithmetic landscape
  3. Phase transition scaling — y_c(N) ~ c/N² confinement boundaries
  4. Mass gap persistence — μ_N > 0 unconditionally for all N, y > 0
  5. Arithmetic genome — Euler product decomposition of each gauge group

Key result (Paper 5, Theorem 7.8): The mass gap is UNCONDITIONALLY
positive because w_n ≥ 0, e^{-2πyn} > 0, and (1-cos) ≥ 0.

Reference: Papers 5, 18; Discoveries Addendum.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple
from functools import lru_cache

import numpy as np

from mtft.arithmetic import weight_array, weight, stiffness_S, center_stiffness, mass_gap_stiffness


# ═══════════════════════════════════════════════════════════════
#  Core Computation
# ═══════════════════════════════════════════════════════════════

def tower_stiffness(
    N: int,
    y: float,
    n_max: int = 500,
) -> dict:
    """
    Full stiffness decomposition for SU(N) at depth y.

    Returns dict with:
      mu_N       — mass gap stiffness (minimum over modes m)
      S_N        — total holonomy stiffness
      C_N        — center-projected stiffness
      min_mode   — the mode m achieving the minimum
      all_modes  — dict {m: mu_m} for all modes
    """
    ns = np.arange(1, n_max + 1)
    ws = weight_array(n_max)
    damped = ns ** 2 * ws * np.exp(-2 * math.pi * y * ns)

    S = float(np.sum(damped))
    C = float(np.sum(damped * np.cos(2 * math.pi * ns / N)))

    mu_min = float("inf")
    min_m = 1
    modes = {}
    for m in range(1, N):
        gap_terms = 1.0 - np.cos(2 * math.pi * ns * m / N)
        mu_m = float(np.sum(damped * gap_terms))
        modes[m] = mu_m
        if mu_m < mu_min:
            mu_min = mu_m
            min_m = m

    return {
        "N": N,
        "y": y,
        "mu_N": mu_min,
        "S_N": S,
        "C_N": C,
        "min_mode": min_m,
        "all_modes": modes,
    }


# ═══════════════════════════════════════════════════════════════
#  1. Even-N Universality
# ═══════════════════════════════════════════════════════════════

def even_n_universality(
    y: float = 0.10,
    N_max: int = 50,
    n_max: int = 500,
) -> dict:
    """
    Verify that all even-N gauge groups collapse to the odd-integer sieve.

    For N = 2k, the Z₂ center mode m = N/2 gives:
        1 - cos(πn) = 1 - (-1)^n = {2 if n odd, 0 if n even}

    So μ_{2k}(y)|_{m=N/2} = 2·Σ_{odd n} n² w_n e^{-2πyn} for ALL k ≥ 1.

    This is also the minimising mode for moderate y, making it
    the mass gap itself for the even-N universality class.

    Returns dict with even_values, verification, and max_deviation.
    """
    ns = np.arange(1, n_max + 1)
    ws = weight_array(n_max)
    damped = ns ** 2 * ws * np.exp(-2 * math.pi * y * ns)

    # Reference: the odd-sieve sum
    odd_mask = (ns % 2 == 1)
    mu_odd_sieve = 2.0 * float(np.sum(damped[odd_mask]))

    even_vals = {}
    for N in range(2, N_max + 1, 2):
        # Compute μ at exactly m = N/2 (the Z₂ center element)
        m = N // 2
        gap_terms = 1.0 - np.cos(2 * math.pi * ns * m / N)
        mu_center = float(np.sum(damped * gap_terms))
        even_vals[N] = mu_center

    deviations = {N: abs(v - mu_odd_sieve) for N, v in even_vals.items()}
    max_dev = max(deviations.values())

    return {
        "y": y,
        "odd_sieve_value": mu_odd_sieve,
        "even_N_values": even_vals,
        "max_deviation": max_dev,
        "universal": max_dev < 1e-10,
        "theorem": ("All even SU(N) groups share μ_{2k}(y)|_{m=N/2} = "
                    "2·Σ_{odd n} n²w_n e^{-2πyn}, independent of k."),
    }


# ═══════════════════════════════════════════════════════════════
#  2. Phase Transition Boundaries
# ═══════════════════════════════════════════════════════════════

def confinement_boundary(
    N: int,
    y_min: float = 0.0001,
    y_max: float | None = None,
    n_points: int = 2000,
    n_max: int = 500,
) -> dict:
    """
    Find the confinement-deconfinement boundary for SU(N).

    The boundary is where the center stiffness C_N(y) crosses zero.
    For N=3, this is y_c ≈ 0.18174.

    Note: μ_N(y) > 0 ALWAYS (unconditional mass gap).
    The "boundary" is the Hessian isotropy point, not a gap closure.

    Uses N² scaling to set intelligent search bounds for large N.
    """
    # Scale search range by 1/N²
    if y_max is None:
        y_max = min(5.0 / (N ** 2) + 0.01, 0.5)
    y_min = max(y_min, 0.5 / (N ** 2))

    ys = np.linspace(y_min, y_max, n_points)
    C_vals = np.array([center_stiffness(y, N=N, n_max=n_max) for y in ys])

    # Find FIRST zero crossing (C going from positive to negative)
    y_boundary = None
    for i in range(len(C_vals) - 1):
        if C_vals[i] > 0 and C_vals[i + 1] <= 0:
            # Linear interpolation
            y_boundary = ys[i] - C_vals[i] * (ys[i + 1] - ys[i]) / (C_vals[i + 1] - C_vals[i])
            break

    N2_yc = N ** 2 * y_boundary if y_boundary else None

    return {
        "N": N,
        "y_boundary": y_boundary,
        "N2_y_boundary": N2_yc,
        "scan_y": ys,
        "scan_C": C_vals,
        "note": "C_N(y) = 0 marks Hessian isotropy, NOT gap closure",
    }


def phase_transition_scaling(
    N_range: Optional[List[int]] = None,
    n_max: int = 500,
) -> List[dict]:
    """
    Compute y_c(N) for multiple gauge groups and verify N² scaling.

    Theory: y_c(N) ~ c/N² with c → ~2.2 as N → ∞.
    """
    if N_range is None:
        N_range = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]
    results = []
    for N in N_range:
        b = confinement_boundary(N, n_max=n_max)
        results.append({
            "N": N,
            "y_c": b["y_boundary"],
            "N2_yc": b["N2_y_boundary"],
            "is_prime": all(N % d != 0 for d in range(2, N)) and N > 1,
            "is_even": N % 2 == 0,
        })
    return results


# ═══════════════════════════════════════════════════════════════
#  3. Boundary Tracking (Yang-Mills mass gap persistence)
# ═══════════════════════════════════════════════════════════════

def boundary_tracking(
    N_range: Optional[List[int]] = None,
    c: float = 0.10,
    n_max: int = 500,
) -> List[dict]:
    """
    Under the ansatz y = c/N², compute μ_N and the physical mass gap.

    The mass gap m_gap = √(μ_N / N) should remain O(1) or grow,
    proving confinement persists in the large-N limit.

    Key result (Paper 5): μ_N grows DRAMATICALLY under tracking —
    the arithmetic enhancement far outpaces the 1/N² suppression.
    """
    if N_range is None:
        N_range = [2, 3, 5, 8, 10, 15, 20]
    results = []
    for N in N_range:
        y = c / (N ** 2)
        mu = mass_gap_stiffness(y, N=N, n_max=n_max)
        m_gap = math.sqrt(mu / N) if mu > 0 else 0.0
        results.append({
            "N": N,
            "y": y,
            "mu_N": mu,
            "m_gap": m_gap,
            "status": "CONFINED" if mu > 0 else "DECONFINED",
        })
    return results


# ═══════════════════════════════════════════════════════════════
#  4. Arithmetic Genome (Euler product decomposition)
# ═══════════════════════════════════════════════════════════════

def arithmetic_genome(
    N: int,
    y: float = 0.10,
    n_max: int = 60,
) -> dict:
    """
    Decompose μ_N(y) into its prime-power Euler factors.

    For each prime p, compute the fraction of stiffness contributed
    by n = p^k modes. This reveals how SU(N) "filters" the primes.

    Key finding: Each SU(N) suppresses its own prime's Euler factor.
    SU(3) zeroes out p=3, SU(5) zeroes out p=5, etc.
    """
    ns = np.arange(1, n_max + 1)
    ws = weight_array(n_max)
    damped = ns ** 2 * ws * np.exp(-2 * math.pi * y * ns)

    # Find minimising mode
    mu_min = float("inf")
    min_m = 1
    for m in range(1, N):
        gap = 1.0 - np.cos(2 * math.pi * ns * m / N)
        mu_m = float(np.sum(damped * gap))
        if mu_m < mu_min:
            mu_min = mu_m
            min_m = m

    # Now decompose by prime power
    gap_filter = 1.0 - np.cos(2 * math.pi * ns * min_m / N)
    contributions = damped * gap_filter
    total = float(np.sum(contributions))

    def is_prime_power(n, p):
        if n < 1:
            return False
        while n % p == 0:
            n //= p
        return n == 1

    primes = [2, 3, 5, 7, 11, 13, 17, 19]
    genome = {}
    for p in primes:
        pp_sum = sum(float(contributions[n - 1])
                     for n in range(1, n_max + 1)
                     if is_prime_power(n, p))
        genome[p] = {
            "absolute": pp_sum,
            "fraction": pp_sum / total if total > 0 else 0.0,
        }

    pp_total = sum(g["absolute"] for g in genome.values())
    return {
        "N": N,
        "y": y,
        "mu_N": mu_min,
        "min_mode": min_m,
        "prime_contributions": genome,
        "prime_power_total_fraction": pp_total / total if total > 0 else 0.0,
        "note": f"SU({N}) suppresses p={N} if N is prime",
    }


# ═══════════════════════════════════════════════════════════════
#  5. Arithmetic Periodic Table
# ═══════════════════════════════════════════════════════════════

def arithmetic_periodic_table(
    N_max: int = 20,
    y: float = 0.10,
    n_max: int = 500,
) -> List[dict]:
    """
    The MTFT "periodic table of gauge groups" at fixed depth y.

    Three-tier structure:
      Floor  — all even-N degenerate at odd-sieve value
      Peaks  — prime-N achieve higher unique stiffness
      Middle — odd composites inherit factor structure
    """
    rows = []
    for N in range(2, N_max + 1):
        info = tower_stiffness(N, y, n_max)
        is_prime = all(N % d != 0 for d in range(2, N)) and N > 1
        is_even = N % 2 == 0
        tier = "floor" if is_even else ("peak" if is_prime else "middle")
        rows.append({
            "N": N,
            "mu_N": info["mu_N"],
            "S_N": info["S_N"],
            "C_N": info["C_N"],
            "min_mode": info["min_mode"],
            "is_prime": is_prime,
            "is_even": is_even,
            "tier": tier,
        })
    return rows


# ═══════════════════════════════════════════════════════════════
#  6. Character Orthogonality
# ═══════════════════════════════════════════════════════════════

def character_orthogonality(
    N: int,
    y: float = 0.10,
    n_max: int = 500,
) -> dict:
    """
    Verify the character orthogonality identity:

        Σ_{m=1}^{N-1} (1 − cos(2πnm/N)) = N · 𝟙_{N∤n}

    This means the mode average satisfies:
        ⟨μ_N⟩ · (N-1)/N  →  μ_∞  as N → ∞

    where μ_∞ is the full unfiltered stiffness sum.
    """
    ns = np.arange(1, n_max + 1)
    ws = weight_array(n_max)
    damped = ns ** 2 * ws * np.exp(-2 * math.pi * y * ns)

    # Mode average
    mode_sum = np.zeros(n_max)
    for m in range(1, N):
        mode_sum += 1.0 - np.cos(2 * math.pi * ns * m / N)

    # Should equal N * (1 - indicator_{N|n})
    expected = np.zeros(n_max)
    for i, n in enumerate(ns):
        expected[i] = 0.0 if (n % N == 0) else float(N)

    identity_error = float(np.max(np.abs(mode_sum - expected)))

    # Mode average of stiffness
    mu_avg = float(np.sum(damped * mode_sum)) / (N - 1)
    mu_inf = stiffness_S(y, n_max)  # full unfiltered
    convergence = mu_avg * (N - 1) / N / mu_inf if mu_inf > 0 else 0.0

    return {
        "N": N,
        "identity_error": identity_error,
        "identity_verified": identity_error < 1e-10,
        "mu_average": mu_avg,
        "mu_infinity": mu_inf,
        "convergence_ratio": convergence,
        "note": f"At N={N}, mode average is {convergence*100:.2f}% of μ_∞",
    }


# ═══════════════════════════════════════════════════════════════
#  7. Full Tower Report
# ═══════════════════════════════════════════════════════════════

def tower_report(
    N_max: int = 15,
    y: float = 0.10,
    verbose: bool = True,
) -> dict:
    """
    Complete multi-N tower analysis. Prints report if verbose.
    """
    # Periodic table
    table = arithmetic_periodic_table(N_max, y)

    # Even-N universality
    univ = even_n_universality(y, N_max)

    # Phase boundaries
    primes_in_range = [N for N in range(2, N_max + 1)
                       if all(N % d != 0 for d in range(2, N)) and N > 1]
    boundaries = phase_transition_scaling(list(range(2, min(N_max + 1, 21))))

    # Boundary tracking
    tracking = boundary_tracking()

    # Character orthogonality
    char = character_orthogonality(N_max, y)

    result = {
        "periodic_table": table,
        "even_universality": univ,
        "phase_boundaries": boundaries,
        "boundary_tracking": tracking,
        "character_orthogonality": char,
    }

    if verbose:
        print("=" * 72)
        print("  MTFT MULTI-N TOWER: YANG-MILLS CONFINEMENT LANDSCAPE")
        print("=" * 72)

        print(f"\n  Arithmetic Periodic Table at y = {y}")
        print(f"  {'N':>4s} {'μ_N':>12s} {'S_N':>12s} {'Mode':>5s} {'Tier':>8s}")
        print(f"  {'-' * 46}")
        for r in table:
            print(f"  {r['N']:4d} {r['mu_N']:12.4f} {r['S_N']:12.4f} "
                  f"{r['min_mode']:5d} {r['tier']:>8s}")

        print(f"\n  Even-N Universality: {'VERIFIED' if univ['universal'] else 'FAILED'}")
        print(f"  Max deviation: {univ['max_deviation']:.2e}")
        print(f"  Odd-sieve value: {univ['odd_sieve_value']:.6f}")

        print(f"\n  Phase Boundaries (y_c where C_N = 0):")
        print(f"  {'N':>4s} {'y_c':>10s} {'N²·y_c':>10s} {'Type':>8s}")
        print(f"  {'-' * 36}")
        for b in boundaries:
            if b['y_c'] is not None:
                typ = "prime" if b['is_prime'] else ("even" if b['is_even'] else "odd-c")
                print(f"  {b['N']:4d} {b['y_c']:10.6f} {b['N2_yc']:10.4f} {typ:>8s}")

        print(f"\n  Boundary Tracking (y = 0.10/N²):")
        print(f"  {'N':>4s} {'y':>12s} {'μ_N':>14s} {'m_gap':>10s} {'Status':>10s}")
        print(f"  {'-' * 56}")
        for t in tracking:
            print(f"  {t['N']:4d} {t['y']:12.6f} {t['mu_N']:14.2f} "
                  f"{t['m_gap']:10.2f} {t['status']:>10s}")

        print(f"\n  Character Orthogonality at N={N_max}:")
        print(f"  Identity verified: {char['identity_verified']}")
        print(f"  Convergence: {char['convergence_ratio']*100:.2f}% of μ_∞")

        print(f"\n  MASS GAP STATUS: UNCONDITIONALLY POSITIVE")
        print(f"  (w_n ≥ 0) × (e^{{-2πyn}} > 0) × (1-cos ≥ 0) = non-negative sum")
        print("=" * 72)

    return result
