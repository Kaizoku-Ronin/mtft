"""
MTFT × Riemann: Explicit Formula and ζ-Zero Connection
========================================================

The MTFT stiffness admits an explicit formula decomposition:

    μ_N(y) = Main(y) + Σ_ρ Φ_N(y; ρ) + Trivial(y)

where ρ = 1/2 + iγ_k runs over the nontrivial zeros of ζ(s).

SUPERSESSION NOTICE (July 2026, audit B4/B5)
--------------------------------------------
The old headline "RH ⟺ κ(y) ≥ 0" is FALSE (κ^Λ(y) < 0 unconditionally,
verified June 2026).  The corrected equivalence — Theorem 1 of
"The Corrected RH Equivalence" (July 2026) — is implemented in the
final section of this module:

    RH  ⟺  limsup_{y→0⁺} |𝒟(y)| < ∞,
    𝒟(y) := Δκ(y)·X^{−3/2},   Δκ := κ^Λ − κ_Main (stable normal form Df 3),

evaluated on the skeleton stiffness μ^Λ(y) = Σ_n w_n e^{−2πny} with
w_n = Σ_{dm=n} d²·mΛ(m).  On true zeros the envelope log-slope of |𝒟|
is ≈ 0 (bounded); a synthetic off-line quadruplet at β₀ gives slope
≈ 1/2 − β₀ (divergence).  See corrected_rh_diagnostic().

The legacy functions below are retained for compatibility; read
rh_diagnostic()'s docstring before using it.

Reference: Papers 5, 18, 23; Bridge Paper 14; July 2026 RH draft.
"""

from __future__ import annotations

import math
import cmath
from typing import List, Optional, Tuple

import numpy as np

from mtft.arithmetic import mass_gap_stiffness, stiffness_S, weight_array


# ═══════════════════════════════════════════════════════════════
#  Riemann Zeros (first 30, from LMFDB to 30+ digits)
# ═══════════════════════════════════════════════════════════════

RIEMANN_ZEROS: List[float] = [
    14.134725141734693790,
    21.022039638771554993,
    25.010857580145688763,
    30.424876125859513210,
    32.935061587739189691,
    37.586178158825671257,
    40.918719012147495187,
    43.327073280914999519,
    48.005150881167159727,
    49.773832477672302181,
    52.970321477714460644,
    56.446247697063394804,
    59.347044002602353079,
    60.831778524609809844,
    65.112544048081606660,
    67.079810529494173714,
    69.546401711173979253,
    72.067157674481907582,
    75.704690699083933168,
    77.144840068874805372,
    79.337375020249367922,
    82.910380854086030183,
    84.735492980517050105,
    87.425274613125229406,
    88.809111207634465423,
    92.491899270558484070,
    94.651344040519837726,
    95.870634228245309758,
    98.831194218193692234,
    101.317851005731391228,
]


# ═══════════════════════════════════════════════════════════════
#  Gamma suppression
# ═══════════════════════════════════════════════════════════════

def gamma_suppression(gamma_k: float) -> float:
    """
    The Gamma-factor suppression for the k-th zero:

        |Γ(3/2 + iγ)| ~ e^{-πγ/2} · γ^{1} · (const)

    This kills the tail extraordinarily fast.
    """
    try:
        # |Γ(3/2 + iγ)| via the reflection formula
        s = complex(1.5, gamma_k)
        # Use Stirling: |Γ(σ+it)| ~ √(2π) |t|^{σ-1/2} e^{-π|t|/2}
        return math.sqrt(2 * math.pi) * abs(gamma_k) ** 1.0 * math.exp(-math.pi * abs(gamma_k) / 2)
    except OverflowError:
        return 0.0


def gamma_suppression_table(n_zeros: int = 30) -> List[dict]:
    """
    Show how fast the Gamma factor kills the zero-sum tail.

    Key result: e^{-πγ/2} suppresses γ₁=14.13 by ~10⁻¹⁸.
    First 5-10 zeros capture everything numerically.
    """
    rows = []
    for i, gamma in enumerate(RIEMANN_ZEROS[:n_zeros]):
        sup = gamma_suppression(gamma)
        rows.append({
            "k": i + 1,
            "gamma_k": gamma,
            "suppression": sup,
            "log10_suppression": math.log10(sup) if sup > 0 else float("-inf"),
        })
    return rows


# ═══════════════════════════════════════════════════════════════
#  Zero Contribution to Stiffness
# ═══════════════════════════════════════════════════════════════

def zero_contribution(
    y: float,
    gamma_k: float,
    N: int = 3,
    beta: float = 0.5,
) -> complex:
    """
    Single zero contribution to the stiffness explicit formula:

        Φ_N(y; ρ) ~ (2πy)^{-(ρ+1)} · Γ(ρ+1) · B_N(ρ)

    where ρ = β + iγ_k and B_N is the SU(N) geometry factor.

    For the SU(3) center projection (3∤n filter), B_3 involves
    the local cosine structure: C_3(ρ+1).

    Returns complex contribution (take real part for μ_N).
    """
    rho = complex(beta, gamma_k)

    # (2πy)^{-(ρ+1)}
    base = 2.0 * math.pi * y
    if base <= 0:
        return 0j
    power = -(rho + 1)
    term1 = base ** power.real * cmath.exp(1j * power.imag * math.log(base))

    # Gamma suppression (Stirling approximation for large γ)
    sup = gamma_suppression(gamma_k)
    phase = gamma_k * math.log(abs(gamma_k) + 1e-30) - gamma_k  # Stirling phase
    term2 = sup * cmath.exp(1j * phase)

    # SU(N) geometry factor: simplified as ~N/(N-1) for the filter
    B_N = N / (N - 1.0) if N > 1 else 1.0

    return term1 * term2 * B_N


def zero_sum(
    y: float,
    N: int = 3,
    n_zeros: int = 30,
    beta: float = 0.5,
) -> float:
    """
    Sum of all zero contributions (oscillatory part):

        Σ_ρ Re[Φ_N(y; ρ)]

    Pairs ρ and ρ̄ contribute conjugate terms, so the sum is real.
    """
    total = 0.0
    for gamma in RIEMANN_ZEROS[:n_zeros]:
        # ρ and ρ̄ pair: Φ(ρ) + Φ(ρ̄) = 2·Re[Φ(ρ)]
        phi = zero_contribution(y, gamma, N, beta)
        total += 2.0 * phi.real
    return total


# ═══════════════════════════════════════════════════════════════
#  Explicit Formula Reconstruction
# ═══════════════════════════════════════════════════════════════

def explicit_formula(
    y: float,
    N: int = 3,
    n_zeros: int = 30,
    beta: float = 0.5,
    n_max: int = 500,
) -> dict:
    """
    Full explicit formula decomposition:

        μ_N(y) = Main(y) + Σ_ρ Φ_N(y; ρ) + Trivial(y)

    Returns dict with:
      direct    — μ_N from direct summation (exact)
      main      — smooth power-law envelope (fitted)
      zero_sum  — oscillatory contribution from ζ zeros
      trivial   — correction from trivial zeros (small)
      total     — main + zero_sum + trivial
      residual  — |direct - total| / |direct|
    """
    direct = mass_gap_stiffness(y, N=N, n_max=n_max)

    # Main term: smooth envelope ~ A·y^{-3} for small y
    # Fit A from the known value at a reference point
    y_ref = 0.5
    mu_ref = mass_gap_stiffness(y_ref, N=N, n_max=n_max)
    A = mu_ref * y_ref ** 3
    main = A * y ** (-3)

    # Zero sum
    zs = zero_sum(y, N, n_zeros, beta)

    # Trivial zeros: contribute at s = -2, -4, ...
    # These are exponentially small for y > 0.01
    trivial = 0.0
    for k in range(1, 6):
        s_triv = -2 * k
        trivial += (-1) ** k * (2 * math.pi * y) ** (2 * k - 1) / math.factorial(2 * k)
    trivial *= 0.01  # scale factor (subdominant)

    total = main + zs + trivial

    return {
        "y": y,
        "N": N,
        "direct": direct,
        "main": main,
        "zero_sum": zs,
        "trivial": trivial,
        "total": total,
        "residual": abs(direct - total) / abs(direct) if direct != 0 else 0.0,
        "n_zeros_used": min(n_zeros, len(RIEMANN_ZEROS)),
    }


# ═══════════════════════════════════════════════════════════════
#  Bakry-Émery Curvature (RH Diagnostic)
# ═══════════════════════════════════════════════════════════════

def bakry_emery_curvature(
    y: float,
    N: int = 3,
    n_max: int = 500,
    dy_frac: float = 0.01,
) -> float:
    """
    Bakry-Émery curvature of the arithmetic vacuum:

        κ(y) = d²/d(ln y)² ln μ_N(y)

    LEGACY FRAMING — the criterion "RH ⟺ κ(y) ≥ 0 for all y > 0" once
    attached to this quantity is FALSE (κ^Λ(y) < 0 unconditionally;
    audit B4, June 2026 verification).  The corrected RH equivalence is
    the boundedness diagnostic of corrected_rh_diagnostic() (Th 1, July
    2026).  This function remains valid as a curvature evaluator.

    Computed numerically via central differences.
    """
    mu = mass_gap_stiffness(y, N=N, n_max=n_max)
    if mu <= 0:
        return float("-inf")

    dy = y * dy_frac
    mu_plus = mass_gap_stiffness(y + dy, N=N, n_max=n_max)
    mu_minus = mass_gap_stiffness(y - dy, N=N, n_max=n_max)

    if mu_plus <= 0 or mu_minus <= 0:
        return float("-inf")

    ln_mu = math.log(mu)
    ln_mu_p = math.log(mu_plus)
    ln_mu_m = math.log(mu_minus)

    ln_y = math.log(y)
    dln_y = math.log((y + dy) / y)

    # d²(ln μ)/d(ln y)²
    kappa = (ln_mu_p - 2 * ln_mu + ln_mu_m) / (dln_y ** 2)
    return kappa


def rh_diagnostic(
    y_min: float = 0.02,
    y_max: float = 1.0,
    n_points: int = 100,
    N: int = 3,
    n_max: int = 500,
) -> dict:
    """
    SUPERSEDED DIAGNOSTIC (audit B4).  The criterion "RH ⟺ κ(y) ≥ 0
    everywhere" is false: κ^Λ(y) < 0 unconditionally for all y > 0
    (verified June 2026).  The RH-sensitive object is the BOUNDEDNESS of
    the normalized oscillation 𝒟(y), not the sign of the curvature —
    use corrected_rh_diagnostic() (July 2026 Th 1) instead.

    Legacy behavior preserved: scans Bakry-Émery curvature across y.

    Returns:
      curvature_values: array of κ(y)
      curvature_positive: True if κ ≥ 0 everywhere
      min_curvature: minimum κ found
      y_at_min: depth where minimum occurs
      superseded: True (always)
    """
    ys = np.linspace(y_min, y_max, n_points)
    kappas = np.array([bakry_emery_curvature(y, N, n_max) for y in ys])

    min_idx = np.argmin(kappas)

    return {
        "y_values": ys,
        "curvature_values": kappas,
        "curvature_positive": bool(np.all(kappas >= -1e-10)),
        "min_curvature": float(kappas[min_idx]),
        "y_at_min": float(ys[min_idx]),
        "N": N,
        "consistent_with_RH": bool(np.all(kappas >= -1e-10)),
        "superseded": True,
        "supersession_note": "kappa >= 0 is NOT an RH criterion (false since "
                             "June 2026 verification). Use corrected_rh_diagnostic().",
    }


# ═══════════════════════════════════════════════════════════════
#  Tower Rigidity
# ═══════════════════════════════════════════════════════════════

def tower_rigidity(
    y: float = 0.1,
    beta_values: Optional[List[float]] = None,
    n_zeros: int = 10,
) -> dict:
    """
    Test how off-critical-line zeros would affect the stiffness.

    On the critical line (β=0.5): |Φ| ~ y^{-3/2} — balanced.
    Off the line (β>0.5): |Φ| ~ y^{-(β+1)} — dominates at small y.

    A single off-line zero creates exponentially growing oscillations
    (the tower rigidity principle).  Note: the "curvature positivity"
    phrasing below is legacy framing — κ ≥ 0 is not an RH criterion
    (audit B4); the live statement is the envelope-slope divergence of
    corrected_rh_diagnostic().
    """
    if beta_values is None:
        beta_values = [0.5, 0.6, 0.75, 0.9]

    results = {}
    for beta in beta_values:
        envelope = sum(
            abs(zero_contribution(y, gamma, N=3, beta=beta))
            for gamma in RIEMANN_ZEROS[:n_zeros]
        )
        scaling = y ** (-(beta + 1))
        results[beta] = {
            "envelope": envelope,
            "scaling_exponent": -(beta + 1),
            "y_scaling": scaling,
            "dominates": beta > 0.5,
        }

    return {
        "y": y,
        "beta_envelopes": results,
        "conclusion": (
            "On-line (β=0.5): balanced oscillations. "
            "Off-line (β>0.5): envelope grows as y^{-(β+1)} — an "
            "unbounded normalized oscillation (the live Th 1 form; the "
            "older 'curvature positivity' phrasing is legacy, audit B4). "
            "A single off-line zero would break the entire N-tower."
        ),
    }


# ═══════════════════════════════════════════════════════════════
#  Corrected RH Equivalence (July 2026 Th 1) — requires mpmath
# ═══════════════════════════════════════════════════════════════
#
# Skeleton stiffness:  μ^Λ(y) = Σ_n w_n e^{−nX}, X = 2πy,
#   w_n = Σ_{dm=n} d²·mΛ(m)   (von Mangoldt Λ)
# Mellin parent: F(s) = ζ(s−2)·G(s−1), G(w) = −ζ′(w)/ζ(w).
#
# Explicit formula (Pr 1):  μ^Λ = M + Z + T,
#   M(X) = 2G(2)X^{−3} + ζ(0)X^{−2}   (zero-free part; trivial tower T
#                                     is O(X³ ln X), negligible for y ≤ 10⁻²)
#   Z(X) = Σ_ρ m_ρ c_ρ X^{−(ρ+1)},    c_ρ = −Γ(ρ+1)ζ(ρ−1)
#
# Diagnostic (Df 2):  Δκ = D²log μ^Λ − D²log M,  𝒟 = Δκ·X^{−3/2},
# evaluated in the stable normal form (Df 3):
#   Δκ = D[P/Q],  P = (DZ)M − Z(DM),  DP = (D²Z)M − Z(D²M),  Q = M(M+Z).
#
# Th 1:  RH ⟺ limsup_{y→0⁺} |𝒟(y)| < ∞.
# Numerical signature: envelope slope of log|𝒟| vs log y is ≈ 0 on the
# true zeros (bounded), ≈ 1/2 − β₀ < 0 for a synthetic off-line
# quadruplet at real part β₀ (divergence as y → 0⁺).

def _mp():
    try:
        import mpmath
        return mpmath
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "The corrected RH diagnostic needs mpmath: pip install mpmath"
        ) from e


def skeleton_weights(n_max: int) -> list:
    """
    w_n = Σ_{dm=n} d²·m·Λ(m) for n = 1..n_max (exact integer arithmetic).
    """
    import math as _math
    w = [0] * (n_max + 1)
    # Λ(m) via prime powers
    lam = [0.0] * (n_max + 1)
    for p in range(2, n_max + 1):
        if lam[p] == 0.0:  # prime candidate (sieve)
            # mark p as prime: check no smaller prime divided — use simple test
            is_prime = all(p % q for q in range(2, int(p ** 0.5) + 1))
            if not is_prime:
                continue
            lp = _math.log(p)
            pk = p
            while pk <= n_max:
                lam[pk] = lp
                pk *= p
    for d in range(1, n_max + 1):
        d2 = d * d
        for m in range(2, n_max // d + 1):
            if lam[m] != 0.0:
                w[d * m] += d2 * m * lam[m]
    return w


def skeleton_stiffness(y: float, n_max: int = 2000) -> float:
    """μ^Λ(y) = Σ w_n e^{−2πny} by direct summation (float)."""
    mp = _mp()
    w = skeleton_weights(n_max)
    X = 2.0 * math.pi * y
    return float(mp.fsum(w[n] * mp.exp(-n * X) for n in range(1, n_max + 1)))


def _G(w):
    mp = _mp()
    return -mp.zeta(w, derivative=1) / mp.zeta(w)


def zero_coefficient(rho) -> complex:
    """c_ρ = −Γ(ρ+1)·ζ(ρ−1)  (Df 4 coefficient formula; ≠ 0 by Pr 2)."""
    mp = _mp()
    return -mp.gamma(rho + 1) * mp.zeta(rho - 1)


def _main_part(X):
    """M, DM, D²M at X (D = d/d ln X).  Trivial tower omitted (O(X³ ln X))."""
    mp = _mp()
    terms = [(2 * _G(2), mp.mpf(3)), (mp.zeta(0), mp.mpf(2))]
    M = DM = D2M = mp.mpf(0)
    for c, a in terms:
        t = c * X ** (-a)
        M += t
        DM += -a * t
        D2M += a * a * t
    return M, DM, D2M


def _zero_part(X, zeros):
    """Z, DZ, D²Z at X.  zeros: list of (rho, c_rho); conjugate pairs → 2 Re."""
    mp = _mp()
    Z = DZ = D2Z = mp.mpf(0)
    for rho, c in zeros:
        a = rho + 1
        t = c * X ** (-a)
        Z += 2 * mp.re(t)
        DZ += 2 * mp.re(-a * t)
        D2Z += 2 * mp.re(a * a * t)
    return Z, DZ, D2Z


def delta_kappa_stable(y: float, zeros) -> float:
    """
    Δκ(y) via the stable normal form (Df 3): exact identity, immune to
    the catastrophic cancellation of the naive two-κ subtraction.
    """
    mp = _mp()
    X = 2 * mp.pi * y
    M, DM, D2M = _main_part(X)
    Z, DZ, D2Z = _zero_part(X, zeros)
    P = DZ * M - Z * DM
    DP = D2Z * M - Z * D2M
    Q = M * (M + Z)
    DQ = DM * (M + Z) + M * (DM + DZ)
    return (DP * Q - P * DQ) / (Q * Q)


def normalized_oscillation(y: float, zeros) -> float:
    """𝒟(y) = Δκ(y)·X^{−3/2} — the Th 1 diagnostic object."""
    mp = _mp()
    X = 2 * mp.pi * y
    return delta_kappa_stable(y, zeros) * X ** (-mp.mpf(3) / 2)


def envelope_slope(zeros, y_min_pow=-7.0, y_max_pow=-1.8, n_points=None,
                   min_bin=10, samples_per_period=10) -> dict:
    """
    Binned-RMS envelope slope of log₁₀|𝒟| vs log₁₀ y over
    y ∈ [10^y_min_pow, 10^y_max_pow].  Predicted: 0 (true zeros),
    1/2 − β₀ (synthetic off-line quadruplet at Re β₀).

    A.7 discipline (v0.9.0, audit S.4/T-E4): the grid density defaults to
    `samples_per_period` samples per γ₁ oscillation period in ln X
    (estimator_standards.recommended_samples_per_decade: γ₁ ≈ 5.18
    periods/decade → 52 samples/decade at the default 10/period), and
    the fit drops every half-decade bin with fewer than
    `min_bin` samples (terminal-bin leverage guard, same rule as
    estimator_standards.binned_log_slope).  At the legacy 91-point grid
    (≈8.8 samples per half-decade bin) a bare min_bin guard discards all
    eleven bins — guard and density therefore ship together.

    Returns: slope, n_bins (used), n_bins_dropped, points (usable bins
    only), n_points, samples_per_decade.
    """
    from mtft.estimator_standards import (binned_log_slope,
                                          recommended_samples_per_decade)
    mp = _mp()
    if n_points is None:
        gamma1 = abs(float(mp.im(zeros[0][0]))) if zeros else 14.134725
        spd = recommended_samples_per_decade(gamma1,
                                             per_period=samples_per_period)
        n_points = int(math.ceil((y_max_pow - y_min_pow) * spd)) + 1
    ys, ds = [], []
    for i in range(n_points):
        ly = y_min_pow + (y_max_pow - y_min_pow) * i / (n_points - 1)
        y = mp.mpf(10) ** ly
        ys.append(y)
        ds.append(abs(normalized_oscillation(y, zeros)))
    slope, n_used, dropped = binned_log_slope(ys, ds, bin_width=0.5,
                                              min_bin=min_bin)
    bins = {}
    for y, d in zip(ys, ds):
        b = math.floor(float(mp.log10(y)) * 2) / 2
        bins.setdefault(b, []).append(float(d * d))
    pts = sorted((b, math.log10(math.sqrt(sum(v) / len(v))))
                 for b, v in bins.items() if len(v) >= min_bin)
    return {"slope": slope, "n_bins": n_used, "n_bins_dropped": dropped,
            "points": pts, "n_points": n_points,
            "samples_per_decade": n_points / (y_max_pow - y_min_pow)}


def on_line_zeros(n_zeros: int = 12) -> list:
    """First n_zeros true zeros as (ρ, c_ρ) pairs, ρ = 1/2 + iγ."""
    mp = _mp()
    half = mp.mpf(1) / 2
    return [(mp.mpc(half, g), zero_coefficient(mp.mpc(half, g)))
            for g in RIEMANN_ZEROS[:n_zeros]]


def offline_quadruplet(beta: float, gamma: float) -> list:
    """
    Synthetic off-line quadruplet (Df 4 admissible multiset):
    (β ± iγ, 1−β ± iγ) with the canonical coefficients.
    """
    mp = _mp()
    r1 = mp.mpc(beta, gamma)
    r2 = mp.mpc(1 - beta, gamma)
    return [(r1, zero_coefficient(r1)), (r2, zero_coefficient(r2))]


def corrected_rh_diagnostic(n_zeros: int = 12,
                            beta_tests=(0.6, 0.75, 0.9),
                            y_range=(-7.0, -1.8)) -> dict:
    """
    The July 2026 diagnostic, replacing "κ ≥ 0".

    Runs the envelope-slope test on:
      (a) the true zeros           — expect slope ≈ 0     (bounded 𝒟: Th 1a)
      (b) synthetic off-line quadruplets at β₀ — expect ≈ 1/2 − β₀ (Th 1b)

    Anchors (draft Appendix A.2/A.3): on-line log-slope of |Δκ| = 1.479
    (pred. 3/2, i.e. |𝒟| slope ≈ 0); off-line slopes −0.089/−0.258/−0.413
    for β₀ = 0.6/0.75/0.9 (pred. −0.1/−0.25/−0.4).
    """
    mp = _mp()
    mp.mp.dps = 30
    out = {"theorem": "RH  <==>  limsup_{y->0+} |D(y)| < infinity",
           "diagnostic": "envelope slope of log|D|: 0 bounded, 1/2-beta0 divergent"}
    on = on_line_zeros(n_zeros)
    res_on = envelope_slope(on, *y_range)
    out["on_line"] = {"slope": res_on["slope"], "predicted": 0.0,
                      "bounded": abs(res_on["slope"]) < 0.05}
    out["off_line"] = {}
    for b in beta_tests:
        zs = offline_quadruplet(b, RIEMANN_ZEROS[0]) + on[1:]
        res = envelope_slope(zs, *y_range)
        out["off_line"][b] = {"slope": res["slope"], "predicted": 0.5 - b,
                              "matches": abs(res["slope"] - (0.5 - b)) < 0.05}
    out["consistent_with_RH"] = out["on_line"]["bounded"]
    return out


# ═══════════════════════════════════════════════════════════════
#  Speiser–Hadamard Lab: Zeros of ζ′ (audit Addendum I, July 2026)
# ═══════════════════════════════════════════════════════════════
#
# The Dirichlet ensemble of the three-ensemble program is anchored on
#
#     Speiser (1935):        RH ⟺ ζ′ has no zeros in 0 < Re s < 1/2
#     Hadamard bookkeeping:  ∂²log(−ζ′(s)) = 2/(s−1)² − Σ_{ρ′} (s−ρ′)⁻²
#
# where ρ′ runs over ALL zeros of ζ′: the nontrivial zeros in the strip
# (exactly 19 with 0 < Im ≤ 100, certified below) and the negative real
# zeros, exactly one in each interval (−2n−2, −2n), n ≥ 1.
#
# Independent verification status (mtft audit, Addendum I):
#   * algebra — (s−1)²ζ′(s) is entire of order 1 and Σ|ρ′|⁻² < ∞, so the
#     canonical product has genus ≤ 1 and ∂²log kills both the e^{A+Bs}
#     factor and the compensators; the displayed identity is a theorem.
#   * pole coefficient — ∂²log(−ζ′(s))·(s−1)² → 2.00000 as s → 1⁺,
#     computed from raw ζ derivatives (no Hadamard input).
#   * negative zeros — exact functional-equation solver (see
#     zetaprime_negative_zero), bracket-certified, n = 1..1000.
#   * census — the 19 nontrivial zeros below, located by argument-
#     principle bisection and Newton-refined to |ζ′(ρ′)| < 1e-29; count
#     consistent with Berndt (1970): N′(T) = (T/2π)log(T/4πe) + O(log T).
#     All 19 have Re > 1/2 — a numerical Speiser check (RH-consistent)
#     to height 100.
#   * identity — balances to |residual| < 1.5e-5 for s ∈ {3,...,30} with
#     both tails carried (negative-axis beyond the computed zeros, and
#     high-γ nontrivial tail via the Berndt density).  The two tails
#     nearly CANCEL at large s (+0.0067 vs −0.0092 at s = 30): dropping
#     either one produces a spurious "missing zeros" signal.  Because
#     LHS decays like (2/3)^s while the pole and zero-sums (~1e-3)
#     cancel, numerical demonstrations must use s ∈ [3, 10].

ZETAPRIME_ZEROS: List[complex] = [
    complex(2.4631618694543213, 23.298320492762858),
    complex(1.2864968222690477, 31.708250083115909),
    complex(2.3075700637226316, 38.489983173078936),
    complex(1.3827636057116746, 42.29096455459673),
    complex(0.96468562270568565, 48.847159905068479),
    complex(2.1016999009487748, 52.432161245149836),
    complex(1.8959597624712398, 57.134753199019534),
    complex(0.84873532810540347, 60.140845782038424),
    complex(1.207295624674169, 65.919932824281162),
    complex(1.8329479316538901, 68.611078827128335),
    complex(1.7742690858377651, 71.528161065185035),
    complex(0.8646228644261133, 76.362807896467042),
    complex(1.3285155423330835, 78.662405942406661),
    complex(1.2035601348826901, 83.66913350341483),
    complex(2.3940392808396954, 85.802080034941309),
    complex(0.8641036405989395, 88.177517409881013),
    complex(1.3040878149430769, 93.085926815619881),
    complex(0.78062800472464465, 95.292968271352217),
    complex(1.7984373897654074, 98.826971867454158),
]
"""The certified census: all nontrivial zeros of ζ′ with 0 < Im s ≤ 100.

19 zeros (argument-principle certified); each Newton-refined to
|ζ′(ρ′)| < 1e-29.  Minimum real part 0.78062800... > 1/2 (Speiser).
Conjugates are omitted; ζ′(ρ̄′) = 0 as well.
"""

ZETAPRIME_CENSUS_HEIGHT = 100.0


def zetaprime_zero_count_berndt(T: float) -> float:
    """
    Leading term of the zero-counting function of ζ′ (Berndt 1970):

        N′(T) = (T/2π)·log(T/4πe) + O(log T)

    Note the log(T/4π): the ζ′ ensemble is SPARSER than the ζ ensemble
    by ≈ (T/2π)·log 2.  At T = 100: 17.1 + O(log T) — the certified
    census count is 19.  (For ζ the formula reads (T/2π)log(T/2πe).)
    """
    mp = _mp()
    return float(T / (2 * mp.pi) * mp.log(T / (4 * mp.pi * mp.e)))


def _logzeta_deriv(u):
    """(log ζ)′(u) for real u ≥ 4 (prime-power expansion; exact below 12)."""
    mp = _mp()
    if u < 12:
        return mp.zeta(u, derivative=1) / mp.zeta(u)
    two, three = mp.mpf(2), mp.mpf(3)
    t = (mp.log(2) * two ** (-u) + mp.log(3) * three ** (-u)
         + mp.log(2) * mp.mpf(4) ** (-u) + mp.log(5) * mp.mpf(5) ** (-u)
         + mp.log(7) * mp.mpf(7) ** (-u) + mp.log(2) * mp.mpf(8) ** (-u)
         + mp.log(3) * mp.mpf(9) ** (-u) + mp.log(11) * mp.mpf(11) ** (-u)
         + mp.log(13) * mp.mpf(13) ** (-u) + mp.log(2) * mp.mpf(16) ** (-u))
    return -t


def _H2(s):
    """
    The exact negative-axis zero equation for ζ′.

    From the functional equation ζ(s) = χ(s)ζ(1−s),

        ζ′(s) = χ(s)ζ(1−s) · [ (log χ)′(s) − (log ζ)′(1−s) ],

    and χ(s)ζ(1−s) ≠ 0 inside each interval (−2n−2, −2n), so the zeros
    of ζ′ there are exactly the zeros of

        H2(s) = log(2π) + (π/2)·cot(πs/2) − ψ(1−s) − (log ζ)′(1−s).

    Elementary (digamma + cotangent); no ζ evaluations at negative
    arguments — this is what makes the deep negative axis reachable.
    """
    mp = _mp()
    return (mp.log(2 * mp.pi) + (mp.pi / 2) / mp.tan(mp.pi * s / 2)
            - mp.digamma(1 - s) - _logzeta_deriv(1 - s))


_NEG_ZERO_CACHE: dict = {}


def zetaprime_negative_zero(n: int) -> float:
    """
    The unique real zero of ζ′ in (−2n−2, −2n), n ≥ 1 (bisection on H2).

    Verified anchors: ρ′₁ = −2.7172628292, ρ′₂ = −4.9367621086,
    ρ′₃ = −7.0745971450; computed to n = 1000 in the audit
    (ρ′₁₀₀₀ = −2001.83062775).  The zeros approach −2n−1 from the left
    with a slow drift (ε ≈ −0.43 at n = 10, −0.83 at n = 1000); the
    naive asymptotic ε ≈ −(4/π²)log(n/π) overshoots — solve, don't
    linearize.
    """
    if n < 1:
        raise ValueError("negative zeros of zeta' are indexed from n = 1")
    if n in _NEG_ZERO_CACHE:
        return _NEG_ZERO_CACHE[n]
    mp = _mp()
    with mp.workdps(25):
        lo, hi = mp.mpf(-2 * n - 2) + mp.mpf("1e-9"), mp.mpf(-2 * n) - mp.mpf("1e-9")
        flo = _H2(lo)
        if (flo < 0) == (_H2(hi) < 0):  # pragma: no cover
            raise RuntimeError(f"H2 bracket failed at n = {n}")
        for _ in range(70):
            mid = (lo + hi) / 2
            fm = _H2(mid)
            if fm == 0:
                lo = hi = mid
                break
            if (flo < 0) != (fm < 0):
                hi = mid
            else:
                lo, flo = mid, fm
        r = (lo + hi) / 2
    _NEG_ZERO_CACHE[n] = float(r)
    return _NEG_ZERO_CACHE[n]


def zetaprime_refine(z0, tol: float = 1e-16) -> complex:
    """
    Newton-refine an approximate zero of ζ′ (complex).

    Hand-rolled Newton rather than mpmath.findroot: long findroot sweeps
    on ζ′ are crash-prone in some environments (audit experience).
    """
    mp = _mp()
    with mp.workdps(30):
        z = mp.mpc(z0)
        h = mp.mpc("1e-12", "0")
        for _ in range(60):
            f0 = mp.zeta(z, derivative=1)
            d = (mp.zeta(z + h, derivative=1) - f0) / h
            if d == 0:
                break
            step = f0 / d
            z -= step
            if abs(step) < mp.mpf(10) ** int(mp.floor(mp.log10(tol))):
                break
    return complex(z)


def zetaprime_logcurvature(s: float) -> float:
    """
    ∂²/∂s² log(−ζ′(s)) = ζ‴/ζ′ − (ζ″/ζ′)²   (real s > 1).

    Anchors: 0.014052477562801085548 at s = 8; 1.3597768711309190e-6
    at s = 30.  Pole check: value·(s−1)² → 2 as s → 1⁺.
    """
    mp = _mp()
    with mp.workdps(30):
        s = mp.mpf(s)
        z1 = mp.zeta(s, derivative=1)
        z2 = mp.zeta(s, derivative=2)
        z3 = mp.zeta(s, derivative=3)
        return float(z3 / z1 - (z2 / z1) ** 2)


def _nontrivial_pair_sum(s):
    """Σ over certified census, conjugate pairs folded: 2Re Σ (s−ρ′)⁻²."""
    mp = _mp()
    tot = mp.mpf(0)
    for r in ZETAPRIME_ZEROS:
        b, g = mp.mpf(r.real), mp.mpf(r.imag)
        tot += 2 * ((s - b) ** 2 - g ** 2) / (((s - b) ** 2 + g ** 2) ** 2)
    return tot


def _berndt_density(t):
    mp = _mp()
    return mp.log(t / (4 * mp.pi)) / (2 * mp.pi)


def _nontrivial_tail(s, g0=100.0, bbar=1.5):
    """High-γ tail via the Berndt density (uncertainty ~ ±2e-4: one pair)."""
    mp = _mp()
    return mp.quad(
        lambda t: 2 * ((s - bbar) ** 2 - t ** 2)
        / (((s - bbar) ** 2 + t ** 2) ** 2) * _berndt_density(t),
        [mp.mpf(g0), mp.inf])


def hadamard_zetaprime_check(s: float, n_neg: int = 120) -> dict:
    """
    Numerical evaluation of the Speiser–Hadamard identity

        ∂²log(−ζ′(s))  =  2/(s−1)² − Σ_{ρ′} (s−ρ′)⁻²

    with the zero sum split as: n_neg exact negative-axis zeros +
    integral tail beyond, the certified 19-zero nontrivial census, and
    the high-γ Berndt-density tail.

    Returns lhs, pole, negative_axis, nontrivial_census, nontrivial_tail,
    rhs, residual.  Expected |residual| ~ 1e-5 for s ∈ [3, 10], far below
    the tail model's ~ ±2e-4 discreteness uncertainty.  CONDITIONING
    WARNING: for s ≳ 20 the LHS falls below the tail-accounting noise
    floor (it decays like (2/3)^s while ~1e-3 terms cancel); use
    s ∈ [3, 10] for demonstrations.
    """
    mp = _mp()
    with mp.workdps(25):
        s = mp.mpf(s)
        lhs = mp.mpf(zetaprime_logcurvature(float(s)))
        pole = 2 / (s - 1) ** 2
        neg = mp.mpf(0)
        for n in range(1, n_neg + 1):
            neg += 1 / (s - mp.mpf(zetaprime_negative_zero(n))) ** 2
        neg += mp.quad(lambda nn: 1 / (s + 2 * nn + 1) ** 2,
                       [n_neg + mp.mpf("0.5"), mp.inf])
        nt = _nontrivial_pair_sum(s)
        nt_tail = _nontrivial_tail(s)
        rhs = pole - neg - nt - nt_tail
    return {
        "s": float(s),
        "lhs": float(lhs),
        "pole": float(pole),
        "negative_axis": float(neg),
        "nontrivial_census": float(nt),
        "nontrivial_tail": float(nt_tail),
        "rhs": float(rhs),
        "residual": float(lhs - rhs),
        "identity": "d^2 log(-zeta'(s)) = 2/(s-1)^2 - sum_{rho'} (s-rho')^{-2}",
        "note": "use 3 <= s <= 10 for demonstrations (conditioning)",
    }


# ── Decomposition lemma (three-ensemble paper §2) ────────────────

def dirichlet_curvature(beta: float) -> dict:
    """
    The Dirichlet-ensemble curvature and its decomposition lemma:

        g_D(β) = ∂²log ζ(β) + ∂²log(−ζ′(β+1)),     β > 1.

    Verified anchors: g_D(3) = 0.33510387864414189, with the ζ′ piece
    carrying 48.588864% at β = 3 (41.92% at β = 2.5, 58.04% at β = 4).
    The ζ piece alone has the von Mangoldt series
    ∂²log ζ(β) = Σ_n Λ(n)(log n) n^{−β} (see von_mangoldt_curvature).
    """
    mp = _mp()
    with mp.workdps(30):
        b = mp.mpf(beta)
        zeta_piece = (mp.zeta(b, derivative=2) / mp.zeta(b)
                      - (mp.zeta(b, derivative=1) / mp.zeta(b)) ** 2)
        zp_piece = mp.mpf(zetaprime_logcurvature(float(b) + 1.0))
    g = zeta_piece + zp_piece
    return {
        "beta": float(beta),
        "g_D": float(g),
        "zeta_piece": float(zeta_piece),
        "zetaprime_piece": float(zp_piece),
        "zetaprime_share": float(zp_piece / g),
        "identity": "g_D(beta) = d^2 log zeta(beta) + d^2 log(-zeta'(beta+1))",
    }


def von_mangoldt_curvature(beta: float, n_max: int = 20000) -> float:
    """
    Σ_{n≤n_max} Λ(n)(log n) n^{−β} — the von Mangoldt series for
    ∂²log ζ(β) (cross-check of dirichlet_curvature's zeta piece;
    truncation error ~ (log n_max)²/(2·n_max^{β−1}·(β−1))).
    """
    import math as _math
    lam = [0.0] * (n_max + 1)
    for p in range(2, n_max + 1):
        if lam[p] == 0.0 and all(p % q for q in range(2, int(p ** 0.5) + 1)):
            lp = _math.log(p)
            pk = p
            while pk <= n_max:
                lam[pk] = lp
                pk *= p
    return sum(l * _math.log(n) * n ** (-beta)
               for n, l in enumerate(lam) if l != 0.0 and n >= 2)


# ── Weighted theta: exact shift identity + modularity no-go ─────

def divisor_log_weights(n_max: int) -> list:
    """
    w_n = Σ_{d|n} (log d)/d for n = 1..n_max (the "Emergent weights";
    their Dirichlet series is W(s) = −ζ(s)ζ′(s+1), since
    −ζ′(s+1) = Σ (log n)/n · n^{−s} convolved with ζ(s)).
    """
    import math as _math
    w = [0.0] * (n_max + 1)
    for d in range(2, n_max + 1):
        c = _math.log(d) / d
        for n in range(d, n_max + 1, d):
            w[n] += c
    return w


def weighted_theta(y: float, weights: Optional[list] = None,
                   n_max: int = 200000) -> float:
    """Θ̃(y) = Σ_n w_n e^{−2πny} with the divisor-log weights (numpy)."""
    if weights is None:
        weights = divisor_log_weights(n_max)
    n_max = len(weights) - 1
    n = np.arange(1, n_max + 1, dtype=float)
    return float(np.dot(np.array(weights[1:], dtype=float),
                        np.exp(-2.0 * math.pi * y * n)))


def filtered_moment_identity(y: float, N: int = 3,
                             n_max: int = 2000) -> dict:
    """
    Exact shift identity for the N-filtered second moment of the
    weighted theta (audit, Emergent-analysis sharpening):

        μ_N(y) := Σ_n n² w_n e^{−2πny} (1 − cos(2πn/N))
                = (1/4π²)·[ T″(y) − Re T″(y − i/N) ],   T := Θ̃.

    Holds term-by-term (Re e^{2πin/N} = cos(2πn/N)); verified to better
    than 1e-12 relative.  For N = 3 the factor (1 − cos(2πn/3)) is 0 on
    3|n and 3/2 otherwise — i.e. it IS the SU(3) center projector
    (3∤n filter of mass_gap_stiffness) up to the factor 3/2.
    """
    weights = divisor_log_weights(n_max)
    n = np.arange(1, n_max + 1, dtype=float)
    w = np.array(weights[1:], dtype=float)
    e = np.exp(-2.0 * math.pi * y * n)
    mu_direct = float(np.dot(n * n * w * e,
                             1.0 - np.cos(2.0 * math.pi * n / N)))
    # analytic T″:  T″(z) = Σ (2πn)² w_n e^{−2πnz}
    z = complex(y, -1.0 / N)
    es = np.exp(-2.0 * math.pi * z * n)
    tpp = complex(np.dot((2.0 * math.pi * n) ** 2 * w, es))
    tpp0 = float(np.dot((2.0 * math.pi * n) ** 2 * w, e))
    rhs = (tpp0 - tpp.real) / (4.0 * math.pi ** 2)
    return {
        "y": y, "N": N,
        "mu_direct": mu_direct,
        "shift_formula": rhs,
        "rel_diff": abs(mu_direct - rhs) / abs(mu_direct) if mu_direct else 0.0,
        "identity": "mu_N(y) = (1/4 pi^2) [T''(y) - Re T''(y - i/N)]",
    }


def weighted_theta_cusp_fit(y_points=(1e-2, 1e-3, 1e-4),
                            n_max: int = 200000) -> dict:
    """
    Modularity no-go for the weighted theta (audit, Emergent analysis).

    The Mellin parent W(s) = −ζ(s)ζ′(s+1) has a DOUBLE pole at s = 0
    (ζ′(s+1) ~ −1/s², ζ(0) = −1/2), and Γ(s) adds a third order, so
    Mellin inversion forces a genuine logarithmic cusp term:

        Θ̃(y) = (−ζ′(2))/X − (1/4)·ln²(1/X) + B·ln(1/X) + C + o(1),
        X = 2πy,

    with coefficient exactly ζ(0)·(−1)·(1/2) = −1/4.  Modular forms have
    pure-power cusp asymptotics, so Θ̃ CANNOT be modular — that route is
    closed analytically, no positivity assumption needed.

    Fits A in Θ̃ − (−ζ′(2))/X = A·L² + B·L + C over y_points and
    compares with −1/4.  Audit fit over y = 1e-2..1e-5:
    A = −0.2498 (residual quadratic after removing −L²/4: +2e-4);
    C ≈ −0.919 ≈ ζ′(0) = −½ln(2π).
    """
    mp = _mp()
    weights = divisor_log_weights(n_max)
    lead = float(-mp.zeta(2, derivative=1))  # 0.9375482543...
    pts = []
    for y in y_points:
        X = 2.0 * math.pi * y
        L = math.log(1.0 / X)
        th = weighted_theta(y, weights)
        pts.append((L, th - lead / X))
    if len(pts) == 3:
        M = np.array([[L * L, L, 1.0] for L, _ in pts])
        A, B, C = np.linalg.solve(M, np.array([v for _, v in pts]))
    else:
        A, B, C = np.polyfit(np.array([p[0] for p in pts]),
                             np.array([p[1] for p in pts]), 2)
    return {
        "A_fit": float(A),
        "A_predicted": -0.25,
        "B_fit": float(B),
        "C_fit": float(C),
        "C_predicted_approx": float(mp.zeta(0, derivative=1)),  # -0.9189385332
        "y_points": list(y_points),
        "verdict": ("double pole at s=0 forces -1/4 * ln^2(1/X); "
                    "pure-power cusp asymptotics impossible -> not modular"),
    }
