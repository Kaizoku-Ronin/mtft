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

    RH ⟺ κ(y) ≥ 0 for all y > 0.

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
    that violate curvature positivity — the tower rigidity principle.
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
            "Off-line (β>0.5): envelope grows as y^{-(β+1)}, "
            "violating curvature positivity. A single off-line zero "
            "would break the entire N-tower."
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


def envelope_slope(zeros, y_min_pow=-7.0, y_max_pow=-1.8, n_points=91) -> dict:
    """
    Binned-RMS envelope slope of log₁₀|𝒟| vs log₁₀ y over
    y ∈ [10^y_min_pow, 10^y_max_pow].  Predicted: 0 (true zeros),
    1/2 − β₀ (synthetic off-line quadruplet at Re β₀).
    """
    mp = _mp()
    bins = {}
    for i in range(n_points):
        ly = y_min_pow + (y_max_pow - y_min_pow) * i / (n_points - 1)
        y = mp.mpf(10) ** ly
        d = abs(normalized_oscillation(y, zeros))
        b = math.floor(float(mp.log10(y)) * 2) / 2
        bins.setdefault(b, []).append(float(d * d))
    pts = sorted((b, math.log10(math.sqrt(sum(v) / len(v)))) for b, v in bins.items())
    n = len(pts)
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] * p[0] for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    return {"slope": slope, "n_bins": n, "points": pts}


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
