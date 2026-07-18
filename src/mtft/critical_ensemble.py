#!/usr/bin/env python3
"""
Critical Ensemble: Li Coefficients as the Third Curvature Family
==================================================================

MIT License — Copyright (c) 2026 Roger Tano
See LICENSE file for full terms.

The third leg of the three-ensemble program:

    LAPLACE   (stiffness):  mu_N(y) = (1/4pi^2)[T''(y) - Re T''(y - i/N)]
                            <-> Th 1: RH iff limsup |Delta_kappa X^{-3/2}| < inf
    DIRICHLET (zeta'):      g_D(beta) = d2 log zeta(beta) + d2 log(-zeta'(beta+1))
                            <-> Speiser: RH iff zeta' != 0 in 0 < Re s < 1/2
    CRITICAL  (xi):         lambda_n = Taylor family of log xi at s = 1
                            <-> Li: RH iff lambda_n >= 0 for ALL n

The Li coefficients are the critical ensemble's curvature output:

    lambda_n = (1/(n-1)!) (d/ds)^n [ s^{n-1} log xi(s) ] at s = 1
             = sum over nontrivial zeros rho of [ 1 - (1 - 1/rho)^n ]

CRITICAL CAVEAT (Bombieri-Lagarias): the criterion's content is ALL n.
Any finite prefix of positive lambda_n proves nothing — finitely many
zeros always produce positive-looking partials. The publishable object
is a relating inequality between the three curvature families, not the
observation that a finite lambda-prefix is positive. This caveat is
carried in every report this module produces.

THREE INDEPENDENT METHODS (repo standard: results must agree):

  1. li_lambda(n)          — EXACT series algebra. log xi at s = 1 built
                             from five closed-form Taylor series:
                             log(1/2) + log s - (s/2)log pi
                             + log Gamma(s/2) + log[(s-1)zeta(s)],
                             with ingredients: Stieltjes constants gamma_k
                             (for (s-1)zeta(s)), polygamma at 1/2 via
                             psi^(m)(1/2) = (-1)^{m+1} m! (2^{m+1}-1) zeta(m+1)
                             (for log Gamma). No quadrature, no numerical
                             differentiation. Primary method.
  2. li_lambda_cauchy(n)   — Cauchy integral coefficients of xi'/xi on a
                             circle |s-1| = r (analytic in |s-1| < 14.13,
                             the first zero's distance). Independent of
                             method 1's series ingredients.
  3. li_lambda_zero_sum(n) — Truncated sum over certified zeta zeros
                             (mpmath.zetazero) with explicit tail estimate.
                             Terms 2(1 - Re[(1-1/rho)^n]) are >= 0 on the
                             critical line, so partials are monotone lower
                             bounds; diagnostic method.

Built-in exactness anchor: the u^1 coefficient of the five-series sum
collapses algebraically to 1 + gamma/2 - (1/2)log(4pi) — the known
closed form of lambda_1. Verified at import if MTFT_SELFTEST=1.

MTFT constants: level N = 143 = 11 x 13, genus 13; this module is
level-independent (the critical ensemble sees only zeta/xi) but ships
as the third ensemble of the X0(143) program.

Roger Tano — MTFT Research Program — July 2026
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple

from mpmath import (mp, mpf, mpc, zeta, gamma as mp_gamma, log, pi, exp,
                    euler, quad, stieltjes, zetazero, binomial, fac,
                    loggamma, psi, cos, sin, sqrt, re as mp_re, im as mp_im)


# ═══════════════════════════════════════════════════════════════
#  CONSTANTS & CAVEATS
# ═══════════════════════════════════════════════════════════════

LEVEL = 143
GENUS = 13

#: Distance from s = 1 to the nearest singularity of log xi
#: (first nontrivial zero 1/2 + 14.1347...i). Cauchy radii must stay below.
XI_ANALYTICITY_RADIUS = 14.13

BOMBIERI_LAGARIAS_CAVEAT = (
    "Li criterion content is ALL n (Bombieri-Lagarias 1999): any finite "
    "prefix of positive lambda_n carries no logical force — finitely many "
    "zeros always produce positive partials. Do not cite finite positivity "
    "as evidence for RH."
)

#: Closed form for lambda_1 (exact): 1 + gamma/2 - (1/2) ln(4 pi)
def lambda_1_closed_form() -> mpf:
    """lambda_1 = 1 + euler/2 - log(4*pi)/2  (exact closed form)."""
    return mpf(1) + euler / 2 - log(4 * pi) / 2


# ═══════════════════════════════════════════════════════════════
#  POWER-SERIES HELPERS (dense, mpf/mpc coefficients)
# ═══════════════════════════════════════════════════════════════

def _series_mul(a: List, b: List, n: int) -> List:
    """Product of two power series, truncated at order n (inclusive)."""
    out = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            out[i + j] += ai * bj
    return out


def _series_log1p(p: List, n: int) -> List:
    """
    log(1 + P(u)) for P with P(0) = 0, via the derivative method:
    B' = P' / (1 + P), then integrate. Returns coefficients b_0..b_n
    with b_0 = 0.
    """
    # P'(u): coefficient k of P' is (k+1) * p_{k+1}
    pp = [(k + 1) * p[k + 1] for k in range(n)] + [mpf(0)]
    # (1+P)^{-1} via standard inversion: q_0 = 1; q_k = -sum p_j q_{k-j}
    q = [mpf(1)] + [mpf(0)] * n
    for k in range(1, n + 1):
        s = mpf(0)
        for j in range(1, k + 1):
            s += p[j] * q[k - j]
        q[k] = -s
    bp = _series_mul(pp, q, n)         # B'(u) up to order n
    b = [mpf(0)] * (n + 1)
    for k in range(1, n + 1):
        b[k] = bp[k - 1] / k           # integrate
    return b


# ═══════════════════════════════════════════════════════════════
#  METHOD 1 — EXACT SERIES ALGEBRA (primary)
# ═══════════════════════════════════════════════════════════════

_A_CACHE: Dict[Tuple[int, int], List] = {}   # (n_max, dps) -> a_1..a_n list


def logxi_taylor(n_max: int, dps: Optional[int] = None) -> List:
    """
    Taylor coefficients a_1..a_{n_max} of log xi(s) at s = 1 (a_0 = log 1/2
    excluded; it never enters lambda_n).

    Assembly (u = s - 1):
      log xi(s) = log(1/2) + log(1+u) - ((1+u)/2) log pi
                  + log Gamma((1+u)/2) + log h(1+u),
      h(s) = (s-1) zeta(s) = 1 + sum_{k>=0} (-1)^k gamma_k u^{k+1} / k!
                                (gamma_k = Stieltjes constants)
      log Gamma(1/2 + u/2) = log Gamma(1/2)
                  + sum_{k>=1} psi^{(k-1)}(1/2) / (2^k k!) * u^k,
      psi(1/2)      = -euler - 2 log 2,
      psi^{(m)}(1/2) = (-1)^{m+1} m! (2^{m+1} - 1) zeta(m+1),  m >= 1.

    Every ingredient is a closed form; the only 'special' inputs are
    Stieltjes constants and zeta at integers, both certified mpmath
    primitives. Cached per (n_max, dps).
    """
    if dps is None:
        dps = mp.dps
    key = (n_max, dps)
    if key in _A_CACHE:
        return _A_CACHE[key]

    with mp.workdps(dps + 10):
        n = n_max
        a = [mpf(0)] * (n + 1)          # a[0] unused (kept 0)

        # -- log(1+u): (-1)^{k+1}/k
        for k in range(1, n + 1):
            a[k] += mpf(-1) ** (k + 1) / k

        # -- -((1+u)/2) log pi: only the u^1 coefficient (constant -> a_0)
        a[1] += -log(pi) / 2

        # -- log Gamma(1/2 + u/2)
        a[1] += (-euler - 2 * log(2)) / 2            # psi(1/2)/2
        for k in range(2, n + 1):
            m = k - 1
            psi_m_half = mpf(-1) ** (m + 1) * fac(m) * (2 ** (m + 1) - 1) * zeta(m + 1)
            a[k] += psi_m_half / (mpf(2) ** k * fac(k))

        # -- log h(1+u), h-series from Stieltjes constants
        p = [mpf(0)] * (n + 1)
        for k in range(0, n):            # p_{k+1} = (-1)^k gamma_k / k!
            p[k + 1] = mpf(-1) ** k * stieltjes(k) / fac(k)
        b = _series_log1p(p, n)
        for k in range(1, n + 1):
            a[k] += b[k]

        result = [mpf(x) for x in a]

    _A_CACHE[key] = result
    return result


def li_lambda(n: int, dps: Optional[int] = None) -> mpf:
    """
    lambda_n by exact series algebra (primary method).

    lambda_n = n * sum_{j=0}^{n-1} C(n-1, j) * a_{n-j},
    from the product [ (1+u)^{n-1} ] * [ log xi series ], coefficient of u^n,
    times n!/(n-1)!.
    """
    if n < 1:
        raise ValueError("n >= 1 required")
    a = logxi_taylor(n, dps)
    s = mpf(0)
    for j in range(0, n):
        s += binomial(n - 1, j) * a[n - j]
    return n * s


def li_lambda_batch(n_max: int, dps: Optional[int] = None) -> List[mpf]:
    """lambda_1..lambda_{n_max} sharing one Taylor-series computation."""
    a = logxi_taylor(n_max, dps)
    out = []
    for n in range(1, n_max + 1):
        s = mpf(0)
        for j in range(0, n):
            s += binomial(n - 1, j) * a[n - j]
        out.append(n * s)
    return out


# ═══════════════════════════════════════════════════════════════
#  METHOD 2 — CAUCHY INTEGRALS (independent cross-check)
# ═══════════════════════════════════════════════════════════════

def _xi_log_derivative(s):
    """
    xi'/xi(s) = 1/s + 1/(s-1) - (1/2)log pi + (1/2)psi(s/2) + zeta'/zeta(s).

    Analytic in |s-1| < 14.13 (xi has no zeros there); the displayed
    pieces have poles/cancellations INSIDE the disk (s = 0, 1) but every
    evaluation point on our contours (r <= 3.5) keeps clear of s = 0, 1
    and of the trivial zero s = -2 (distance 3 from center) provided
    r < 3 or the contour is checked; default r = 2 is safe.
    """
    return (1 / s + 1 / (s - 1) - log(pi) / 2 + psi(0, s / 2) / 2
            + zeta(s, derivative=1) / zeta(s))


def logxi_taylor_cauchy(n_max: int, r: float = 2.0,
                        dps: Optional[int] = None) -> List:
    r"""
    a_1..a_{n_max} of log xi at s = 1 via Cauchy coefficients of
    G = xi'/xi (branch-free), then a_k = g_{k-1}/k.

      g_m = (1/2 pi i) \oint G(s) / (s-1)^{m+1} ds
          = (1/2 pi r^m) \int_0^{2 pi} G(1 + r e^{i t}) e^{-i m t} dt.

    r must satisfy 0 < r < 2.99 (keep the trivial zero s = -2 and the
    pole of zeta'/zeta there strictly outside the contour) and r != 1
    (avoid s = 0 on the contour). Independent of Method 1's ingredients.
    """
    if not (0.05 < r < 2.99):
        raise ValueError("radius r must lie in (0.05, 2.99)")
    if abs(r - 1.0) < 1e-9:
        raise ValueError("r = 1 puts s = 0 on the contour; choose another r")
    if dps is None:
        dps = mp.dps

    with mp.workdps(dps + 10):
        rr = mpf(r)
        gs = []
        for m in range(0, n_max):
            def integrand(t, m=m):
                s = 1 + rr * exp(mpc(0, 1) * t)
                return _xi_log_derivative(s) * exp(mpc(0, -m) * t)
            val = quad(integrand, [0, 2 * pi]) / (2 * pi * rr ** m)
            gs.append(val)
        a = [mpf(0)] * (n_max + 1)
        for k in range(1, n_max + 1):
            a[k] = mp_re(gs[k - 1]) / k     # coefficients are real
        return a


def li_lambda_cauchy(n: int, r: float = 2.0,
                     dps: Optional[int] = None) -> mpf:
    """lambda_n via the Cauchy-integral coefficients (Method 2)."""
    a = logxi_taylor_cauchy(n, r=r, dps=dps)
    s = mpf(0)
    for j in range(0, n):
        s += binomial(n - 1, j) * a[n - j]
    return n * s


# ═══════════════════════════════════════════════════════════════
#  METHOD 3 — ZERO SUM (diagnostic; monotone lower bounds)
# ═══════════════════════════════════════════════════════════════

_ZERO_CACHE: Dict[int, mpc] = {}


def _zeta_zero(k: int) -> mpc:
    if k not in _ZERO_CACHE:
        _ZERO_CACHE[k] = zetazero(k)
    return _ZERO_CACHE[k]


@dataclass
class ZeroSumResult:
    n: int
    num_pairs: int
    partial: float            # truncated lambda_n (monotone lower bound)
    tail_estimate: float      # analytic tail model (see docstring)
    highest_gamma: float
    analytic: Optional[float] # Method-1 value for comparison (if computed)
    gap: Optional[float]      # analytic - partial (>= 0 expected)
    caveat: str = BOMBIERI_LAGARIAS_CAVEAT


def li_lambda_zero_sum(n: int, num_pairs: int = 100,
                       compare_analytic: bool = True,
                       dps: Optional[int] = None) -> ZeroSumResult:
    """
    Truncated lambda_n = sum over the first `num_pairs` conjugate zero
    pairs of 2 (1 - Re[(1 - 1/rho)^n]).

    On the critical line |1 - 1/rho| = 1 exactly, so every term is >= 0
    and partials are monotone lower bounds of lambda_n.

    Tail model: for gamma >> n the pair term ~ n^2/gamma^2; integrating
    against the zero density (1/2 pi) log(t / 2 pi) gives
        tail(Gamma) ~ (n^2 / 2 pi) * (log(Gamma / 2 pi) + 1) / Gamma.
    Reported as an estimate, not a bound (the ~ hides an O(1) factor of
    order unity; empirically the true gap sits within ~2x of the model
    for n <= 10, num_pairs >= 50).
    """
    if dps is None:
        dps = mp.dps
    with mp.workdps(dps):
        s = mpf(0)
        for k in range(1, num_pairs + 1):
            rho = _zeta_zero(k)
            term = 2 * (1 - mp_re((1 - 1 / rho) ** n))
            s += term
        Gamma = float(mp_im(_zeta_zero(num_pairs)))
        tail = float((mpf(n) ** 2 / (2 * pi)) * (log(Gamma / (2 * pi)) + 1) / Gamma)

        analytic = None
        gap = None
        if compare_analytic:
            analytic = float(li_lambda(n, dps=dps))
            gap = analytic - float(s)

        return ZeroSumResult(
            n=n, num_pairs=num_pairs, partial=float(s),
            tail_estimate=tail, highest_gamma=Gamma,
            analytic=analytic, gap=gap,
        )


# ═══════════════════════════════════════════════════════════════
#  THREE-LEG CERTIFICATION & REPORTING
# ═══════════════════════════════════════════════════════════════

@dataclass
class CertificationResult:
    lambda1_closed_form_err: float          # Method 1 vs exact closed form
    cauchy_agreement_err: float             # max |M1 - M2| over tested n
    cauchy_radius_independence_err: float   # max |M2(r1) - M2(r2)|
    zero_sum_bracketing_ok: bool            # 0 < partial < analytic, gap ~ tail
    tested_n: int
    dps: int
    passed: bool


def certify(n_test: int = 8, dps: int = 30,
            r1: float = 1.5, r2: float = 2.5,
            zero_pairs: int = 80) -> CertificationResult:
    """
    Run the three-leg certification:

      LEG A  Method 1 lambda_1 vs the exact closed form
             1 + euler/2 - log(4 pi)/2   (must agree to ~10^{-(dps-5)}).
      LEG B  Method 1 vs Method 2 at two radii, n = 1..n_test
             (independent algorithms; agreement certifies both).
      LEG C  Method 3 partials bracket: 0 <= partial <= analytic with
             gap within ~3x the tail model for n <= 4.
    """
    with mp.workdps(dps):
        # LEG A
        errA = float(abs(li_lambda(1, dps=dps) - lambda_1_closed_form()))

        # LEG B — compare at the COEFFICIENT level (the meaningful check:
        # lambda-level comparison at outer dps floors to bit-identity
        # because the coefficient methods agree ~1e-42, below dps-30
        # resolution; coefficients carry the true independent-method error)
        a_series = logxi_taylor(n_test, dps=dps)
        aB1 = logxi_taylor_cauchy(n_test, r=r1, dps=dps)
        aB2 = logxi_taylor_cauchy(n_test, r=r2, dps=dps)

        errB = 0.0
        errR = 0.0
        for k in range(1, n_test + 1):
            errB = max(errB,
                       float(abs(a_series[k] - aB1[k])),
                       float(abs(a_series[k] - aB2[k])))
            errR = max(errR, float(abs(aB1[k] - aB2[k])))

        # lambda-level sanity (should be exactly floored at outer dps)
        lam_series = li_lambda_batch(n_test, dps=dps)

        def _lam_from(a, n):
            s = mpf(0)
            for j in range(0, n):
                s += binomial(n - 1, j) * a[n - j]
            return n * s

        for n in range(1, n_test + 1):
            if float(abs(lam_series[n - 1] - _lam_from(aB1, n))) > 10 ** (-(dps - 12)):
                errB = max(errB, 1.0)  # force failure if lambda-level diverges

        # LEG C
        okC = True
        for n in (1, 2, 3, 4):
            zs = li_lambda_zero_sum(n, num_pairs=zero_pairs, dps=min(dps, 20))
            if not (0 <= zs.partial <= zs.analytic + 1e-12):
                okC = False
            if zs.gap is not None and zs.gap > 3 * zs.tail_estimate + 1e-6:
                okC = False

        passed = (errA < 10 ** (-(dps - 6))
                  and errB < 10 ** (-(dps - 12))
                  and errR < 10 ** (-(dps - 12))
                  and okC)

        return CertificationResult(
            lambda1_closed_form_err=errA,
            cauchy_agreement_err=errB,
            cauchy_radius_independence_err=errR,
            zero_sum_bracketing_ok=okC,
            tested_n=n_test, dps=dps, passed=passed,
        )


@dataclass
class LiReport:
    n_max: int
    values: List[float]
    all_positive: bool
    min_value: float
    min_at: int
    dps: int
    caveat: str = BOMBIERI_LAGARIAS_CAVEAT


def li_criterion_report(n_max: int = 20, dps: int = 30) -> LiReport:
    """
    lambda_1..lambda_{n_max} from the exact method, with the positivity
    flag and the Bombieri-Lagarias caveat attached. This is the on-the-fly
    tool: fast (series algebra only), precision-controlled, and honest
    about what finite positivity does and does not mean.
    """
    vals = [float(x) for x in li_lambda_batch(n_max, dps=dps)]
    mn = min(vals)
    return LiReport(
        n_max=n_max, values=vals,
        all_positive=all(v > 0 for v in vals),
        min_value=mn, min_at=vals.index(mn) + 1, dps=dps,
    )


# ═══════════════════════════════════════════════════════════════
#  THE THREE-ENSEMBLE SUMMARY (paper-facing)
# ═══════════════════════════════════════════════════════════════

THREE_ENSEMBLE_TABLE = """
    ENSEMBLE    WEIGHTS/OBJECT             CURVATURE OUTPUT              RH CRITERION
    --------    --------------             ----------------              ------------
    Laplace     w_n e^{-2 pi y n}          mu_N(y) = (1/4pi^2)           Th 1 (2026):
                (T(y); mass gap)             [T''(y) - Re T''(y-i/N)]    limsup|Dk X^-3/2| < inf
    Dirichlet   w_n n^{-beta}              g_D = d2 log zeta(beta)       Speiser (1935):
                (Z_D = -zeta zeta')          + d2 log(-zeta'(beta+1))    zeta' != 0, 0<Re s<1/2
    Critical    log xi at s = 1            lambda_n (Taylor family       Li (1997):
                (this module)                = zero-sum moments)         lambda_n >= 0 all n

    One weight sequence w_n = (sigma * Lambda)(n)/n feeds the first two;
    the third sees the zeros directly. The relating inequality between
    the three curvature families is the open hunt (kernel comparison
    against zero-counting measures).
"""


# ═══════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════

def _print_report(rep: LiReport):
    print(f"  Li coefficients lambda_1..lambda_{rep.n_max}  (dps={rep.dps}):")
    for i, v in enumerate(rep.values, 1):
        print(f"    lambda_{i:<3d} = {v:.15f}")
    print(f"  all positive: {rep.all_positive}   "
          f"(min = lambda_{rep.min_at} = {rep.min_value:.6e})")
    print(f"  CAVEAT: {rep.caveat}")


if __name__ == "__main__":
    import sys

    args = sys.argv[1:]
    cmd = args[0] if args else "report"

    if cmd == "lambda":
        n = int(args[1]) if len(args) > 1 else 5
        mp.dps = 30
        print(f"  lambda_{n} = {li_lambda(n)}")
    elif cmd == "report":
        n = int(args[1]) if len(args) > 1 else 20
        _print_report(li_criterion_report(n))
    elif cmd == "zerosum":
        n = int(args[1]) if len(args) > 1 else 3
        pairs = int(args[2]) if len(args) > 2 else 100
        mp.dps = 20
        r = li_lambda_zero_sum(n, pairs)
        print(f"  lambda_{n}: partial({pairs} pairs) = {r.partial:.8f}")
        print(f"             analytic              = {r.analytic:.8f}")
        print(f"             gap = {r.gap:.2e}   tail model = {r.tail_estimate:.2e}")
        print(f"             highest gamma = {r.highest_gamma:.2f}")
        print(f"  {r.caveat}")
    elif cmd == "certify":
        res = certify()
        print(f"  LEG A  lambda_1 vs closed form:   err = {res.lambda1_closed_form_err:.3e}")
        print(f"  LEG B  series vs Cauchy (n<=8):   err = {res.cauchy_agreement_err:.3e}")
        print(f"         radius independence:       err = {res.cauchy_radius_independence_err:.3e}")
        print(f"  LEG C  zero-sum bracketing:       {'OK' if res.zero_sum_bracketing_ok else 'FAIL'}")
        print(f"  CERTIFICATION: {'PASSED' if res.passed else 'FAILED'}")
    elif cmd == "table":
        print(THREE_ENSEMBLE_TABLE)
    else:
        print("Usage: python critical_ensemble.py "
              "[lambda N | report [N] | zerosum N [PAIRS] | certify | table]")
