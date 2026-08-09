#!/usr/bin/env python3
"""
w2_susceptibility.py — the Tano susceptibility and the (E, w) Fisher metric
============================================================================

MIT License — Copyright (c) 2026 Roger Tano

The 2024 Gibbs machine (mtft.combinatorial, v0.13.0) cranked one order
higher on the primon gas.  v0.13.0 landed the mean, <w>_beta =
-zeta'(beta+1).  This study lands the SECOND moment, the energy-weight
covariance, and the full Fisher metric of the two-parameter exponential
family  p_n ∝ exp(lambda w_n - beta log n).

CLAIMS AND CLASSES
------------------
I1 (Pr / EXACT).  For beta > 1,

    <w^2>_beta = T(beta),   sum_n w_n^2 n^{-s} = zeta(s) * T(s),

  with, writing sigma = s+1, u = 2s+2, v = s+2,

    T(s) = zeta''(v) C0 - 2 zeta'(v) C1 + zeta(v) C2,
    C0 = zeta(sigma)^2/zeta(u),
    C1 = zeta(sigma)^2 zeta'(u)/zeta(u)^2 - zeta(sigma) zeta'(sigma)/zeta(u),
    C2 = zeta(sigma)^2 (2 zeta'(u)^2/zeta(u)^3 - zeta''(u)/zeta(u)^2)
         - 2 zeta(sigma) zeta'(sigma) zeta'(u)/zeta(u)^2
         + zeta'(sigma)^2/zeta(u).

  Proof route: w_n^2 = sum_{d|n} sum_{e|n} (log d/d)(log e/e); summing
  n over multiples of lcm(d,e) factors out zeta(s); parametrize d = ga,
  e = gb with (a,b) = 1 so lcm = gab; the coprime sums C0, C1, C2
  evaluate by Moebius removal (sum mu(k) k^{-u} = 1/zeta(u) and its two
  u-derivatives); the g-sum contributes zeta(v), -zeta'(v), zeta''(v).
  Algebra gate below: brute lcm double sum agrees to 2e-10 (s=3).

I2 (Pr / EXACT).  Cov_beta(log n, w_n) = zeta''(beta+1).
  TWO DISJOINT PROOF ROUTES: (a) Gibbs calculus on the proved mean,
  d<w>/dbeta = -Cov(E, w), with <w> = -zeta'(beta+1); (b) direct
  Dirichlet series, sum w_n log n n^{-s} = -(d/ds)[-zeta(s)zeta'(s+1)]
  = zeta'(s)zeta'(s+1) + zeta(s)zeta''(s+1), divided by zeta and
  centered.  Both certified numerically below.

I3 (Pr).  HAGEDORN TRANSPARENCY OF THE WEIGHT SECTOR.  Every moment
  <w^k>_beta stays finite as beta -> 1+, although the primon gas itself
  diverges there (Z = zeta(beta) -> inf, <E> -> inf, Var(E) -> inf).
  Proof: the k-fold divisor expansion gives sum_n w_n^k n^{-s} =
  zeta(s) T_k(s) with T_k(s) = sum over k-tuples of
  prod_i (log d_i / d_i) * lcm(d_1..d_k)^{-s}; since
  lcm >= max >= (prod d_i)^{1/k},
      T_k(1) <= (sum_d (log d) d^{-1-1/k})^k = (-zeta'(1+1/k))^k < inf.
  The single zeta(s) prefactor always cancels the partition function,
  so <w^k>_{beta->1+} = T_k(1).  The weight observable does not see the
  Hagedorn wall.

I4 (Df + Pr + DIAGNOSTIC).  The two-parameter family
  Z(beta, lambda) = sum_n exp(lambda w_n) n^{-beta} has Fisher matrix

      g = [[ (log zeta)''(beta),   -zeta''(beta+1)          ],
           [ -zeta''(beta+1),       T(beta) - zeta'(beta+1)^2 ]],

  the covariance matrix of the sufficient statistics (-E, w).
  det g > 0 strictly for every beta > 1 (Pr: Cauchy-Schwarz is strict
  because w is not affine in log n — witness n = 1, 2, 3: affinity
  through n=1,2 forces w_n = (1/2) log n, but w_3 = (log 3)/3 != 
  (log 3)/2).  Scan below is confirmation (DIAGNOSTIC), including the
  near-cancellation of det g at leading 4^{-beta} order in the cold
  limit — the gas freezes onto {1, 2}, and two points are always
  affinely dependent, so strict positivity at large beta is carried by
  the n = 3 witness.

I5 (Pr + Cert).  DECORRELATION AT THE WALL, LOCK AT FREEZE-OUT.
  rho(beta) = Cov/sqrt(Var E Var w) = zeta''(beta+1) /
  sqrt((log zeta)''(beta) (T(beta) - zeta'(beta+1)^2)).
  As beta -> 1+:  (log zeta)'' ~ (beta-1)^{-2}, everything else finite,
  so rho ~ [zeta''(2)/sqrt(Var_cold)] * (beta-1) -> 0: at the Hagedorn
  wall the weight FORGETS the energy.  As beta -> inf, rho -> 1: the
  frozen two-state gas locks them.  Wall-rate coefficient
  zeta''(2)/sqrt(Var_cold) = 2.19175584954... (EXACT closed form).

COLD CONSTANTS (beta -> 1+ endpoints; EXACT, closed form, dps 30)
  <w^2>_cold = T(1)              = 1.70276979154901697001...
  chi_w,cold = T(1) - zeta'(2)^2 = 0.82377306237833093427...
  Dw_cold    = sqrt(chi_w,cold)  = 0.90761944799476996127...
  Cov_cold   = zeta''(2)         = 1.98928023429890102342...
  (Companions to <w>_cold = -zeta'(2) = 2*T_INF, the cold-gas alpha.)

PROPOSED NAME (Df): chi_w(beta) := Var_beta(w), "the Tano
susceptibility" — the fluctuation of the arithmetic weight in the
primon heat bath.

FILED, NOT CLAIMED (AG-D5): Dw_cold = 0.90762 sits near <w>_cold =
0.93755 (ratio 0.9681).  Proximity, no mechanism; the honest statement
is only that the cold weight remains order-unity relatively
fluctuating (non-self-averaging at the wall).  Dismissed as evidence.

GATES (all must pass; E2 route pairs share no computational steps)
  G1  algebra gate: T_closed vs brute lcm double sum (s = 3, 2.5)
  G2  D2 identity: divisor sieve of w^2 vs zeta(s)T(s), gap <= tail
      majorant (w_n <= (log n)^2/2 + 1 squared, integral comparison)
  G3  Cov series route: sieve of w log n vs zeta'zeta'_+ + zeta zeta''_+
  G4  Cov ensemble route: truncated Gibbs covariance vs zeta''(beta+1)
  G5  Gibbs consistency: numerical d<w>/dbeta + zeta''(beta+1) = 0
  G6  det g > 0 on the scan grid; rho in (0, 1) and increasing on grid
  G7  cold constants stable at dps 30 vs dps 45 (evaluation integrity)

Run:  python studies/w2_susceptibility.py
Writes w2_susceptibility_ledger.json next to itself.  ~40 s.

PROMOTION PATH (for a future minor bump into mtft.combinatorial):
  weight_second_moment_exact, weight_susceptibility_exact,
  weight_energy_covariance_exact, weight_fisher_matrix  (defined below,
  lift-ready; tests would mirror gates G2-G6).
"""

from __future__ import annotations

import json
import math
import os
import sys

import numpy as np
import mpmath as mp

mp.mp.dps = 30


def _z(s, k=0):
    return mp.zeta(s, derivative=k) if k else mp.zeta(s)


# ── promotion-ready closed forms ────────────────────────────────────

def weight_second_moment_exact(beta):
    """<w^2>_beta = T(beta).  I1, Pr/EXACT (beta > 1; finite at 1+)."""
    s = mp.mpf(beta)
    sig, u, v = s + 1, 2 * s + 2, s + 2
    zs, zps = _z(sig), _z(sig, 1)
    zu, zpu, zppu = _z(u), _z(u, 1), _z(u, 2)
    zv, zpv, zppv = _z(v), _z(v, 1), _z(v, 2)
    C0 = zs ** 2 / zu
    C1 = zs ** 2 * zpu / zu ** 2 - zs * zps / zu
    C2 = (zs ** 2 * (2 * zpu ** 2 / zu ** 3 - zppu / zu ** 2)
          - 2 * zs * zps * zpu / zu ** 2 + zps ** 2 / zu)
    return zppv * C0 - 2 * zpv * C1 + zv * C2


def weight_susceptibility_exact(beta):
    """chi_w(beta) = Var_beta(w) = T(beta) - zeta'(beta+1)^2.  Pr/EXACT."""
    return weight_second_moment_exact(beta) - _z(mp.mpf(beta) + 1, 1) ** 2


def weight_energy_covariance_exact(beta):
    """Cov_beta(log n, w) = zeta''(beta+1).  I2, Pr/EXACT."""
    return _z(mp.mpf(beta) + 1, 2)


def weight_fisher_matrix(beta):
    """Fisher matrix of the (beta, lambda) family at lambda = 0.  I4."""
    b = mp.mpf(beta)
    g11 = mp.diff(lambda x: mp.log(mp.zeta(x)), b, 2)
    g12 = -weight_energy_covariance_exact(b)
    g22 = weight_susceptibility_exact(b)
    return g11, g12, g22


# ── independent numeric routes ──────────────────────────────────────

def sieve_w(n_max: int) -> np.ndarray:
    w = np.zeros(n_max + 1)
    for d in range(2, n_max + 1):
        w[d::d] += math.log(d) / d
    return w


def brute_T(s: float, D: int) -> float:
    d = np.arange(1, D + 1)
    lcm = np.multiply.outer(d, d) // np.gcd.outer(d, d)
    ld = np.log(d)
    return float((np.outer(ld / d, ld / d) * lcm.astype(float) ** (-s)).sum())


def integral_tail(k: int, s: float, N: int) -> float:
    """int_N^inf x^{-s} (log x)^k dx, closed incomplete-gamma form."""
    L, a = math.log(N), s - 1.0
    tot, ff = 0.0, 1.0
    for j in range(k + 1):
        tot += ff * L ** (k - j) / a ** (j + 1)
        ff *= (k - j)
    return N ** (1 - s) * tot


def tail_w2(s: float, N: int) -> float:
    """Majorant for sum_{n>N} w_n^2 n^{-s} from w_n <= (log n)^2/2 + 1,
    with a factor-2 safety margin and first-term slack."""
    L = math.log(N)
    first = ((L * L / 2 + 1) ** 2) * N ** (-s)
    return 2.0 * (integral_tail(4, s, N) / 4 + integral_tail(2, s, N)
                  + integral_tail(0, s, N) + first)


def tail_wlog(s: float, N: int) -> float:
    L = math.log(N)
    first = (L ** 3 / 2 + L) * N ** (-s)
    return 2.0 * (integral_tail(3, s, N) / 2 + integral_tail(1, s, N) + first)


# ── the gates ───────────────────────────────────────────────────────

def main() -> int:
    ledger = {"study": "w2_susceptibility", "dps": mp.mp.dps, "gates": {}}
    ok = True

    def gate(name, passed, **info):
        nonlocal ok
        ok &= bool(passed)
        ledger["gates"][name] = {"passed": bool(passed), **info}
        print(f"[{'PASS' if passed else 'FAIL'}] {name}  " +
              "  ".join(f"{k}={v}" for k, v in info.items()))

    # G1 algebra gate
    for s, D, tol in ((3.0, 3000, 1e-8), (2.5, 3000, 1e-6)):
        tb, tc = brute_T(s, D), float(weight_second_moment_exact(s))
        gate(f"G1_algebra_s{s}", abs(tb - tc) < tol,
             gap=f"{abs(tb-tc):.2e}", tol=tol)

    # G2 D2 identity, sieve vs closed, tail-majorized
    for s, N in ((2.5, 400_000), (3.0, 200_000)):
        w = sieve_w(N)
        n = np.arange(1, N + 1, dtype=float)
        lhs = float((w[1:] ** 2) @ n ** (-s))
        rhs = float(_z(s) * weight_second_moment_exact(s))
        t = tail_w2(s, N)
        gate(f"G2_D2_identity_s{s}", abs(lhs - rhs) <= t,
             gap=f"{abs(lhs-rhs):.3e}", tail=f"{t:.3e}")

    # G3 Cov series route
    s, N = 2.5, 400_000
    w = sieve_w(N)
    n = np.arange(1, N + 1, dtype=float)
    lhs = float((w[1:] * np.log(n)) @ n ** (-s))
    rhs = float(_z(s, 1) * _z(s + 1, 1) + _z(s) * _z(s + 1, 2))
    gate("G3_cov_series", abs(lhs - rhs) <= tail_wlog(s, N),
         gap=f"{abs(lhs-rhs):.3e}", tail=f"{tail_wlog(s, N):.3e}")

    # G4 Cov ensemble route
    beta = 2.5
    p = n ** (-beta); p /= p.sum()
    E = np.log(n)
    cov_ens = float(p @ (E * w[1:]) - (p @ E) * (p @ w[1:]))
    cov_ex = float(weight_energy_covariance_exact(beta))
    gate("G4_cov_ensemble", abs(cov_ens - cov_ex) < 1e-5,
         gap=f"{abs(cov_ens-cov_ex):.2e}")

    # G5 Gibbs consistency (route (a) of I2)
    d_mean = float(mp.diff(lambda b: -_z(b + 1, 1), mp.mpf(beta)))
    gate("G5_gibbs_consistency", abs(d_mean + cov_ex) < 1e-12,
         gap=f"{abs(d_mean+cov_ex):.2e}")

    # G6 Fisher scan: det > 0, rho in (0,1) increasing on grid
    grid = [1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    rhos, det_min = [], None
    for b in grid:
        g11, g12, g22 = weight_fisher_matrix(b)
        det = g11 * g22 - g12 ** 2
        rho = -g12 / mp.sqrt(g11 * g22)
        rhos.append(float(rho))
        det_min = det if det_min is None or det < det_min else det_min
    increasing = all(r2 > r1 for r1, r2 in zip(rhos, rhos[1:]))
    gate("G6_fisher_scan", det_min > 0 and increasing
         and all(0 < r < 1 for r in rhos),
         det_min=mp.nstr(det_min, 6), rho_first=f"{rhos[0]:.4f}",
         rho_last=f"{rhos[-1]:.4f}")
    ledger["rho_curve"] = dict(zip(map(str, grid), rhos))

    # G7 precision integrity of the cold constants
    cold30 = {
        "w2_cold": weight_second_moment_exact(1),
        "chi_cold": weight_susceptibility_exact(1),
        "cov_cold": weight_energy_covariance_exact(1),
    }
    mp.mp.dps = 45
    cold45 = {
        "w2_cold": weight_second_moment_exact(1),
        "chi_cold": weight_susceptibility_exact(1),
        "cov_cold": weight_energy_covariance_exact(1),
    }
    mp.mp.dps = 30
    drift = max(abs(cold30[k] - cold45[k]) for k in cold30)
    gate("G7_cold_precision", drift < mp.mpf(10) ** (-25),
         max_drift=mp.nstr(drift, 3))

    ledger["cold_constants_EXACT"] = {
        "w2_cold_T1": mp.nstr(cold45["w2_cold"], 30),
        "chi_w_cold": mp.nstr(cold45["chi_cold"], 30),
        "Dw_cold": mp.nstr(mp.sqrt(cold45["chi_cold"]), 30),
        "cov_cold_zetapp2": mp.nstr(cold45["cov_cold"], 30),
        "wall_rate_rho_coeff": mp.nstr(
            cold45["cov_cold"] / mp.sqrt(cold45["chi_cold"]), 30),
        "mean_cold_minus_zetap2": mp.nstr(-_z(2, 1), 30),
    }
    ledger["all_passed"] = ok

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "w2_susceptibility_ledger.json")
    with open(out, "w") as f:
        json.dump(ledger, f, indent=2)
    print(f"\nledger -> {out}")
    print("ALL GATES PASS" if ok else "GATE FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
