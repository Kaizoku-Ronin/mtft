# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
op3_spectral_study.py — The OP3 Spectral Study (third rung)
============================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

Question (note v0.1.1, OP3): is the spectral measure of the correlator
F_p(t) = omega(mu_p* alpha_t(mu_p)) sharp at the edge omega* = (beta+1)log p
— an isolated integrable singularity — or a smooth edge?

The measure is purely atomic:  sigma_p = sum_n rho_n delta(omega - dE_n),
dE_n = (beta+1)log p - g_n,  g_n = log(1 + log p/log n),  strictly
increasing to omega*.  The study therefore certifies four statements:

  Pr A (lines):   low-n atoms are isolated with O(1) weight; the n = 2
                  line alone carries rho_2 = log2 * 2^{-(b+1)}/(-zeta'(b+1))
                  (43.73% at beta = 2).  Detector: pi*eta*S_eta(dE_2) ->
                  rho_2 + O(eta^2).
  Pr B (edge):    sigma_p({omega*}) = 0 and the density in the gap
                  variable is C-infinity-flat:
                  f(g) = L^2 e^{-beta L} e^g / ((e^g - 1) Z2),
                  L = log p/(e^g - 1)  ~ (log p)^2/(Z2 g^3) e^{-beta log p/g}.
                  Detectors: pi*eta*S_eta(omega*) = C* eta^2 + o(eta^2) with
                  C* = sum rho_n/g_n^2 computed independently (the rate is
                  Lorentzian-tail dominated — resolution-limited spectroscopy
                  CANNOT see the essential singularity; the windowed edge law,
                  gate SG0, is the true edge probe); plus binned line density
                  vs the exact f(g).
  Pr C (count):   shell counts #{eps1 <= g_n < eps2} = floor(M(eps1)) -
                  floor(M(eps2)), M(eps) = p^{1/(e^eps - 1)} — lines
                  proliferate super-exponentially while shell weight dies
                  super-exponentially.
  Pr D (Wiener / no clustering):  all dE_n distinct  =>
                  lim T^-1 int_0^T |F(t)|^2 dt = sum_n rho_n^2
                  = zeta''(2b+2)/zeta'(b+1)^2  > 0.
                  Correlations NEVER decay: the marked-gas algebra is not
                  mixing, OS4-style clustering is impossible without spatial
                  extension.  This is the reconstruction-theory payload.

Gates: SG0 admission (edge law, module instrument), SG1 atom detector,
SG2 edge emptiness + C* two-leg, SG3 density law vs binned data,
SG4 Wiener three legs (exact finite-T pair sum / time sampling / closed
form), SG5 shell count law, SG6 recurrence (DIAGNOSTIC).

Instrument: mtft >= 0.9.1 (edge_mass with predicted_em / mass_plus_tail).
Run:  py op3_spectral_study.py
"""

from __future__ import annotations
import math
import numpy as np
import mpmath as mp
import mtft.marked_gas as mg

BETA, P = 2.0, 2
REPORT = []


def rec(name, gtype, value, cls, ok, note=""):
    REPORT.append((name, gtype, value, cls, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name:<30} {gtype:<12} "
          f"{value:<30} {cls:<18} {note}")


def levels(nmax, beta=BETA, p=P):
    n = np.arange(2, nmax + 1, dtype=np.float64)
    with mp.workdps(30):
        Z2 = float(-mp.zeta(beta + 1, derivative=1))
    rho = np.log(n) * n ** (-(beta + 1.0)) / Z2
    g = np.log1p(math.log(p) / np.log(n))
    return n, rho, g, Z2


# ----------------------------------------------------------------------
# SG0 — admission control: the edge law on the shipped instrument
# ----------------------------------------------------------------------
def sg0():
    worst = 0.0
    for eps in (0.3, 0.2, 0.1, 0.07):
        d = mg.edge_mass(P, BETA, eps, 2_000_000).detail
        worst = max(worst, abs(d["ratio_em_corrected"] - 1.0))
    rec("SG0 edge-law admission", "Instrument", f"max|r_em_corr-1|={worst:.1e}",
        "CERTIFIED(1e-6)", worst < 1e-6, "mtft 0.9.1 edge_mass, 4 eps values")


# ----------------------------------------------------------------------
# SG1 — atom detector at the dominant line
# ----------------------------------------------------------------------
def sg1(nmax=200_000):
    n, rho, g, Z2 = levels(nmax)
    dE = (BETA + 1.0) * math.log(P) - g
    with mp.workdps(30):
        rho2_exact = float(mp.log(2) * mp.mpf(2) ** (-(BETA + 1))
                           / (-mp.zeta(BETA + 1, derivative=1)))
    target = dE[0]                        # the n = 2 line
    vals = {}
    for eta in (1e-2, 1e-3, 1e-4):
        S = np.sum(rho * (eta / math.pi) / ((target - dE) ** 2 + eta ** 2))
        vals[eta] = math.pi * eta * S
    # quadratic approach: (val - rho2)/eta^2 stable across eta
    c1 = (vals[1e-2] - rho2_exact) / 1e-4
    c2 = (vals[1e-3] - rho2_exact) / 1e-6
    ok = (abs(vals[1e-4] - rho2_exact) < 1e-7
          and abs(c1 / c2 - 1) < 0.05)
    rec("SG1 atom at dE_2", "Spectral", f"pi*eta*S={vals[1e-4]:.9f}",
        "CERTIFIED(1e-7)", ok,
        f"rho_2={rho2_exact:.9f} (43.73%), O(eta^2) coeff stable "
        f"{c1:.4f}/{c2:.4f}")


# ----------------------------------------------------------------------
# SG2 — edge emptiness: pi*eta*S_eta(omega*) = C* eta^2, C* independent
# ----------------------------------------------------------------------
def sg2(nmax=2_000_000):
    n, rho, g, Z2 = levels(nmax)
    Cstar = float(np.sum(rho / g ** 2))   # converges; tail beyond nmax tiny
    ok_all, worst = True, 0.0
    for eta in (1e-2, 1e-3):
        S = np.sum(rho * (eta / math.pi) / (g ** 2 + eta ** 2))
        ratio = (math.pi * eta * S) / (Cstar * eta ** 2)
        worst = max(worst, abs(ratio - 1.0))
        ok_all &= abs(ratio - 1.0) < 2e-2
    rec("SG2 no atom at edge", "Spectral", f"C*={Cstar:.6f}",
        "CERTIFIED(2e-2)", ok_all,
        f"pi*eta*S(omega*)/(C* eta^2) within {worst:.1e} of 1; "
        "rate is Lorentzian-tail (essential singularity invisible to eta)")


# ----------------------------------------------------------------------
# SG3 — the density law f(g) against binned line data
# ----------------------------------------------------------------------
def sg3(nmax=2_000_000):
    n, rho, g, Z2 = levels(nmax)

    def f_exact(gv):
        L = math.log(P) / (math.exp(gv) - 1.0)
        return (L ** 2 * math.exp(-BETA * L) * math.exp(gv)
                / ((math.exp(gv) - 1.0) * Z2))

    edges = [0.05, 0.08, 0.12, 0.18, 0.25]
    worst = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        direct = float(rho[(g >= lo) & (g < hi)].sum())
        with mp.workdps(25):
            pred = float(mp.quad(lambda gv: f_exact(float(gv)), [lo, hi]))
        worst = max(worst, abs(direct / pred - 1.0))
    rec("SG3 density law f(g)", "Structural", f"max bin dev={worst:.2e}",
        "CERTIFIED(2e-2)", worst < 2e-2,
        f"4 bins in [0.05,0.25], change-of-variables exact form")


# ----------------------------------------------------------------------
# SG4 — Wiener / no clustering: three legs
# ----------------------------------------------------------------------
def sg4(nbasis=4000):
    n, rho, g, Z2 = levels(nbasis)
    dE = (BETA + 1.0) * math.log(P) - g
    with mp.workdps(30):
        closed = float(mp.zeta(2 * BETA + 2, derivative=2)
                       / mp.zeta(BETA + 1, derivative=1) ** 2)
    direct = float(np.sum(rho ** 2))      # leg 2: frequency-domain sum
    # leg 1: exact finite-T pair sum  A(T) = sum rho_m rho_n phi((dEm-dEn)T)
    diffs = np.subtract.outer(dE, dE)
    rr = np.outer(rho, rho)
    A = {}
    for T in (50.0, 400.0, 3200.0):
        x = diffs * T
        phi = np.where(np.abs(x) < 1e-12, 1.0,
                       (np.exp(1j * x) - 1.0) / (1j * np.where(x == 0, 1, x)))
        A[T] = float(np.real(np.sum(rr * phi)))
    # leg 3: time sampling Cesaro mean of |F|^2
    t = np.arange(0.0, 3200.0, 0.05)
    F = np.zeros(len(t), dtype=complex)
    for chunk in range(0, len(dE), 500):
        F += np.exp(1j * np.outer(t, dE[chunk:chunk + 500])) \
             @ rho[chunk:chunk + 500].astype(complex)
    leg3 = float(np.mean(np.abs(F) ** 2))
    ok = (abs(direct - closed) < 5e-6
          and abs(A[3200.0] - closed) < 2e-3
          and abs(A[3200.0] - A[400.0]) < abs(A[400.0] - A[50.0])
          and abs(leg3 - closed) < 2e-3)
    rec("SG4 Wiener (no clustering)", "Theorem", f"lim={closed:.9f}",
        "CERTIFIED(2e-3)", ok,
        f"legs: closed {closed:.6f} | sum {direct:.6f} | "
        f"A(3200)={A[3200.0]:.6f} | time-avg {leg3:.6f}  — |F|^2 never decays")


# ----------------------------------------------------------------------
# SG5 — shell count law
# ----------------------------------------------------------------------
def sg5(nmax=2_000_000):
    n, rho, g, Z2 = levels(nmax)
    ok = True
    dets = []
    for e1, e2 in ((0.1, 0.2), (0.07, 0.1), (0.06, 0.07)):
        M1 = math.exp(math.log(P) / (math.exp(e1) - 1.0))
        M2 = math.exp(math.log(P) / (math.exp(e2) - 1.0))
        pred = math.floor(M1) - math.floor(M2)
        count = int(np.sum((g >= e1) & (g < e2)))
        ok &= (count == pred)
        dets.append(f"[{e1},{e2}):{count}={pred}")
    rec("SG5 shell count law", "Identity", "; ".join(dets), "EXACT", ok,
        "floor(M(e1))-floor(M(e2)); proliferation vs weight starvation")


# ----------------------------------------------------------------------
# SG6 — recurrence (DIAGNOSTIC): almost-periodic, |F| revisits 1
# ----------------------------------------------------------------------
def sg6(nbasis=2000):
    n, rho, g, Z2 = levels(nbasis)
    dE = (BETA + 1.0) * math.log(P) - g
    t = np.arange(5.0, 4000.0, 0.05)
    F = np.zeros(len(t), dtype=complex)
    for chunk in range(0, len(dE), 500):
        F += np.exp(1j * np.outer(t, dE[chunk:chunk + 500])) \
             @ rho[chunk:chunk + 500].astype(complex)
    m = float(np.max(np.abs(F)))
    tm = float(t[int(np.argmax(np.abs(F)))])
    rec("SG6 recurrence", "Diagnostic", f"max|F| on [5,4000] = {m:.4f}",
        "DIAGNOSTIC", True, f"at t={tm:.2f}; |F(0)|=1, no decay ever")


if __name__ == "__main__":
    print("=" * 100)
    print("  OP3 SPECTRAL STUDY — marked primon gas, p = 2, beta = 2  "
          "(instrument: mtft 0.9.1)")
    print("=" * 100)
    sg0(); sg1(); sg2(); sg3(); sg4(); sg5(); sg6()
    print("-" * 100)
    n_pass = sum(1 for r in REPORT if r[4])
    print(f"  {n_pass}/{len(REPORT)} gates green")
    print("  VERDICT: pure point spectrum; isolated O(1) lines (Rydberg-"
          "like series); edge C^inf-flat, no atom;")
    print("  correlations never decay (Wiener limit > 0) => clustering "
          "impossible without spatial extension (OP5 forced).")
    print("=" * 100)
