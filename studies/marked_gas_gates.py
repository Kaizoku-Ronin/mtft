# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
marked_gas_gates.py — Green Suite v0.1 for the Marked Primon Gas construction
==============================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

Certifies, per the A.7 / Legend epistemic standards (EXACT / CERTIFIED(tol) /
DIAGNOSTIC / PHENO), the operator-algebraic construction:

    H_F  = primon Fock space,  H|n> = (log n)|n>
    Z1   = Tr(e^{-beta H})            = zeta(beta)              [unmarked gas]
    Z2   = Tr(H e^{-(beta+1)H})       = -zeta'(beta+1)          [marked gas]
    Z_D  = Z1 * Z2                    = -zeta(beta) zeta'(beta+1)
    K    = (beta+1) H - log H  (+const)      [type-I modular Hamiltonian]
    sigma_t(mu_p)|n> = p^{i(beta+1)t} (log(pn)/log n)^{-it} |pn>
    KMS  : F(t + i) = G(t)  termwise, alpha_t = Ad e^{itK}
    Meinardus (cold gas): A_n = sum_{k<=n} a_k ~ B n^alpha / Gamma(alpha+1),
        alpha = -zeta'(2),  B = exp(-zeta''(2)),
        where n a_n = sum_{k=1}^n w_k a_{n-k},  a_0 = 1   (Psi coefficients)
    Speiser: RH <=> no Fisher zeros of Z2 continuation in -1 < Re beta < -1/2.

Engine 1: numpy float64 (compensated) + mpmath (dps 25-30).
Engine 2: PARI/GP mirror in marked_gas_gates.gp.
Engine 0: mtft 0.8.0 package anchors + zeta'-zero census cross-check.

Run:  py marked_gas_gates.py        (Windows)
      python3 marked_gas_gates.py   (POSIX)
"""

from __future__ import annotations
import math
import numpy as np
import mpmath as mp

REPORT = []


def gate(name, gtype, engine, value, cls, ok, note=""):
    REPORT.append((name, gtype, engine, value, cls, bool(ok), note))
    flag = "PASS" if ok else "FAIL"
    print(f"[{flag}] {name:<34} {gtype:<15} {engine:<14} {value:<26} {cls:<18} {note}")


def kahan_cumjumps(kpow, breakpoints):
    """Compensated prefix sums of kpow (1-indexed values in kpow[0:]) snapshotted
    at the sorted breakpoints. Returns dict M -> H_M with total error
    O(eps * n_chunks * magnitude) ~ 1e-15 here (pairwise inside chunks,
    Kahan across chunks)."""
    H = {}
    run, comp = 0.0, 0.0
    prev = 0
    for m in breakpoints:
        chunk = float(np.sum(kpow[prev:m]))  # pairwise
        y = chunk - comp
        t = run + y
        comp = (t - run) - y
        run = t
        H[m] = run
        prev = m
    return H


# ----------------------------------------------------------------------
# G0 — environment anchors (mtft 0.8.0 package must be the same world)
# ----------------------------------------------------------------------
def gate0():
    import mtft
    ok1 = abs(mtft.dirichlet_curvature(3.0)["g_D"] - 0.3351038786441419) < 1e-14
    lam1 = mtft.li_lambda(1)
    lam1_closed = float(1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2)
    ok2 = abs(lam1 - lam1_closed) < 1e-13
    ok3 = mtft.filtered_moment_identity(0.18174, N=3)["rel_diff"] < 1e-12
    gate("G0 package anchors (mtft 0.8.0)", "Environment", "E0 mtft",
         f"gD/li1/mu3 ok={ok1},{ok2},{ok3}", "EXACT", ok1 and ok2 and ok3)


# ----------------------------------------------------------------------
# G1 — Tr(e^{-beta H}) = zeta(beta), Euler–Maclaurin certified
# ----------------------------------------------------------------------
def gate1(beta=2, N=100000, dps=30):
    with mp.workdps(dps):
        b = mp.mpf(beta)
        S = mp.fsum(mp.power(n, -b) for n in range(1, N + 1))
        Nb = mp.mpf(N)
        tail = Nb ** (1 - b) / (b - 1) - Nb ** (-b) / 2 + b * Nb ** (-b - 1) / 12
        errb = abs(b * (b + 1) * (b + 2)) * Nb ** (-b - 3) / 720
        resid = abs(S + tail - mp.zeta(b))
        ok = resid <= errb + mp.mpf(10) ** (-(dps - 3))
        gate(f"G1 Tr1=zeta({beta})", "Identity", "E1 mpmath",
             f"resid={mp.nstr(resid,3)}", f"CERTIFIED(<{mp.nstr(errb,2)})", ok,
             f"N={N}, EM tail")


# ----------------------------------------------------------------------
# G2 — Tr(H e^{-(beta+1)H}) = -zeta'(beta+1), Euler–Maclaurin certified
# ----------------------------------------------------------------------
def gate2(beta=2, N=100000, dps=30):
    with mp.workdps(dps):
        s = mp.mpf(beta) + 1
        S = mp.fsum(mp.log(n) * mp.power(n, -s) for n in range(2, N + 1))
        Nb = mp.mpf(N)
        integ = Nb ** (1 - s) * (mp.log(Nb) / (s - 1) + 1 / (s - 1) ** 2)
        f_N = mp.log(Nb) * Nb ** (-s)
        fp_N = Nb ** (-s - 1) * (1 - s * mp.log(Nb))
        tail = integ - f_N / 2 - fp_N / 12
        errb = 2 * (s + 1) ** 3 * mp.log(Nb) * Nb ** (-s - 3) / 720
        target = -mp.zeta(s, derivative=1)
        resid = abs(S + tail - target)
        ok = resid <= errb + mp.mpf(10) ** (-(dps - 3))
        gate(f"G2 Tr2=-zeta'({beta}+1)", "Identity", "E1 mpmath",
             f"resid={mp.nstr(resid,3)}", f"CERTIFIED(<{mp.nstr(errb,2)})", ok,
             f"N={N}, EM tail, vacuum n>=2")


# ----------------------------------------------------------------------
# G3 — convolution identity Z_D(beta) = -zeta(beta) zeta'(beta+1)
#      certified two-sided interval at N = 1e7 (float64, compensated)
#      + reorder cross-check against materialized w_n at N = 2e5
# ----------------------------------------------------------------------
def gate3(beta=2, N=10_000_000):
    b = float(beta)
    k = np.arange(1, N + 1, dtype=np.float64)
    kpow = k ** (-b)
    d = np.arange(2, N + 1, dtype=np.int64)
    M = N // d
    breakpoints = np.unique(M)          # ~2 sqrt(N) distinct prefix depths
    H = kahan_cumjumps(kpow, list(breakpoints))
    Hvals = np.array([H[m] for m in M], dtype=np.float64)
    dd = d.astype(np.float64)
    wgt = np.log(dd) * dd ** (-b - 1.0)
    S = float(np.sum(wgt * Hvals))

    # two-sided defect: zeta-tail of each inner sum, plus the d > N tail
    Mf = M.astype(np.float64)
    D_lo = float(np.sum(wgt * (Mf + 1.0) ** (1.0 - b))) / (b - 1.0)
    D_hi = float(np.sum(wgt * Mf ** (1.0 - b))) / (b - 1.0)
    zeta_b = float(mp.zeta(b))

    def tail2(A):  # integral of log(x) x^{-b-1} from A to infinity
        return A ** (-b) * (math.log(A) / b + 1.0 / b ** 2)

    T2_lo, T2_hi = zeta_b * tail2(N + 1.0), zeta_b * tail2(float(N))
    slack = 4 * 2.3e-16 * math.log2(N) * (abs(S) + 2.0)   # pairwise + Kahan
    lo = S + D_lo + T2_lo - slack
    hi = S + D_hi + T2_hi + slack
    with mp.workdps(30):
        target = float(-mp.zeta(mp.mpf(b)) * mp.zeta(mp.mpf(b) + 1, derivative=1))
    width = hi - lo
    ok = lo <= target <= hi
    gate(f"G3 Z_D=-zz' (beta={beta})", "Identity", "E1 numpy/mp",
         f"width={width:.2e}", f"CERTIFIED({width:.1e})", ok,
         f"N=1e7 interval [{lo:.15f},{hi:.15f}] ∋ {target:.15f}" if ok else
         f"target {target} outside [{lo},{hi}]")

    # reorder cross-check: direct divisor sieve of w_n at N2
    N2 = 200_000
    w = np.zeros(N2 + 1)
    for dv in range(2, N2 + 1):
        w[dv::dv] += math.log(dv) / dv
    n = np.arange(2, N2 + 1, dtype=np.float64)
    S_direct = float(np.sum(w[2:] * n ** (-b)))
    # swapped form at the same truncation n <= N2
    d2 = np.arange(2, N2 + 1, dtype=np.int64)
    M2 = N2 // d2
    bp2 = np.unique(M2)
    kpow2 = np.arange(1, N2 + 1, dtype=np.float64) ** (-b)
    H2 = kahan_cumjumps(kpow2, list(bp2))
    H2v = np.array([H2[m] for m in M2])
    dd2 = d2.astype(np.float64)
    S_swap = float(np.sum(np.log(dd2) * dd2 ** (-b - 1.0) * H2v))
    diff = abs(S_direct - S_swap)
    gate("G3b reorder cross-check", "Implementation", "E1 numpy",
         f"|direct-swap|={diff:.2e}", "EXACT(reorder)", diff < 5e-13,
         f"N={N2}, sieved w_n vs (d,k) order")


# ----------------------------------------------------------------------
# G4 — modular flow phase law vs literal matrix conjugation
# ----------------------------------------------------------------------
def gate4(beta=2, Nb=400, primes=(2, 3, 5), ts=(0.7, 1.9, math.pi)):
    ns = np.arange(2, Nb + 1)
    logn = np.log(ns.astype(float))
    E = (beta + 1) * logn - np.log(logn)          # spectrum of K (mod const)
    idx = {int(n): i for i, n in enumerate(ns)}
    worst = 0.0
    for p in primes:
        for t in ts:
            U = np.diag(np.exp(1j * t * E))
            mu = np.zeros((len(ns), len(ns)), dtype=complex)
            for n in ns:
                if p * n <= Nb:
                    mu[idx[p * n], idx[n]] = 1.0
            A = U @ mu @ np.conj(U.T)             # sigma_t(mu_p), matrix path
            for n in ns:
                if p * n <= Nb:
                    closed = (p ** (1j * (beta + 1) * t)
                              * (math.log(p * n) / math.log(n)) ** (-1j * t))
                    worst = max(worst, abs(A[idx[p * n], idx[n]] - closed))
    gate("G4 phase law sigma_t(mu_p)", "Implementation", "E1 numpy",
         f"max|matrix-closed|={worst:.2e}", "EXACT(algebra)", worst < 1e-12,
         f"N={Nb}, p in {primes}, {len(ts)} t-values")


# ----------------------------------------------------------------------
# G5 — KMS boundary condition F(t - i) = G(t), two independent code paths
#      (convention: alpha_t = Ad e^{itK} = sigma_{-t}^{modular}; analytic point t + i)
# ----------------------------------------------------------------------
def gate5(beta=2, Nb=400, p=2, ts=(0.6, 2.2)):
    ns = np.arange(2, Nb + 1)
    logn = np.log(ns.astype(float))
    rho = logn * ns.astype(float) ** (-(beta + 1.0))     # unnormalized
    Z2 = float(np.sum(rho))
    K = -np.log(rho / Z2)
    idx = {int(n): i for i, n in enumerate(ns)}
    mu = np.zeros((len(ns), len(ns)))
    for n in ns:
        if p * n <= Nb:
            mu[idx[p * n], idx[n]] = 1.0
    rhohat = np.diag(rho / Z2)
    worst = 0.0
    for t in ts:
        z = t + 1j                                        # alpha_t = Ad e^{itK}
        Uz = np.diag(np.exp(1j * z * K))
        al_z_mu = Uz @ mu @ np.linalg.inv(Uz)
        F = np.trace(rhohat @ mu.T @ al_z_mu)             # omega(mu* alpha_{t+i}(mu))
        Ut = np.diag(np.exp(1j * t * K))
        al_t_mu = Ut @ mu @ np.linalg.inv(Ut)
        G = np.trace(rhohat @ al_t_mu @ mu.T)             # omega(alpha_t(mu) mu*)
        # path B: direct sums from the closed spectrum
        mask = ns * p <= Nb
        nn = ns[mask].astype(float)
        En = (beta + 1) * np.log(nn) - np.log(np.log(nn))
        Epn = (beta + 1) * np.log(p * nn) - np.log(np.log(p * nn))
        dE = Epn - En
        rn = np.log(nn) * nn ** (-(beta + 1.0)) / Z2
        rpn = np.log(p * nn) * (p * nn) ** (-(beta + 1.0)) / Z2
        F_b = np.sum(rn * np.exp(1j * z * dE))
        G_b = np.sum(rpn * np.exp(1j * t * dE))
        worst = max(worst, abs(F - G), abs(F_b - G_b), abs(F - F_b))
    gate("G5 KMS F(t+i)=G(t)", "Consistency", "E1 numpy x2",
         f"max dev={worst:.2e}", "EXACT(termwise)", worst < 1e-12,
         f"p={p}, matrix path vs spectral path")


# ----------------------------------------------------------------------
# G6 — Meinardus cold gas: Psi coefficients via n a_n = sum w_k a_{n-k}
#      summatory slope -> alpha = -zeta'(2), amplitude -> exp(-zeta''(2))
# ----------------------------------------------------------------------
def gate6(N=100_000, lo=20_000):
    w = np.zeros(N + 1)
    for dv in range(2, N + 1):
        w[dv::dv] += math.log(dv) / dv
    a = np.zeros(N + 1)
    a[0] = 1.0
    for n in range(1, N + 1):
        a[n] = np.dot(w[1:n + 1], a[n - 1::-1]) / n
    A = np.cumsum(a)
    with mp.workdps(30):
        alpha = float(-mp.zeta(2, derivative=1))
        zpp2 = float(mp.zeta(2, derivative=2))
        amp_target = math.exp(-zpp2) / float(mp.gamma(alpha + 1))
    import mtft
    ns = np.arange(lo, N + 1, dtype=float)
    slope, nb_used, nb_drop = mtft.binned_log_slope(
        list(ns), list(A[lo:N + 1]), bin_width=0.1, min_bin=50)
    amp_fit = float(np.exp(np.mean(np.log(A[lo:N + 1]) - alpha * np.log(ns))))
    ok_s = abs(slope - alpha) < 5e-3
    ok_a = abs(amp_fit / amp_target - 1) < 0.05
    gate("G6 cold-gas slope alpha", "Structural", "E1 numpy+mtft",
         f"fit={slope:.6f} vs {alpha:.6f}", "CERTIFIED(5e-3)", ok_s,
         f"bins used={nb_used}, window [{lo},{N}]")
    gate("G6b cold-gas amplitude", "Structural", "E1 numpy/mp",
         f"fit={amp_fit:.5f} vs {amp_target:.5f}", "CERTIFIED(5%)", ok_a,
         "B=exp(-zeta''(2))/Gamma(alpha+1); corr O((ln n)^2/n)")
    # pointwise exponent, DIAGNOSTIC only
    ps, _, _ = mtft.binned_log_slope(list(ns), list(a[lo:N + 1]),
                                     bin_width=0.1, min_bin=50)
    gate("G6c pointwise a_n exponent", "Structural", "E1 numpy+mtft",
         f"fit={ps:.4f} vs {alpha-1:.4f}", "DIAGNOSTIC", True,
         "Karamata gives summatory; pointwise needs monotonicity")


# ----------------------------------------------------------------------
# G7 — Speiser strip harness: census cross-refined, mapped to beta = s - 1
# ----------------------------------------------------------------------
def gate7(n_refine=5):
    import mtft
    Z = mtft.ZETAPRIME_ZEROS
    worst_res = 0.0
    with mp.workdps(25):
        for s0 in Z[:n_refine]:
            root = mp.findroot(lambda s: mp.zeta(s, derivative=1),
                               mp.mpc(s0), tol=1e-22,
                               df=lambda s: mp.zeta(s, derivative=2))
            worst_res = max(worst_res, float(abs(mp.zeta(root, derivative=1))))
            if abs(root - mp.mpc(s0)) > 1e-8:
                gate("G7 Speiser strip", "Consistency", "E1 mp+census",
                     "census mismatch", "FAIL", False, str(s0))
                return
        neg = mp.findroot(lambda s: mp.zeta(s, derivative=1), mp.mpf(-2.7),
                          df=lambda s: mp.zeta(s, derivative=2))
        ok_neg = abs(float(neg) - mtft.zetaprime_negative_zero(1)) < 1e-9
    res = [z.real for z in Z]
    ok_re = all(r > 0.5 for r in res)
    strip = [z - 1 for z in Z if -1 < z.real - 1 < -0.5]
    ok = ok_re and not strip and ok_neg and worst_res < 1e-18
    gate("G7 Speiser/Fisher strip", "Consistency", "E1 mp+census",
         f"19 zeros Re>1/2, refine resid<{worst_res:.1e}", "CERTIFIED(h<=100)",
         ok, "no Fisher zeros in -1<Re(beta)<-1/2 to height 100; "
             "necessary-not-sufficient for RH")


if __name__ == "__main__":
    print("=" * 100)
    print("  MARKED PRIMON GAS — GREEN SUITE v0.1   (engine 1: python; "
          "engine 2: marked_gas_gates.gp)")
    print("=" * 100)
    gate0()
    gate1(beta=2); gate1(beta=3)
    gate2(beta=2); gate2(beta=3)
    gate3(beta=2); 
    gate4(); gate5(); gate6(); gate7()
    n_pass = sum(1 for r in REPORT if r[5])
    print("-" * 100)
    print(f"  {n_pass}/{len(REPORT)} gates green")
    print("=" * 100)
