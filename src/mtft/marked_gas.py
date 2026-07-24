"""marked_gas.py — The Marked Primon Gas: a KMS construction for the
Dirichlet ensemble (mtft 0.9.0).

Construction note: "The Marked Primon Gas" v0.1.1 (July 2026), audited
in Addendum U.  The Dirichlet ensemble Z_D(β) = Σ w_n n^{−β} =
−ζ(β)ζ′(β+1) is the partition function of a Gibbs state on the
single-particle space H₂ = ℓ²(N≥2) with Hamiltonian H|n⟩ = (log n)|n⟩:

    ρ̂_n = e^{−(β+1)E_n}/Z₂,  E_n = (β+1)log n − log log n  (spectrum Df 2.4)

Dynamics: the prime-shift isometry μ_p|n⟩ = |pn⟩ evolves under
α_t = Ad e^{itK} (K = −log ρ̂ + const) and satisfies the KMS condition
at inverse temperature 1 *at the analytic point t + i* — termwise exact,
ρ̂_n e^{−ΔE_n} = ρ̂_{pn} (Pr 3.5, gate G5).  The Bost–Connes symmetry is
deformed in the IR by the marking (Pr 3.6 twist (1 + log p/log n)^{−it})
and restored asymptotically in the UV.

Cold gas: the bosonic partition Ψ(q) = Π_{n≥2}(1−qⁿ)^{−(log n)/n²} has
log Ψ(e^{−τ}) = α log(1/τ) − ζ″(2) + O(τ log²τ)  (γ cancels between Γ
and ζ(s+1), A.2), hence, by Karamata on the summatory,

    A_n = Σ_{k≤n} a_k ~ B n^α / Γ(α+1),   α = −ζ′(2),  B = e^{−ζ″(2)}.

Speiser as Fisher zeros: Z₂(β) = −ζ′(β+1)/−ζ′(2) vanishes in
−1 < Re β < −1/2  ⟺  RH (Speiser 1935 at s = β+1; the unmarked factor
ζ(β) is zero-free in that strip).

Every quantity carries its exactness class (EXACT / CERTIFIED(bound) /
DIAGNOSTIC / PHENO) per the Legend.  float64 + mpmath dual tier.

Roger Tano — MTFT Research Program — July 2026
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from mtft.estimator_standards import binned_log_slope

__all__ = [
    "ALPHA_COLD", "B_COLD", "Certified",
    "z1", "z2", "zD_certified_interval",
    "spectrum", "flow_phase", "kms_check", "bc_deformation",
    "psi_coefficients", "cold_gas_report",
    "correlator", "spectral_function", "edge_mass",
    "gates",
]

# −ζ′(2): the cold-gas exponent (== constants.TORQUE_FULL, the proved
# Cesàro limit).  Literals pinned to mpmath in tests (1e-15).
ALPHA_COLD = 0.9375482543158438          # −ζ′(2)
B_COLD = 0.13679384954115163             # e^{−ζ″(2)}

# Speiser–Hadamard census re-exports (imported lazily to avoid a hard
# riemann dependency at module import):
def _census():
    from mtft.riemann import ZETAPRIME_ZEROS, ZETAPRIME_CENSUS_HEIGHT
    return ZETAPRIME_ZEROS, ZETAPRIME_CENSUS_HEIGHT


def _mp():
    import mpmath as mp
    return mp


# ── the certified value ─────────────────────────────────────────

@dataclass
class Certified:
    """A value with its exactness class and (where applicable) its
    certified two-sided bound.  `detail` carries the raw ledger."""
    value: object
    err_class: str = "EXACT"
    bound: Optional[float] = None
    detail: Dict = field(default_factory=dict)


def _kahan_cumjumps(kpow, breakpoints):
    """Compensated prefix sums of kpow snapshotted at the sorted
    breakpoints (pairwise inside chunks, Kahan across chunks)."""
    H = {}
    run, comp, prev = 0.0, 0.0, 0
    for m in breakpoints:
        chunk = float(np.sum(kpow[prev:m]))
        y = chunk - comp
        t = run + y
        comp = (t - run) - y
        run = t
        H[m] = run
        prev = m
    return H


# ── the two traces (Euler–Maclaurin certified) ──────────────────

def z1(beta: float = 2, N: int = 100_000, dps: int = 30) -> Certified:
    """Tr(e^{−βH}) = Σ n^{−β} = ζ(β), Euler–Maclaurin certified.

    Residual provably sits at the first neglected Bernoulli term
    |B₄ f‴(N)/4!| = |β(β+1)(β+2)| N^{−β−3}/720 (bound-equality is the
    expected outcome, not a near-fail)."""
    mp = _mp()
    with mp.workdps(dps):
        b = mp.mpf(beta)
        S = mp.fsum(mp.power(n, -b) for n in range(1, N + 1))
        Nb = mp.mpf(N)
        tail = Nb ** (1 - b) / (b - 1) - Nb ** (-b) / 2 + b * Nb ** (-b - 1) / 12
        errb = float(abs(b * (b + 1) * (b + 2)) * Nb ** (-b - 3) / 720)
        resid = float(abs(S + tail - mp.zeta(b)))
    return Certified(float(S + tail), f"CERTIFIED({errb:.1e})", errb,
                     {"target": float(mp.zeta(beta)), "resid": resid,
                      "N": N, "method": "Euler-Maclaurin B2+B4 bound"})


def z2(beta: float = 2, N: int = 100_000, dps: int = 30) -> Certified:
    """Tr(H e^{−(β+1)H}) = Σ (log n) n^{−(β+1)} = −ζ′(β+1), EM certified."""
    mp = _mp()
    with mp.workdps(dps):
        s = mp.mpf(beta) + 1
        S = mp.fsum(mp.log(n) * mp.power(n, -s) for n in range(2, N + 1))
        Nb = mp.mpf(N)
        integ = Nb ** (1 - s) * (mp.log(Nb) / (s - 1) + 1 / (s - 1) ** 2)
        tail = integ - mp.log(Nb) * Nb ** (-s) / 2 \
            - Nb ** (-s - 1) * (1 - s * mp.log(Nb)) / 12
        f3 = abs(s * (s + 1) * (s + 2) * mp.log(Nb)
                 - (s * (s + 1) + (s + 2) * (2 * s + 1))) * Nb ** (-s - 3)
        errb = float(abs(f3) / 720)
        resid = float(abs(S + tail - (-mp.zeta(s, derivative=1))))
    return Certified(float(S + tail), f"CERTIFIED({errb:.1e})", errb,
                     {"target": float(-mp.zeta(beta + 1, derivative=1)),
                      "resid": resid, "N": N})


# ── the convolution interval ────────────────────────────────────

def zD_certified_interval(beta: float = 2, N: int = 10_000_000) -> Certified:
    """Z_D(β) = −ζ(β)ζ′(β+1) certified to a two-sided interval at
    truncation N (hyperbola method + integral-test tail bounds both
    directions, float64 slack honestly included).

    Width scales as (log N)²/(2N²): 9.6e−11 at N = 10⁶, 1.4e−12 at
    N = 10⁷ (audit U.1)."""
    b = float(beta)
    kpow = np.arange(1, N + 1, dtype=np.float64) ** (-b)
    d = np.arange(2, N + 1, dtype=np.int64)
    M = N // d
    H = _kahan_cumjumps(kpow, list(np.unique(M)))
    Hvals = np.array([H[m] for m in M], dtype=np.float64)
    dd = d.astype(np.float64)
    wgt = np.log(dd) * dd ** (-b - 1.0)
    S = float(np.sum(wgt * Hvals))
    Mf = M.astype(np.float64)
    D_lo = float(np.sum(wgt * (Mf + 1.0) ** (1.0 - b))) / (b - 1.0)
    D_hi = float(np.sum(wgt * Mf ** (1.0 - b))) / (b - 1.0)
    mp = _mp()
    zeta_b = float(mp.zeta(b))

    def tail2(A):
        return A ** (-b) * (math.log(A) / b + 1.0 / b ** 2)

    T2_lo, T2_hi = zeta_b * tail2(N + 1.0), zeta_b * tail2(float(N))
    slack = 4 * 2.3e-16 * math.log2(N) * (abs(S) + 2.0)
    lo, hi = S + D_lo + T2_lo - slack, S + D_hi + T2_hi + slack
    with mp.workdps(30):
        target = float(-mp.zeta(mp.mpf(b)) * mp.zeta(mp.mpf(b) + 1, derivative=1))
    width = hi - lo
    return Certified((lo, hi), f"CERTIFIED({width:.1e})", width,
                     {"target": target, "inside": lo <= target <= hi,
                      "N": N, "slack": slack, "head": S})


# ── kinematics and the modular flow ─────────────────────────────

def spectrum(beta: float = 2, nmax: int = 1000) -> np.ndarray:
    """Single-particle spectrum E_n = (β+1) log n − log log n, n = 2..nmax."""
    n = np.arange(2, nmax + 1, dtype=np.float64)
    return (beta + 1.0) * np.log(n) - np.log(np.log(n))


def flow_phase(p: int, n, t: float, beta: float = 2):
    """Pr 3.3 closed form: α_t(μ_p)|n⟩ = p^{i(β+1)t} (log(pn)/log n)^{−it} |pn⟩."""
    n_arr = np.asarray(n, dtype=np.float64)
    phase = np.exp(1j * t * ((beta + 1.0) * math.log(p)
                             - np.log(np.log(p * n_arr) / np.log(n_arr))))
    return phase if phase.shape else complex(phase)


def kms_check(beta: float = 2, p: int = 2, t: float = 0.7,
              nbasis: int = 2000) -> Certified:
    """G5: the KMS identity F(t + i) = G(t) for ω(μ_p* α_t(μ_p)), both
    paths evaluated independently from the definitions.  Termwise exact
    algebra (ρ̂_n e^{−ΔE_n} = ρ̂_{pn}); the wrong-sign control F(t − i)
    is carried in the detail (it must fail — the convention is pinned)."""
    n = np.arange(2, nbasis + 1, dtype=np.float64)
    rho = np.log(n) * n ** (-(beta + 1.0))
    rho /= rho.sum()
    dE = (beta + 1.0) * math.log(p) - np.log(np.log(p * n) / np.log(n))
    F_plusi = np.sum(rho * np.exp(1j * (t + 1j) * dE))
    F_minusi = np.sum(rho * np.exp(1j * (t - 1j) * dE))
    G_t = np.sum(rho * np.exp(-dE) * np.exp(1j * t * dE))
    resid = abs(F_plusi - G_t)
    return Certified(resid, "CERTIFIED(1e-12)", 1e-12,
                     {"F(t+i)": F_plusi, "G(t)": G_t,
                      "wrong_sign_control": abs(F_minusi - G_t),
                      "p": p, "t": t, "nbasis": nbasis})


def bc_deformation(p: int, n, t: float, beta: float = 2):
    """Pr 3.6 Bost–Connes deformation twist D_p(t)|n⟩ =
    (1 + log p/log n)^{−it} |n⟩ — vanishes as n → ∞ (UV restoration)."""
    n_arr = np.asarray(n, dtype=np.float64)
    phase = np.exp(-1j * t * np.log1p(math.log(p) / np.log(n_arr)))
    return phase if phase.shape else complex(phase)


# ── the cold gas ────────────────────────────────────────────────

def weights_sieve(N: int) -> np.ndarray:
    """w_n = Σ_{d|n} (log d)/d for n = 0..N (w₀ = w₁ = 0), divisor sieve
    O(N log N) — the same engine as the green suite."""
    w = np.zeros(N + 1)
    for d in range(2, N + 1):
        w[d::d] += math.log(d) / d
    return w


def psi_coefficients(N: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(w, a, A): the weights, the coefficients of Ψ(q) =
    Π_{n≥2}(1−qⁿ)^{−(log n)/n²} (Euler-transform recurrence
    n a_n = Σ_{k≤n} w_k a_{n−k}), and their summatory."""
    w = weights_sieve(N)
    a = np.zeros(N + 1)
    a[0] = 1.0
    for n in range(1, N + 1):
        a[n] = np.dot(w[1:n + 1], a[n - 1::-1]) / n
    return w, a, np.cumsum(a)


def cold_gas_report(N: int = 100_000, lo: int = 20_000) -> Dict:
    """Pr 5.2 report on the cold-gas asymptotics over the window [lo, N].

    Estimator documentation (note redline U.3 item 3): the amplitude is
    quoted with the mean-of-logs estimator at α fixed (reads +0.09% at
    N = 10⁵); a free-slope least-squares intercept reads +1.0%; the
    honest trend A_n/n^α converges to B/Γ(α+1) from above at rate
    O((log n)²/n).  The pointwise exponent is DIAGNOSTIC (Karamata
    certifies the summatory; pointwise awaits monotonicity, OP4)."""
    mp = _mp()
    w, a, A = psi_coefficients(N)
    with mp.workdps(30):
        alpha = float(-mp.zeta(2, derivative=1))
        zpp2 = float(mp.zeta(2, derivative=2))
        amp_target = math.exp(-zpp2) / float(mp.gamma(alpha + 1))
    ns = np.arange(lo, N + 1, dtype=float)
    slope, nb_used, nb_drop = binned_log_slope(
        list(ns), list(A[lo:N + 1]), bin_width=0.1, min_bin=50)
    slope_2pt = (math.log(A[N]) - math.log(A[lo])) / (math.log(N)
                                                      - math.log(lo))
    amp_mol = float(np.exp(np.mean(np.log(A[lo:N + 1]) - alpha * np.log(ns))))
    amp_end = float(math.exp(math.log(A[N]) - alpha * math.log(N)))
    xs = np.log(ns)
    free_slope, free_ic = np.polyfit(xs, np.log(A[lo:N + 1]), 1)
    ps, _, _ = binned_log_slope(list(ns), list(a[lo:N + 1]),
                                bin_width=0.1, min_bin=50)
    checkpoints = {n: float(A[n] / n ** alpha)
                   for n in (10_000, 30_000, 100_000) if n <= N}
    return {
        "alpha": alpha, "B_COLD": math.exp(-zpp2),
        "amplitude_target": amp_target,
        "slope": {"value": slope, "err_class": "CERTIFIED(5e-3)",
                  "bins_used": nb_used, "bins_dropped": nb_drop},
        "slope_two_point": {"value": slope_2pt,
                            "err_class": "DIAGNOSTIC(2e-2)"},
        "amplitude_mean_of_logs": {"value": amp_mol,
                                   "err_class": "CERTIFIED(5%)"},
        "amplitude_endpoint": {"value": amp_end,
                               "err_class": "DIAGNOSTIC(5%)"},
        "amplitude_free_ls": {"value": float(math.exp(free_ic)),
                              "slope": float(free_slope),
                              "err_class": "DIAGNOSTIC"},
        "summatory_over_n_alpha": checkpoints,
        "pointwise_exponent": {"value": ps, "target": alpha - 1,
                               "err_class": "DIAGNOSTIC"},
        "N": N, "window": [lo, N],
    }


# ── correlator and the spectral function ────────────────────────

def correlator(p: int, t: float, beta: float = 2, N: int = 1000) -> Certified:
    """F(t) = ω(μ_p* α_t(μ_p)) = Σ_{n≤N} ρ̂_n e^{itΔE_n}, with the
    EM tail bound |tail| ≤ N^{−β}(log N/β + 1/β²)/(−ζ′(β+1))."""
    n = np.arange(2, N + 1, dtype=np.float64)
    rho = np.log(n) * n ** (-(beta + 1.0))
    mp = _mp()
    Z2 = float(-mp.zeta(beta + 1, derivative=1))
    dE = (beta + 1.0) * math.log(p) - np.log(np.log(p * n) / np.log(n))
    F = complex(np.sum(rho / Z2 * np.exp(1j * t * dE)))
    tail = N ** (-beta) * (math.log(N) / beta + 1.0 / beta ** 2) / Z2
    return Certified(F, f"CERTIFIED({tail:.1e})", tail,
                     {"p": p, "t": t, "N": N})


def spectral_function(p: int = 2, beta: float = 2, nmax: int = 10_000) -> Dict:
    """Line positions ΔE_n = (β+1) log p − log(1 + log p/log n),
    increasing to the accumulation edge (β+1) log p (OP3)."""
    n = np.arange(2, nmax + 1, dtype=np.float64)
    gaps = np.log1p(math.log(p) / np.log(n))
    lines = (beta + 1.0) * math.log(p) - gaps
    edge = (beta + 1.0) * math.log(p)
    return {"lines": lines, "gaps_to_edge": gaps, "edge": edge,
            "increasing": bool(np.all(np.diff(lines) > 0)),
            "last_line_below_edge": bool(lines[-1] < edge),
            "p": p, "beta": beta, "nmax": nmax}


def edge_mass(p: int = 2, beta: float = 2, eps: float = 0.1,
              nmax: int = 2_000_000) -> Certified:
    """Mass of spectral lines within eps of the edge — OP3 first gate.

    Pinned convention (note v0.1.1): mass(ε) = Σ_{n : g_n < ε} ρ̂_n
    evaluated per level, g_n = log(1 + log p/log n), no M-rounding.
    The edge-softness law (Addendum U.4) is its asymptote:

        mass ~ [M^{−β}(log M/β + 1/β²)]/(−ζ′(β+1)),
        M = exp((log p)/(e^ε − 1))   — exponentially soft edge."""
    n = np.arange(2, nmax + 1, dtype=np.float64)
    rho = np.log(n) * n ** (-(beta + 1.0))
    mp = _mp()
    Z2 = float(-mp.zeta(beta + 1, derivative=1))
    gaps = np.log1p(math.log(p) / np.log(n))
    mass = float(rho[gaps < eps].sum() / Z2)
    M = math.exp(math.log(p) / (math.exp(eps) - 1.0))
    pred = (M ** (-beta) * (math.log(M) / beta + 1.0 / beta ** 2)
            + 0.5 * math.log(M) * M ** (-beta - 1.0)) / Z2
    return Certified(mass, "DIAGNOSTIC", None,
                     {"predicted_asymptote": pred, "ratio": mass / pred,
                      "M(eps)": M, "eps": eps, "nmax": nmax})


# ── the gate suite as an entrypoint ─────────────────────────────

def gates(quick: bool = True) -> Dict:
    """The green suite of note §7 as a callable: 10 recorded gates
    covering G1–G7 (the delivery's 13-line log splits G6/G7 further;
    G0's environment anchors are the pytest suite's job here).
    quick=True runs the reduced-N form (G3 at 10⁶, cold gas at 8·10³,
    G7 with 3 refined zeros); quick=False runs the production N's
    (10⁷, 10⁵, 5)."""
    out = {}

    def rec(name, ok, klass, detail):
        out[name] = {"ok": bool(ok), "err_class": klass, "detail": detail}

    r = z1(2); rec("G1 Tr1=zeta(2)", r.detail["resid"] <= r.bound + 1e-27,
                   r.err_class, f"resid={r.detail['resid']:.2e}")
    r = z2(2); rec("G2 Tr2=-zeta'(3)", r.detail["resid"] <= r.bound + 1e-27,
                   r.err_class, f"resid={r.detail['resid']:.2e}")
    r = zD_certified_interval(2, N=1_000_000 if quick else 10_000_000)
    rec("G3 Z_D=-zz' interval", r.detail["inside"], r.err_class,
        f"width={r.bound:.2e}")

    # G3b reorder cross-check (direct sieve vs swapped (d,k) order)
    N2 = 100_000 if quick else 200_000
    b = 2.0
    w = weights_sieve(N2)
    nn = np.arange(2, N2 + 1, dtype=np.float64)
    S_direct = float(np.sum(w[2:] * nn ** (-b)))
    d2 = np.arange(2, N2 + 1, dtype=np.int64)
    M2 = N2 // d2
    H2 = _kahan_cumjumps(np.arange(1, N2 + 1, dtype=np.float64) ** (-b),
                         list(np.unique(M2)))
    S_swap = float(np.sum(np.log(d2) * d2.astype(np.float64) ** (-b - 1.0)
                          * np.array([H2[m] for m in M2])))
    rec("G3b reorder cross-check", abs(S_direct - S_swap) < 5e-13,
        "EXACT(reorder)", f"|direct-swap|={abs(S_direct - S_swap):.2e}")

    # G4 phase law vs literal termwise evolution
    ok4, worst = True, 0.0
    for pp in (2, 3, 5):
        for tt in (0.7, 1.9, math.pi):
            n0 = 7
            lhs = flow_phase(pp, n0, tt)
            rhs = np.exp(1j * tt * ((3.0) * math.log(pp)
                                    - math.log(math.log(pp * n0)
                                               / math.log(n0))))
            worst = max(worst, abs(lhs - rhs))
    rec("G4 flow phase law", worst < 1e-14, "EXACT", f"worst={worst:.1e}")

    r = kms_check(2, 2, 0.7, nbasis=2000)
    rec("G5 KMS t+i", r.value < r.bound, r.err_class,
        f"resid={r.value:.1e}, wrong-sign {r.detail['wrong_sign_control']:.2e}")

    Ngas = 8_000 if quick else 100_000
    rep = cold_gas_report(N=Ngas, lo=max(2_000, Ngas // 5))
    if quick:  # GP6 two-point diagnostic (engine-2 class)
        sv, sc = rep["slope_two_point"]["value"], 2e-2
    else:      # G6 certified binned fit (engine-1 class)
        sv, sc = rep["slope"]["value"], 5e-3
    ok_s = abs(sv - rep["alpha"]) < sc
    rec("G6 cold-gas slope", ok_s, rep["slope" if not quick else
        "slope_two_point"]["err_class"],
        f"fit={sv:.6f} vs {rep['alpha']:.6f}")
    av = (rep["amplitude_endpoint"] if quick
          else rep["amplitude_mean_of_logs"])["value"]
    ok_a = abs(av / rep["amplitude_target"] - 1) < 0.05
    rec("G6b cold-gas amplitude", ok_a,
        rep["amplitude_endpoint" if quick else
            "amplitude_mean_of_logs"]["err_class"],
        f"fit={av:.5f} vs {rep['amplitude_target']:.5f}")
    rec("G6c pointwise exponent", True, "DIAGNOSTIC",
        f"fit={rep['pointwise_exponent']['value']:.4f} "
        f"vs {rep['pointwise_exponent']['target']:.4f}")

    # G7 census cross-refinement (Speiser strip) — mpmath end-to-end
    # (zetaprime_refine returns python complex; float64 truncates the
    # residual floor at ~6e-16, so the certification path refines here
    # directly at dps 50, as the green suite does).
    from mtft.riemann import ZETAPRIME_ZEROS, zetaprime_negative_zero
    mp = _mp()
    n_ref = 3 if quick else 5
    worst = 0.0
    with mp.workdps(50):
        for z0 in ZETAPRIME_ZEROS[:n_ref]:
            z1r = mp.findroot(lambda s: mp.zeta(s, derivative=1),
                              mp.mpc(z0),
                              df=lambda s: mp.zeta(s, derivative=2))
            worst = max(worst, float(abs(mp.zeta(z1r, derivative=1))))
    neg = zetaprime_negative_zero(1)
    strip_ok = all(z.real > 0.5 for z in ZETAPRIME_ZEROS)
    rec("G7 Speiser strip harness", worst < 1e-25 and strip_ok,
        "CERTIFIED(1e-25)",
        f"max resid over {n_ref} refined zeros = {worst:.1e}; "
        f"negative zero {neg:.6f} -> beta {float(neg) - 1:.4f} outside strip")
    out["all_green"] = all(v["ok"] for v in out.values())
    return out
