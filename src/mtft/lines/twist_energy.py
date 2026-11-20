r"""mtft.lines.twist_energy — the twist landscape: the relative heat trace of a twist, and the energy and pin it gives.

EXACT within VE-01's product approximation; TE-01, 2026-11-20.  From records FF-14 (the heat trace and the pin) and FF-15 (the
energy and the transform), with FF-20's bands, in the MTFT Python project under claude/experiments/.

A twist ζ ∈ H¹(X₀(143), ℤ/2) (``twists``) lifts one fermion zero mode of every adjoint component on the curve (MATTER-04).  On
the Selberg side, the bosons' heat trace minus the twisted fermions' is

    D_ζ(t̃) = Σ_{γ₀} Σ_{k odd} 2 [ζ(γ₀) = −1] ℓ₀ e^{−t̃/4} e^{−(kℓ₀)²/4t̃} / (√(4πt̃) · 2 sinh(kℓ₀/2)),    t̃ = t/R_C²,

a sum over the primitive oriented closed geodesics γ₀, of length ℓ₀, that ζ makes antiperiodic.  The identity, cusp and
continuum terms cancel, because ζ is trivial on parabolic elements.  Every term is ≥ 0, so a truncation at L_max is a lower
bound at every t̃.  The central value adds a ζ-independent tail for the geodesics beyond L_max, with ζ = −1 on half of them (the
prime geodesic theorem): ½ erfc((L_max − t̃)/(2√t̃)).  The central value tends to 1 as t̃ → ∞: the lifted zero mode.

Through VE-01's sum over the windings w of the ledger's torus (n even), D_ζ gives

  energy   V(ζ) = −32 (4π/1536 r R⁶) ∫ t̃⁻⁴ (S₀ + 1) D_ζ dt̃: the one-loop energy at the pattern, in ledger units (FF-15; the
           "+ 1" is the w = 0 winding)
  pin      pin_a(ζ) = 64π² (4π/1536 r R⁶) ∫ t̃⁻⁴ S_a D_ζ dt̃: the twelve's curvature on cycle a = 1 or τ (FF-14, in FF-12's
           per-cycle convention)

with R² = r/56π, S_a = Σ_{n even} w_a² e^{−|w|²/4R²t̃} and S₀ = Σ_{n even, w ≠ 0} e^{−|w|²/4R²t̃}.  Both are linear in which
classes are twisted.  So each class carries a weight, and every twist's value at once is a Walsh–Hadamard transform over F₂²⁶:
Σ_y μ(y) [ζ·y odd] = (Σμ − F(ζ))/2, with F the transform of the class measure μ.  Both depend only on the lengths of the twisted
geodesics, so they are the same on a twist and on its images under W₁₁, W₁₃, W₁₄₃ and the mirror, which permute the geodesics
and keep their lengths: the landscape's extremes come in orbits (``summary`` lists them).

Layers:

  classes, per_trace_check          the geodesics up to L_max, with their 26 class bits.  Second route: Σ ℓ₀ per trace against
                                    the frozen weights of ``mtft.maass.trace`` (FF-14 T2)
  twisted                           which classes a twist makes antiperiodic
  heat_terms, tail, relative_heat_trace, mean_relative_heat_trace, J
                                    D_ζ on the t̃ grid, as a lower bound and a central value
  lattice_moments, energy_weights, pin_weights
                                    the torus sums, the per-class weights and their tails at a radius r
  energy, pin                       one twist, by direct sums
  fwht, landscape, summary          every twist, by the transform (2²⁶ doubles: 512 MB)
  bands                             the depth V/V_min by the number of systoles twisted (FF-20: five bands that do not overlap)
  pin_matrix, hold                  from a D curve: the twelve's curvature eigenvalues, and the hold of τ (V, its gradient in τ
                                    and sinh d; FF-14 T4)

Conventions.  The energy includes the w = 0 winding (FF-15).  FF-14's hold printed V without it, because
``twist_costs.windings`` has no (0, 0); ``hold`` returns V both ways.  The gradient and sinh d do not depend on it.  The t̃ grid
stops at 600; FF-14 found that the omitted part only raises the values (by 0.25% for the mean twist).

This is not the glueball fit's "constant twist energy" (``mtft.maass.glueballs.twist_energy``, GB-02), which is a different
quantity.  Not here: the twelve's masses, which need FF-13's Wilson-line masses (not in the package).
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from . import twists as _tw
from . import walls as _wl

R_WINDOW = (230.0, 280.0)              # the radii of the records
LMAX = 14.0                            # the records' cut on geodesic length
TG = np.exp(np.linspace(math.log(0.02), math.log(600.0), 2400))      # the t̃ grid (FF-14)
TG.setflags(write=False)
_DLT = np.log(TG[1]) - np.log(TG[0])
KODD = (1, 3, 5, 7)                    # the odd multiples kept in each class's heat term
SWITCH = 4.0                           # above this t̃ the torus sums take their continuum values
NW = 70                                # windings with |m|, |n| ≤ NW


# ------------------------------------------------------------------ the classes
@lru_cache(maxsize=4)
def classes(lmax: float = LMAX) -> tuple:
    """(lengths ℓ₀, class bits Y, cycle lengths ℓ over the level-1 class) of every primitive oriented closed geodesic up to lmax."""
    P = _wl.primitive_classes(lmax)
    L0 = np.array([p[0] for p in P], dtype=float)
    Y = np.array([_tw.class_bits(p[2]) for p in P], dtype=np.int64)
    ELL = np.array([p[1] for p in P], dtype=np.int64)
    for a in (L0, Y, ELL):
        a.setflags(write=False)
    return L0, Y, ELL


def _chebyshev(t0: int, k: int) -> list:
    """Traces of A, A², …, A^k for tr A = t0."""
    out, a, b = [t0], 2, t0
    for _ in range(k - 1):
        a, b = b, t0 * b - a
        out.append(b)
    return out


def per_trace_check(lmax: float = LMAX) -> dict:
    """Second route to the classes: per trace t, Σ ℓ₀ over every class and its powers, against the frozen weights w(t)√(t² − 4)
    of ``mtft.maass.trace`` at level 143 (FF-14 T2)."""
    from mtft.maass import trace as TR
    L0, _, ELL = classes(lmax)
    tcap = math.floor(2 * math.cosh(lmax / 2))
    acc = {}
    for l0, ell in zip(L0.tolist(), ELL.tolist()):
        t = _chebyshev(round(2 * math.cosh(l0 / (2 * ell))), ell)[-1]
        for tk in _chebyshev(t, int(math.acosh(tcap / 2) * 2 / l0) + 1):
            if tk <= tcap:
                acc[tk] = acc.get(tk, 0.0) + l0
    s, _, w = TR.classes(143, 1, "hyp")
    ratios, both, one = [], 0, 0
    for tt, ww in zip(s.tolist(), w.tolist()):
        tt = int(tt)
        if tt > tcap:
            continue
        got = acc.get(tt, 0.0)
        if ww <= 0:
            both, one = (both + 1, one) if got == 0 else (both, one + 1)
            continue
        ratios.append(got / (ww * math.sqrt(tt * tt - 4)))
    return {"trace bound": tcap, "traces with classes": len(ratios), "empty traces (both routes)": both,
            "empty in one route only": one, "traces found here but absent from the frozen table": len(set(acc) - {int(x) for x in s}),
            "max |ratio - 1|": max(abs(x - 1) for x in ratios) if ratios else math.nan}


def _popcount32(v: np.ndarray) -> np.ndarray:
    v = v - ((v >> np.uint32(1)) & np.uint32(0x55555555))
    v = (v & np.uint32(0x33333333)) + ((v >> np.uint32(2)) & np.uint32(0x33333333))
    v = (v + (v >> np.uint32(4))) & np.uint32(0x0F0F0F0F)
    return (v * np.uint32(0x01010101)) >> np.uint32(24)


def twisted(z: int, Y) -> np.ndarray:
    """Boolean mask over the classes: those the twist z makes antiperiodic (ζ(γ) = −1)."""
    return (_popcount32(np.asarray(Y).astype(np.uint32) & np.uint32(z)) & np.uint32(1)).astype(bool)


# ------------------------------------------------------------------ the relative heat trace
def heat_terms(L0, kodd=KODD) -> np.ndarray:
    """h(t̃; ℓ₀) = 2 Σ_{k odd} ℓ₀ e^{−t̃/4} e^{−(kℓ₀)²/4t̃} / (√(4πt̃) 2 sinh(kℓ₀/2)) on the grid, one row per length."""
    L0 = np.asarray(L0, float)[:, None]
    out = np.zeros((L0.shape[0], len(TG)))
    pref = np.exp(-TG / 4) / np.sqrt(4 * math.pi * TG)
    for k in kodd:
        out += 2 * L0 * pref * np.exp(-(k * L0) ** 2 / (4 * TG)) / (2 * np.sinh(k * L0 / 2))
    return out


def tail(lmax: float = LMAX) -> np.ndarray:
    """The ζ-independent tail of D from the geodesics longer than lmax: ½ erfc((lmax − t̃)/(2√t̃))."""
    from scipy.special import erfc
    return 0.5 * erfc((lmax - TG) / (2 * np.sqrt(TG)))


def integrate(f):
    """∫ f dt̃ on the log grid (the grid is the last axis of f)."""
    g = f * TG
    return (g[..., 1:] + g[..., :-1]).sum(-1) * _DLT / 2


def _sum_heat(L0) -> np.ndarray:
    """Σ h(t̃; ℓ₀) over the given lengths, computed once per distinct length."""
    uL, cnt = np.unique(np.round(L0, 12), return_counts=True)
    out = np.zeros(len(TG))
    for i0 in range(0, len(uL), 4000):
        out += (cnt[i0:i0 + 4000, None] * heat_terms(uL[i0:i0 + 4000])).sum(0)
    return out


def relative_heat_trace(z: int, lmax: float = LMAX) -> tuple:
    """(lower bound, central value) of D_ζ on the grid for the twist z.  The untwisted O has D = 0 exactly."""
    if z == 0:
        return np.zeros(len(TG)), np.zeros(len(TG))
    L0, Y, _ = classes(lmax)
    lower = _sum_heat(L0[twisted(z, Y)])
    return lower, lower + tail(lmax)


def mean_relative_heat_trace(lmax: float = LMAX) -> tuple:
    """(lower bound, central value) of D averaged over the 2²⁶ − 1 twists other than O: each class with Y ≠ 0 is twisted by a
    fraction 2²⁵/(2²⁶ − 1) of them."""
    L0, Y, _ = classes(lmax)
    d = _tw.dim()
    lower = (1 << (d - 1)) / ((1 << d) - 1) * _sum_heat(L0[Y != 0])
    return lower, lower + tail(lmax)


def J(D: np.ndarray) -> float:
    """J = ∫ t̃⁻² D dt̃ (FF-14)."""
    return float(integrate(TG ** -2 * D))


# ------------------------------------------------------------------ the torus sums and the weights
def R2(r: float) -> float:
    """R_C² = r/56π."""
    return r / (56 * math.pi)


def _cv(r: float) -> float:
    return 4 * math.pi / (1536 * r * R2(r) ** 3)


@lru_cache(maxsize=8)
def lattice_moments(r: float) -> tuple:
    """(S, S₀) on the grid: S[:, j] = Σ_{n even} (m², n², mn)_j e^{−|w|²/4R²t̃} and S₀ = Σ_{n even, w ≠ 0} e^{−|w|²/4R²t̃}, summed
    directly up to t̃ = 4 and taken at their continuum values above (FF-14)."""
    from mtft.matter import twist_costs as TC
    tau = TC.tau()
    m, n, w2 = TC.windings(NW, tau)
    ev = np.rint(n).astype(int) % 2 == 0
    m, n, w2 = m[ev], n[ev], w2[ev]
    s = 4 * R2(r) * TG
    S, S0 = np.zeros((len(TG), 3)), np.zeros(len(TG))
    ginv = np.array([[abs(tau) ** 2, -tau.real], [-tau.real, 1.0]]) / tau.imag
    for i, si in enumerate(s):
        if TG[i] <= SWITCH:
            e = np.exp(-w2 / si)
            S[i] = ((m * m * e).sum(), (n * n * e).sum(), (m * n * e).sum())
            S0[i] = e.sum()
        else:
            S[i] = (math.pi * si * si / 4) * np.array([ginv[0, 0], ginv[1, 1], ginv[0, 1]])
            S0[i] = math.pi * si / 2 - 1
    S.setflags(write=False)
    S0.setflags(write=False)
    return S, S0


def switch_check(r: float) -> float:
    """The largest relative gap between the direct sums S₁₁, S₂₂ and their continuum values at the switch t̃ = 4."""
    from mtft.matter import twist_costs as TC
    tau = TC.tau()
    m, n, w2 = TC.windings(NW, tau)
    ev = np.rint(n).astype(int) % 2 == 0
    si = 4 * R2(r) * SWITCH
    e = np.exp(-w2[ev] / si)
    ginv = np.array([[abs(tau) ** 2, -tau.real], [-tau.real, 1.0]]) / tau.imag
    d = np.array([(m[ev] ** 2 * e).sum(), (n[ev] ** 2 * e).sum()])
    c = (math.pi * si * si / 4) * np.array([ginv[0, 0], ginv[1, 1]])
    return float(np.abs(d / c - 1).max())


def _class_integrals(lmax: float, kernels) -> np.ndarray:
    """∫ t̃⁻⁴ K(t̃) h(t̃; ℓ₀) dt̃ for every class and every kernel K, computed once per distinct length."""
    L0, _, _ = classes(lmax)
    uL, inv = np.unique(np.round(L0, 12), return_inverse=True)
    out = np.zeros((len(uL), len(kernels)))
    for i0 in range(0, len(uL), 4000):
        hh = heat_terms(uL[i0:i0 + 4000])
        for j, K in enumerate(kernels):
            out[i0:i0 + 4000, j] = integrate(TG ** -4 * K * hh)
    return out[inv.ravel()]


@lru_cache(maxsize=8)
def energy_weights(r: float = 230.0, lmax: float = LMAX) -> tuple:
    """(weight of every class, tail): −32 (4π/1536 r R⁶) ∫ t̃⁻⁴ (S₀ + 1) h dt̃, and the same over the tail of D (FF-15)."""
    _, S0 = lattice_moments(r)
    k = -32 * _cv(r)
    w = k * _class_integrals(lmax, (S0 + 1,))[:, 0]
    w.setflags(write=False)
    return w, float(k * integrate(TG ** -4 * (S0 + 1) * tail(lmax)))


@lru_cache(maxsize=8)
def pin_weights(r: float = 230.0, lmax: float = LMAX) -> tuple:
    """(weights, one column per cycle (1, τ), and the two tails): 64π² (4π/1536 r R⁶) ∫ t̃⁻⁴ S_a h dt̃ (FF-14)."""
    S, _ = lattice_moments(r)
    c = 64 * math.pi ** 2 * _cv(r)
    P = c * _class_integrals(lmax, (S[:, 0], S[:, 1]))
    P.setflags(write=False)
    return P, tuple(float(c * integrate(TG ** -4 * S[:, a] * tail(lmax))) for a in (0, 1))


# ------------------------------------------------------------------ one twist
def energy(z: int, r: float = 230.0, lmax: float = LMAX, central: bool = True) -> float:
    """V(ζ) in ledger units: the sum of the weights of the classes z twists, plus the tail when central (0 for O)."""
    if z == 0:
        return 0.0
    w, t = energy_weights(r, lmax)
    v = float(w[twisted(z, classes(lmax)[1])].sum())
    return v + t if central else v


def pin(z: int, r: float = 230.0, cycle: str = "1", lmax: float = LMAX, central: bool = True) -> float:
    """pin_a(ζ) on cycle "1" or "tau" (FF-12's per-cycle convention; 0 for O)."""
    if z == 0:
        return 0.0
    a = {"1": 0, "tau": 1}[cycle]
    P, tails = pin_weights(r, lmax)
    v = float(P[twisted(z, classes(lmax)[1]), a].sum())
    return v + tails[a] if central else v


# ------------------------------------------------------------------ every twist
def fwht(a: np.ndarray) -> np.ndarray:
    """The Walsh–Hadamard transform in place: F[x] = Σ_y (−1)^{x·y} a[y] (len(a) a power of 2)."""
    n, h = a.size, 1
    if n & (n - 1):
        raise ValueError("the length must be a power of 2")
    while h < n:
        v = a.reshape(-1, 2, h)
        x = v[:, 0, :].copy()
        v[:, 0, :] += v[:, 1, :]
        v[:, 1, :] = x - v[:, 1, :]
        h *= 2
    return a


def landscape(r: float = 230.0, quantity: str = "energy", central: bool = True, lmax: float = LMAX) -> np.ndarray:
    """The energy ("energy") or a pin ("pin 1", "pin tau") of every twist at once, indexed by the twist's 26 bits.  Entry 0, the
    untwisted O, is nan.  2²⁶ doubles: 512 MB, and about 30 s."""
    _, Y, _ = classes(lmax)
    if quantity == "energy":
        w, t = energy_weights(r, lmax)
    elif quantity in ("pin 1", "pin tau"):
        P, tails = pin_weights(r, lmax)
        a = 0 if quantity == "pin 1" else 1
        w, t = P[:, a], tails[a]
    else:
        raise ValueError('quantity is "energy", "pin 1" or "pin tau"')
    mu = np.zeros(1 << _tw.dim())
    np.add.at(mu, Y, w)
    total = mu.sum()
    v = fwht(mu)
    np.subtract(total, v, out=v)
    v *= 0.5
    if central:
        v += t
    v[0] = np.nan
    return v


def summary(v: np.ndarray, tie: float = 1e-12) -> dict:
    """Min, mean, max, the twists where they are, and quantiles of a landscape (nan entries ignored).

    The curve's involutions permute its closed geodesics without changing their lengths, so the energy and the pins are the
    same on a twist and on its images, and the extremes come in orbits.  Values within a relative ``tie`` of an extreme are
    listed (up to 64); "argmin" and "argmax" are the smallest twists among them, not an accident of summation order."""
    mn, mx = float(np.nanmin(v)), float(np.nanmax(v))
    at_min = np.flatnonzero(np.abs(v - mn) <= tie * abs(mn))
    at_max = np.flatnonzero(np.abs(v - mx) <= tie * abs(mx))
    q = np.nanpercentile(v, [0.1, 1, 10, 50, 90, 99, 99.9])
    return {"min": mn, "mean": float(np.nanmean(v)), "max": mx, "argmin": int(at_min[0]), "argmax": int(at_max[0]),
            "twists at the min": at_min[:64].tolist(), "twists at the max": at_max[:64].tolist(),
            "quantiles (0.1, 1, 10, 50, 90, 99, 99.9 %)": q.tolist()}


def bands(V: np.ndarray, systole_classes=None, chunk: int = 1 << 22) -> dict:
    """The depth V/V_min of every twist, by the number of systoles it twists (FF-20).  V is a landscape (nan entries skipped);
    the systoles default to those of ``walls``."""
    S = _wl.systoles() if systole_classes is None else tuple(systole_classes)
    vmin = float(np.nanmin(V))
    acc = {j: [0, math.inf, -math.inf] for j in range(len(S) + 1)}
    for i0 in range(0, V.size, chunk):
        z = np.arange(i0, min(i0 + chunk, V.size), dtype=np.uint32)
        k = np.zeros(z.size, dtype=np.int8)
        for y in S:
            k += (_popcount32(z & np.uint32(y)) & np.uint32(1)).astype(np.int8)
        d = V[i0:i0 + z.size] / vmin
        for j, a in acc.items():
            sel = (k == j) & np.isfinite(d)
            if sel.any():
                a[0] += int(sel.sum())
                a[1] = min(a[1], float(d[sel].min()))
                a[2] = max(a[2], float(d[sel].max()))
    out = {j: {"twists": c, "depth min": lo, "depth max": hi} for j, (c, lo, hi) in acc.items()}
    srt = sorted((b["depth min"], b["depth max"]) for b in out.values() if b["twists"])
    return {"bands": out, "overlap": any(srt[i][1] > srt[i + 1][0] for i in range(len(srt) - 1))}


# ------------------------------------------------------------------ the model's other outputs, from a D curve
def pin_matrix(D: np.ndarray, r: float = 230.0) -> dict:
    """The per-cycle pins and the twelve's curvature eigenvalues 64π² eig(M_A), M_A = Σ_{n even} (κ/2) w wᵀ (FF-14 T5)."""
    S, _ = lattice_moments(r)
    c = _cv(r)
    M = np.array([[c * integrate(TG ** -4 * S[:, 0] * D), c * integrate(TG ** -4 * S[:, 2] * D)],
                  [c * integrate(TG ** -4 * S[:, 2] * D), c * integrate(TG ** -4 * S[:, 1] * D)]])
    return {"per-cycle pins": [64 * math.pi ** 2 * M[0, 0], 64 * math.pi ** 2 * M[1, 1]],
            "twelve: 64π² eig(M_A)": (64 * math.pi ** 2 * np.linalg.eigvalsh(M / 2)).tolist()}


def _V_and_grad(D: np.ndarray, r: float, tau: complex, gradient: bool = True) -> tuple:
    """V without the w = 0 winding (FF-14's convention) and its gradient in (Re τ, Im τ)."""
    from mtft.matter import twist_costs as TC
    m, n, w2 = TC.windings(NW, tau)
    ev = np.rint(n).astype(int) % 2 == 0
    m, n, w2 = m[ev], n[ev], w2[ev]
    x, y = tau.real, tau.imag
    dx = 2 * n * (m + n * x) / y                                    # ∂|w|²/∂Re τ, with |w|² = ((m + n x)² + n² y²)/y
    dy = (2 * n * n * y * y - ((m + n * x) ** 2 + n * n * y * y)) / (y * y)
    s = 4 * R2(r) * TG
    S0, Sx, Sy = np.zeros(len(TG)), np.zeros(len(TG)), np.zeros(len(TG))
    for i, si in enumerate(s):
        if TG[i] <= SWITCH:
            e = np.exp(-w2 / si)
            S0[i] = e.sum()
            if gradient:
                Sx[i] = (-dx / si * e).sum()
                Sy[i] = (-dy / si * e).sum()
        else:
            S0[i] = math.pi * si / 2 - 1
    k = -32 * _cv(r)
    return float(k * integrate(TG ** -4 * S0 * D)), float(k * integrate(TG ** -4 * Sx * D)), float(k * integrate(TG ** -4 * Sy * D))


def hold(D: np.ndarray, r: float = 230.0, eps: float = 1e-4) -> dict:
    """The hold of τ at the pattern (FF-14 T4): V both ways, |∇V| in the hyperbolic metric (analytic, and by central
    differences) and sinh d = c |∇V| / V_γ0 for g₄² = 0.31 and 0.37 (``twist_costs.c_factor``, ``well``)."""
    from mtft.matter import twist_costs as TC
    tau = TC.tau()
    V, gx, gy = _V_and_grad(D, r, tau)
    ghyp = tau.imag * math.hypot(gx, gy)
    Vx = (_V_and_grad(D, r, tau + eps, False)[0] - _V_and_grad(D, r, tau - eps, False)[0]) / (2 * eps)
    Vy = (_V_and_grad(D, r, tau + 1j * eps, False)[0] - _V_and_grad(D, r, tau - 1j * eps, False)[0]) / (2 * eps)
    zero = -32 * _cv(r) * integrate(TG ** -4 * D)
    return {"V (w = 0 included)": V + float(zero), "V (without w = 0, FF-14's print)": V, "|grad V|_hyp (analytic)": ghyp,
            "|grad V|_hyp (differences)": tau.imag * math.hypot(Vx, Vy),
            "sinh d": {g2: TC.c_factor(r) * ghyp / TC.well(r, g2) for g2 in (0.31, 0.37)}}
