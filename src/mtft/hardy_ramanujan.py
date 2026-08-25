"""mtft.hardy_ramanujan — an orthodox end-to-end benchmark (v0.20.0).

    modular form -> partition function -> cumulant geometry -> saddle
    -> integer combinatorics

This module claims NOTHING new.  It exists to exercise the whole MTFT stack on
textbook ground where the answer is independently known, so that any later
arithmetic density-of-states claim has a calibrated precedent.

The ensemble is the bosonic oscillator / integer-partition ensemble

    Z(beta) = prod_{m>=1} (1 - e^{-beta m})^{-1} = q^{1/24}/eta(tau),  q = e^{-beta}

Two independent routes to log Z:

  psi_direct    the convergent series sum_k 1/(k (e^{k beta} - 1))
  psi_modular   the eta modular transformation, exact:
                psi(b) = pi^2/(6b) - b/24 + (1/2) log(b/2pi) + psi(4 pi^2 / b)

These agree to full precision — modularity is an IDENTITY here, not an
asymptotic — which is the E2 pair for this module.

The four thermodynamic layers are the insertion calculus in one variable:
    psi = log Z,  E = -psi',  g = psi'' = Var,  T = psi''' = kappa_3.

Inverting by saddle point recovers Hardy-Ramanujan, and the notorious
prefactor is certified symbolically to be

    sqrt(beta_*/2pi) / sqrt(2 pi psi'')  ==  1/(4 sqrt(3) n),

i.e. (the modular half-log term) x (the Gaussian fluctuation determinant).
"""
from __future__ import annotations

import mpmath as mp
import sympy as sp

__all__ = [
    "psi_direct", "psi_modular", "modularity_residual", "layers",
    "saddle_partition", "hardy_ramanujan_asymptotic", "prefactor_identity",
]


def psi_direct(beta, dps=40):
    """log Z(beta) by the convergent series sum_k 1/(k (e^{k beta} - 1))."""
    beta = mp.mpf(beta)
    if beta <= 0:
        raise ValueError("beta must be positive")
    s, k = mp.mpf(0), 1
    tol = mp.mpf(10) ** (-(mp.mp.dps + 5))
    while True:
        t = 1 / (k * mp.expm1(k * beta))
        s += t
        if t < tol * max(1, s) or k > 200000:
            break
        k += 1
    return s


def psi_modular(beta, dps=40):
    """log Z(beta) via the eta modular transformation (exact, not asymptotic)."""
    beta = mp.mpf(beta)
    dual = 4 * mp.pi ** 2 / beta
    return (mp.pi ** 2 / (6 * beta) - beta / 24
            + mp.log(beta / (2 * mp.pi)) / 2
            + psi_direct(dual, dps))


def modularity_residual(beta, dps=40):
    """E2 check: |psi_direct - psi_modular|.  Should be at machine level."""
    with mp.workdps(dps):
        return abs(psi_direct(beta, dps) - psi_modular(beta, dps))


def layers(beta, dps=30):
    """The four thermodynamic layers (psi, E, g, T) at a given beta."""
    with mp.workdps(dps + 10):
        b = mp.mpf(beta)
        f = lambda x: psi_modular(x, dps)
        return {
            "psi": f(b),
            "E": -mp.diff(f, b),
            "g": mp.diff(f, b, 2),
            "T": mp.diff(f, b, 3),
        }


def saddle_partition(n, dps=30):
    """p(n) by saddle point on the exact psi.  Returns (value, beta_star, layers).

    The saddle solves n = -psi'(beta); the Gaussian fluctuation gives
        p(n) ~ exp(psi(b*) + n b*) / sqrt(2 pi psi''(b*)).
    Because psi is exact rather than truncated, this resums the Hardy-Ramanujan
    closed form and is roughly 3x more accurate at the n tested.
    """
    with mp.workdps(dps + 10):
        n = mp.mpf(n)
        f = lambda x: psi_modular(x, dps)
        b0 = mp.pi / mp.sqrt(6 * n)
        # bracket the saddle, then bisect-then-polish: findroot's secant can
        # step to negative beta where psi is undefined.
        E = lambda b: -mp.diff(f, b)
        lo, hi = b0 / 4, b0 * 4
        while E(lo) < n:
            lo /= 2
        while E(hi) > n:
            hi *= 2
        for _ in range(80):
            mid = (lo + hi) / 2
            if E(mid) > n:
                lo = mid
            else:
                hi = mid
        bs = (lo + hi) / 2
        g = mp.diff(f, bs, 2)
        val = mp.e ** (f(bs) + n * bs) / mp.sqrt(2 * mp.pi * g)
        return +val, +bs, {"psi": f(bs), "E": n, "g": +g,
                           "T": +mp.diff(f, bs, 3)}


def hardy_ramanujan_asymptotic(n):
    """The closed form  p(n) ~ e^{2 pi sqrt(n/6)} / (4 n sqrt 3)."""
    n = mp.mpf(n)
    return mp.e ** (2 * mp.pi * mp.sqrt(n / 6)) / (4 * n * mp.sqrt(3))


def prefactor_identity():
    """Symbolic certificate that the HR prefactor is the Gaussian determinant.

    With beta_* = pi/sqrt(6n) and the leading psi'' = pi^2/(3 beta^3),

        sqrt(beta_*/2pi) / sqrt(2 pi psi'')  ==  1/(4 sqrt(3) n).

    Returns (expression, residual); residual is 0 on success.
    """
    n = sp.Symbol("n", positive=True)
    bstar = sp.pi / sp.sqrt(6 * n)
    pref = sp.sqrt(bstar / (2 * sp.pi)) / sp.sqrt(2 * sp.pi * (sp.pi ** 2 / (3 * bstar ** 3)))
    return sp.simplify(pref), sp.simplify(pref - 1 / (4 * sp.sqrt(3) * n))
