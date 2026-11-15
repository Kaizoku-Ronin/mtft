"""mtft.maass.kbessel — the Bessel function K_{iR}(x) of imaginary order (CERTIFIED against mpmath).

    K_{iR}(x) = ∫₀^∞ exp(−x cosh t) cos(R t) dt,        R real, x > 0.

``kir`` is the fast route: with t = s/√(1 + x) the integrand is even and analytic in s, so the
trapezoid rule on a fixed grid converges geometrically.  It is validated against mpmath for
0 ≤ R ≤ 3 and 0.02 ≤ x ≤ 30 (2.1e-14 in the lab's record MP15, check V0; the gate here asks
1e-12) and returns 0 for x > 60, where K is below e⁻⁶⁰ of its values at the points that matter.
Above R = 3 the true value is exponentially small against the integrand and the rule loses
digits, so ``kir`` refuses R > 3; ``kir_mp`` (mpmath) serves there.  Every Maass form of
Γ₀(143) in the frozen table has R < 2.
"""
from __future__ import annotations

import numpy as np

R_MAX_FAST = 3.0
_S_GRID = np.arange(0.0, 11.0 + 1e-9, 0.2)
_W_S = np.full(_S_GRID.size, 0.2)
_W_S[0] = 0.1                                   # the trapezoid rule's half weight at s = 0


def kir(R: float, x) -> np.ndarray:
    """K_{iR}(x) for 0 ≤ R ≤ 3 and an array x > 0, by the trapezoid rule (double precision)."""
    if not 0.0 <= R <= R_MAX_FAST:
        raise ValueError(f"kir is validated for 0 <= R <= {R_MAX_FAST}; use kir_mp for R = {R}")
    x = np.asarray(x, dtype=float)
    out = np.zeros(x.shape)
    keep = x <= 60.0
    xs = x[keep]
    sc = 1.0 / np.sqrt(1.0 + xs)
    t = sc[:, None] * _S_GRID[None, :]
    out[keep] = (np.exp(-xs[:, None] * np.cosh(t)) * np.cos(R * t)) @ _W_S * sc
    return out


def kir_mp(R: float, x) -> np.ndarray:
    """K_{iR}(x) from mpmath, entry by entry: the reference route, and the one for R > 3."""
    import mpmath as mp
    x = np.asarray(x, dtype=float)
    return np.vectorize(lambda v: float(mp.besselk(1j * R, v).real), otypes=[float])(x)


def k_bessel(R: float, x) -> np.ndarray:
    """``kir`` where it is validated (R ≤ 3), ``kir_mp`` above."""
    return kir(R, x) if R <= R_MAX_FAST else kir_mp(R, x)
