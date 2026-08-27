"""Bergman harmonic channels on X_0(143).

Writing the alpha-normalized differentials as omega(z) = sum a_n q^n, the
bilinear coefficients B_{n,m} = a_n^dagger (Im tau)^{-1} a_m organize the
canonical Bergman density into horizontal Fourier channels

    C_k(y) = sum_{m>=1} B_{m+k,m} e^{-2 pi (2m+k) y},

with  B(x+iy) = C_0(y) + 2 sum_{k>=1} Re[C_k(y) e^{2 pi i k x}].

This is an exact restructuring of ``forms.bergman_density`` (no FFT), and it
explains the observed mode phenomenology: k=4 has a clean single-sign leading
pair (B_{5,1} = 3.9497 vs B_{6,2} = -0.1760), while k=1 suffers destructive
interference (B_{2,1} = -0.3002 against B_{3,2} = +2.3015, opposite signs,
ratio ~7.7), making it strongly pre-asymptotic.  The k=4 -> k=1 dominance
crossover is at y*/y0 = 2.302140221833918907 (y0 = 1/sqrt(143)).

Classes: CERTIFIED(tol) at the stated dps; mode-dominance statements are
properties of the canonical density, not physical claims.
"""
from __future__ import annotations

from functools import lru_cache

import mpmath as mp

from .core import GENUS, omega_symplectic, riemann_matrix
from .forms import raw_qexpansions

CROSSOVER_RATIO_41 = "2.302140221833918907"


@lru_cache(maxsize=4)
def _setup(nmax: int = 140, dps: int = 45):
    with mp.workdps(dps):
        raw = raw_qexpansions(nmax, dps)
        O = omega_symplectic(dps)
        A = O[:, :GENUS]
        t = riemann_matrix(dps)
        Y = mp.matrix([[mp.im(t[i, j]) for j in range(GENUS)]
                       for i in range(GENUS)])
        Yi = Y ** -1
        Ai = A ** -1
        a = tuple(+(Ai * mp.matrix([raw[i, n] for i in range(GENUS)]))
                  for n in range(nmax + 1))
        return a, Yi


def bergman_bilinear(n: int, m: int, nmax: int = 140, dps: int = 45):
    """B_{n,m} = a_n^dagger (Im tau)^{-1} a_m (mpc; real to tolerance)."""
    a, Yi = _setup(nmax, dps)
    with mp.workdps(dps):
        return +(a[n].conjugate().T * Yi * a[m])[0]


def bergman_channel(k: int, y, nmax: int = 140, dps: int = 45):
    """C_k(y) = sum_m B_{m+k,m} e^{-2 pi (2m+k) y}."""
    with mp.workdps(dps):
        y = mp.mpf(y)
        return +sum(bergman_bilinear(m + k, m, nmax, dps)
                    * mp.e ** (-2 * mp.pi * (2 * m + k) * y)
                    for m in range(1, nmax + 1 - k))


def channel_density(z, kmax: int = 80, nmax: int = 140, dps: int = 45):
    """Bergman density via the channel series; equals forms.bergman_density.

    Truncation in k contributes O(|B| e^{-2 pi kmax y}) with polynomially
    growing bilinear coefficients; kmax=80 suffices for ~1e-24 at natural
    heights, and kmax = nmax-1 recovers the exact rearrangement.
    """
    with mp.workdps(dps):
        z = mp.mpc(z)
        x, y = mp.re(z), mp.im(z)
        out = mp.re(bergman_channel(0, y, nmax, dps))
        for k in range(1, kmax + 1):
            out += 2 * mp.re(bergman_channel(k, y, nmax, dps)
                             * mp.e ** (2j * mp.pi * k * x))
        return +out


def mode_crossover(k_hi: int = 4, k_lo: int = 1, bracket=("2.28", "2.32"),
                   nmax: int = 140, dps: int = 45, iters: int = 140):
    """Solve |C_{k_hi}(y)| = |C_{k_lo}(y)| by bisection; returns y*/y0.

    Default (4,1) anchor: 2.302140221833918907, inside the observed
    dominance bracket [2.28, 2.32].
    """
    with mp.workdps(dps):
        y0 = 1 / mp.sqrt(143)
        f = lambda y: (abs(bergman_channel(k_hi, y, nmax, dps))
                       - abs(bergman_channel(k_lo, y, nmax, dps)))
        lo, hi = mp.mpf(bracket[0]) * y0, mp.mpf(bracket[1]) * y0
        flo, fhi = f(lo), f(hi)
        if not (flo > 0 > fhi or flo < 0 < fhi):
            raise ValueError("bracket does not straddle the crossover")
        for _ in range(iters):
            mid = (lo + hi) / 2
            if (f(mid) > 0) == (flo > 0):
                lo = mid
            else:
                hi = mid
        return +(((lo + hi) / 2) / y0)


__all__ = ["bergman_bilinear", "bergman_channel", "channel_density",
           "mode_crossover", "CROSSOVER_RATIO_41"]
