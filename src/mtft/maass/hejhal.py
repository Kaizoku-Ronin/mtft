"""mtft.maass.hejhal — Hejhal's method for Maass cusp forms on Γ₀*(N) (CERTIFIED numerics).

A cusp form with spectral parameter R (eigenvalue λ = ¼ + R²), parity and Atkin–Lehner signs is

    f(x + iy) = Σ_n c_n √y K_{iR}(2π ν_n y) trig(2π ν_n x),      c_1 = 1,

trig = cos (even) or sin (odd), with frequencies ν_n = n, or n − ½ under the η-twist.  On a closed
horocycle y = Y below the Ford domain the coefficients satisfy the linear system

    c_n √Y K_{iR}(2π ν_n Y) = (2/Q) Σ_m χ_m f(z_m*) trig(2π ν_n x_m),

where z_m* is the sample point z_m carried into the domain and χ_m the character of that motion
(``mtft.maass.domain``).  R is an eigenvalue when the solutions from two horocycles agree; the
Hecke relations, which are never imposed, then hold as a check.

Two ways to look for eigenvalues are here, as in the lab's records MP15 and MP18:
``coefficients`` solves the normalised system (c_1 = 1) and ``refine`` finds the R at which a
coefficient agrees on two horocycles; ``sigma_min`` is the pole-free detector, the smallest
singular value of the homogeneous system on one horocycle, which dips at an eigenvalue.

A full scan of a range of R takes minutes and lives in the studies tier; the package refines from
the frozen table (``mtft.maass.spectrum``).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .domain import Horocycle, prime_factors
from .ford import lowest_vertex
from .kbessel import k_bessel

# Level 143, as in the records: the lowest point of the domain met by the lab's 4000 samples,
# the two horocycle heights 0.96 and 0.88 of it, and the truncation M = ⌈40/(2π Y0)⌉.
Y0_143 = 0.010083833810204635
FACTORS = (0.96, 0.88)
DETECTOR_SIZE = 300


@dataclass(frozen=True)
class Setup:
    N: int
    M: int
    horocycles: tuple

    @property
    def heights(self):
        return tuple(h.Y for h in self.horocycles)


def sample_count(M: int) -> int:
    """Sample points on a horocycle for M coefficients: the lab's rule ⌈1.08 M / 0.88⌉ + 12."""
    return math.ceil(1.08 * M / 0.88) + 12


@lru_cache(maxsize=8)
def setup(N: int = 143, M: int | None = None, heights: tuple | None = None) -> Setup:
    """Horocycles and truncation for the one-cusp group Γ₀*(N), N squarefree.

    By default the two horocycles sit at 0.96 and 0.88 of a height Y0 and the truncation is
    M = ⌈40/(2π Y0)⌉, so that K_{iR}(2π M Y) is of size e⁻⁴⁰.  Y0 is the lowest vertex of the Ford
    domain; at N = 143 it is the records' value instead (the lowest point their 4000 samples
    met, 0.010084 against the vertex's 0.009851), which gives the records' M = 632 and 788
    sample points.
    """
    if heights is None:
        Y0 = Y0_143 if N == 143 else lowest_vertex(N)[1]
        heights = tuple(f * Y0 for f in FACTORS)
    else:
        Y0 = max(heights) / max(FACTORS)
    if M is None:
        M = math.ceil(40.0 / (2 * math.pi * Y0))
    Q = sample_count(M)
    return Setup(N, M, tuple(Horocycle(N, Y, Q, True) for Y in heights))


def frequencies(M: int, twist: bool = False) -> np.ndarray:
    """ν_n = n, or n − ½ under the η-twist, n = 1 … M."""
    return np.arange(1, M + 1) - (0.5 if twist else 0.0)


def kcache(R: float, hc: Horocycle, fr: np.ndarray):
    """The two Bessel arrays a horocycle needs at R: at the reduced points, and on the horocycle."""
    return k_bessel(R, 2 * np.pi * np.outer(hc.ys, fr)), k_bessel(R, 2 * np.pi * fr * hc.Y)


def _raw_system(R, hc, fr, parity, chi, Kc):
    Kst, Kdiag = Kc if Kc is not None else kcache(R, hc, fr)
    M = fr.size
    trig = np.cos if parity == 0 else np.sin
    P = np.sqrt(hc.ys)[:, None] * Kst * trig(2 * np.pi * np.outer(hc.xs, fr))
    Tn = trig(2 * np.pi * np.outer(fr, hc.xm)) * (2.0 / hc.Q)
    V = (Tn * chi[None, :]) @ P
    V[np.arange(M), np.arange(M)] -= np.sqrt(hc.Y) * Kdiag
    return V


def system(R: float, hc: Horocycle, fr: np.ndarray, parity: int, chi: np.ndarray, Kc=None) -> np.ndarray:
    """The homogeneous system V c = 0 on one horocycle, columns scaled to unit size."""
    V = _raw_system(R, hc, fr, parity, chi, Kc)
    sc = np.abs(V).max(axis=0)
    sc[sc == 0] = 1.0
    return V / sc[None, :]


def coefficients(R: float, hc: Horocycle, fr: np.ndarray, parity: int, chi: np.ndarray, Kc=None) -> np.ndarray:
    """The coefficient vector with c_1 = 1 on one horocycle (the first equation is dropped)."""
    V = _raw_system(R, hc, fr, parity, chi, Kc)
    A = V[1:, 1:]
    sc = np.abs(A).max(axis=0)
    sc[sc == 0] = 1.0
    cr = np.linalg.solve(A / sc[None, :], -V[1:, 0]) / sc
    return np.concatenate([[1.0], cr])


def sigma_min(R: float, hc: Horocycle, fr: np.ndarray, parity: int, chi: np.ndarray, Kc=None) -> float:
    """σ_min/σ_max of the homogeneous system: the pole-free detector, small at an eigenvalue."""
    sv = np.linalg.svd(system(R, hc, fr, parity, chi, Kc), compute_uv=False)
    return float(sv[-1] / sv[0])


def detector(N: int = 143, factor: float = 0.96, size: int = DETECTOR_SIZE, twist: bool = False):
    """One horocycle and the frequencies for the detector: a 300-coefficient system at N = 143."""
    if N != 143:
        raise ValueError("the detector's defaults are for level 143")
    return Horocycle(N, factor * Y0_143, sample_count(size)), frequencies(size, twist)


def agreement(c1: np.ndarray, c2: np.ndarray, upto: int = 20) -> float:
    """How well two horocycles agree: the largest difference of c_2 … c_upto, relative to their size."""
    n = min(upto, c1.size, c2.size)
    return float(np.max(np.abs(c1[1:n] - c2[1:n])) / max(1.0, np.max(np.abs(c1[1:n]))))


def hecke_residual(c: np.ndarray, twist: bool = False) -> float:
    """The largest failure of five Hecke relations that the method never imposes.

    Untwisted (c[k] is c_{k+1}): c₂c₃ = c₆, c₂c₅ = c₁₀, c₄ = c₂² − 1, c₃c₅ = c₁₅, c₉ = c₃² − 1.
    Twisted (b_m is the coefficient of frequency m/2, m odd, so b_m = c[(m − 1)/2]):
    b₃b₅ = b₁₅, b₃b₇ = b₂₁, b₉ = b₃² − 1, b₂₅ = b₅² − 1, b₅b₇ = b₃₅.
    """
    if c.size < 18:
        raise ValueError("at least 18 coefficients are needed for the Hecke relations")
    if twist:
        return float(max(abs(c[1] * c[2] - c[7]), abs(c[1] * c[3] - c[10]), abs(c[4] - (c[1] ** 2 - 1)),
                         abs(c[12] - (c[2] ** 2 - 1)), abs(c[2] * c[3] - c[17])))
    return float(max(abs(c[1] * c[2] - c[5]), abs(c[1] * c[4] - c[9]), abs(c[3] - (c[1] ** 2 - 1)),
                     abs(c[2] * c[4] - c[14]), abs(c[8] - (c[2] ** 2 - 1))))


def solve_pair(R: float, parity: int, eps=(), twist: bool = False, st: Setup | None = None):
    """The coefficient vectors from the two horocycles of a setup at one R."""
    st = st or setup()
    fr = frequencies(st.M, twist)
    return [coefficients(R, hc, fr, parity, hc.chi(eps, twist)) for hc in st.horocycles]


def refine(R0: float, parity: int, eps=(), twist: bool = False, st: Setup | None = None, half_width: float = 1e-5,
           index: int = 1, tol: float = 2e-13, max_steps: int = 30, upto: int = 20) -> dict:
    """The eigenvalue near R0: the R at which coefficient ``index`` agrees on the two horocycles.

    Looks for a sign change of c[index](horocycle 1) − c[index](horocycle 2) on
    [R0 − half_width, R0 + half_width] and closes on it by the secant rule kept inside the
    bracket (Illinois), until a step is smaller than ``tol``.  Returns R, the coefficient vectors
    of both horocycles, their agreement on c_2 … c_upto and the Hecke residual (None when fewer
    than 18 coefficients are compared: small systems do not determine the high coefficients).
    A missing sign change raises ValueError: R0 is then not within half_width of an eigenvalue
    of this sector, or the bracket holds a pole of the normalised system.
    """
    st = st or setup()

    def gap(R):
        c1, c2 = solve_pair(R, parity, eps, twist, st)
        return c1[index] - c2[index], c1, c2

    a, b = R0 - half_width, R0 + half_width
    fa, fb = gap(a)[0], gap(b)[0]
    if not fa * fb < 0:
        raise ValueError(f"no sign change on [{a}, {b}] for this sector: {fa:.3e}, {fb:.3e}")
    c, side, steps = b, 0, 0
    for steps in range(1, max_steps + 1):
        c_new = (a * fb - b * fa) / (fb - fa)
        fc, c1, c2 = gap(c_new)
        done = abs(c_new - c) < tol or fc == 0
        c = c_new
        if done:
            break
        if fc * fb < 0:
            a, fa, side = b, fb, 0
        else:
            if side == 1:
                fa *= 0.5
            side = 1
        b, fb = c, fc
    return {"R": float(c), "lambda": 0.25 + float(c) ** 2, "coefficients": c1, "coefficients_2": c2,
            "agreement": agreement(c1, c2, upto), "hecke": hecke_residual(c1, twist) if min(upto, c1.size) >= 18 else None,
            "heights": st.heights, "M": st.M, "steps": steps, "primes": prime_factors(st.N)}
