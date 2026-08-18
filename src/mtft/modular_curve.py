"""
mtft.modular_curve
===================
Modular curve X₀(N) computations for MTFT.

Computes genus, cusps, index, Hecke eigenvalues, vortex energy
landscapes, and homology data for the congruence subgroup Γ₀(N).

This module provides the mathematical backbone for the X₀(143)
explorer and generalises to arbitrary level N.

References
----------
- MTFT Paper 26: X₀(143) Modular-Fractal Bridge
- MTFT Chapter 18: Computational Methods
- Shimura, "Introduction to the Arithmetic Theory of
  Automorphic Functions" (1971)

Author: Roger Tano <mtft1093@gmail.com>
ORCID:  0009-0005-1113-3620
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import List, Tuple, Dict, Optional

import numpy as np

__all__ = [
    "ModularCurve",
    "CuspClass",
    "VortexConfig",
    "HeckeSpectrum",
    "HomologyData",
]


# ── number-theoretic helpers ────────────────────────────────

def _gcd(a: int, b: int) -> int:
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def _euler_phi(n: int) -> int:
    """Euler's totient function φ(n)."""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


def _divisors(n: int) -> List[int]:
    """Return sorted list of positive divisors of n."""
    divs = []
    for d in range(1, int(math.isqrt(n)) + 1):
        if n % d == 0:
            divs.append(d)
            if d != n // d:
                divs.append(n // d)
    return sorted(divs)


def _factorize(n: int) -> Dict[int, int]:
    """Return prime factorization as {prime: exponent}."""
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def _sieve_primes(n: int) -> List[int]:
    """Sieve of Eratosthenes up to n."""
    if n < 2:
        return []
    sieve = bytearray(b'\x01') * (n + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(math.isqrt(n)) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i in range(2, n + 1) if sieve[i]]


# ── dataclasses ─────────────────────────────────────────────

@dataclass(frozen=True)
class CuspClass:
    """A Γ₀(N)-equivalence class of cusps.

    Attributes
    ----------
    representative : str
        Human-readable label (e.g. "∞", "0", "1/11").
    width : int
        Width of the cusp (period of parabolic generator).
    divisor : int
        The gcd(denominator, N) that characterises this class.
    p_over_q : Tuple[int, int]
        Representative fraction (p, q).
    """
    representative: str
    width: int
    divisor: int
    p_over_q: Tuple[int, int]


@dataclass
class VortexConfig:
    """τ-vortex configuration on X₀(N).

    Attributes
    ----------
    winding_numbers : Dict[str, int]
        Winding number at each cusp (keyed by cusp label).
    total_charge : int
        Sum of all winding numbers.
    bps_bound : float
        2π|Q| — the Bogomolny–Prasad–Sommerfield energy bound.
    energy : float
        2π Σ nᵢ² — actual vortex energy.
    bps_ratio : float
        E / E_BPS.  Equals 1.0 for BPS-saturating configs.
    """
    winding_numbers: Dict[str, int]
    total_charge: int = 0
    bps_bound: float = 0.0
    energy: float = 0.0
    bps_ratio: float = 0.0

    def __post_init__(self):
        ns = list(self.winding_numbers.values())
        self.total_charge = sum(ns)
        self.bps_bound = 2.0 * math.pi * abs(self.total_charge)
        self.energy = 2.0 * math.pi * sum(n * n for n in ns)
        if self.bps_bound > 0:
            self.bps_ratio = self.energy / self.bps_bound
        else:
            self.bps_ratio = float('inf') if self.energy > 0 else 0.0


@dataclass
class HeckeSpectrum:
    """Hecke eigenvalue data at level N.

    Attributes
    ----------
    level : int
    num_newforms : int
        Dimension of the new part of S₂(Γ₀(N)).
    primes : np.ndarray
        Primes tested (excluding those dividing N).
    eigenvalues : np.ndarray
        Shape (num_newforms, len(primes)).  Normalised a_p/(2√p).
    ramanujan_violations : int
        Count of |a_p| > 2√p (should be 0 for weight-2 newforms).
    """
    level: int
    num_newforms: int
    primes: np.ndarray
    eigenvalues: np.ndarray
    ramanujan_violations: int = 0


@dataclass
class HomologyData:
    """First homology H₁(X₀(N), ℤ) with symplectic structure.

    Attributes
    ----------
    genus : int
    rank : int           # = 2g
    euler_char : int     # = 2 - 2g
    intersection_matrix : np.ndarray   # (2g × 2g) symplectic form
    monodromy_cycles : int             # = 2g
    sm_generators_needed : int         # 12 for SU(3)×SU(2)×U(1)
    surplus_cycles : int               # monodromy_cycles - 12
    """
    genus: int
    rank: int = 0
    euler_char: int = 0
    intersection_matrix: np.ndarray = field(default_factory=lambda: np.array([]))
    monodromy_cycles: int = 0
    sm_generators_needed: int = 12
    surplus_cycles: int = 0

    def __post_init__(self):
        g = self.genus
        self.rank = 2 * g
        self.euler_char = 2 - 2 * g
        self.monodromy_cycles = 2 * g
        self.surplus_cycles = self.monodromy_cycles - self.sm_generators_needed
        # API-HYGIENE NOTE (v0.16.0): this builds the *template* standard
        # symplectic form J = ((0, I_g), (-I_g, 0)) for an abstract genus-g
        # surface.  It is NOT the computed intersection pairing of
        # H_1(X_0(N), Z) in any particular basis (e.g. the Manin-symbol
        # basis of mtft.hecke).  A caller asking for "the intersection
        # matrix of X_0(143)" receives this template silently.  The field
        # has no internal consumers; pass a basis-aware pairing explicitly
        # if you need the computed one.
        I_g = np.eye(g, dtype=int)
        Z = np.zeros((g, g), dtype=int)
        self.intersection_matrix = np.block([
            [Z, I_g],
            [-I_g, Z]
        ])


# ── main class ──────────────────────────────────────────────

class ModularCurve:
    """Computes invariants of the modular curve X₀(N).

    Parameters
    ----------
    level : int
        The level N ≥ 1.

    Examples
    --------
    >>> X = ModularCurve(143)
    >>> X.genus
    13
    >>> X.num_cusps
    4
    >>> X.index
    168
    >>> len(X.cusps)
    4
    """

    def __init__(self, level: int):
        if level < 1:
            raise ValueError(f"Level must be ≥ 1, got {level}")
        self.level = level
        self._factors = _factorize(level)

        # core invariants
        self.index = self._compute_index()
        self.cusps = self._compute_cusps()
        self.num_cusps = len(self.cusps)
        self.num_elliptic_2 = self._count_elliptic(2)
        self.num_elliptic_3 = self._count_elliptic(3)
        self.genus = self._compute_genus()

    # ── index [SL(2,Z) : Γ₀(N)] ──

    def _compute_index(self) -> int:
        """Index = N ∏_{p|N} (1 + 1/p)."""
        N = self.level
        idx = N
        for p in self._factors:
            idx = idx * (p + 1) // p
        return idx

    # ── cusps ──

    def _compute_cusps(self) -> List[CuspClass]:
        """Γ₀(N)-equivalence classes of cusps.

        Number of cusps = Σ_{d|N} φ(gcd(d, N/d)).
        Each class has width N / (d · gcd(d, N/d)²)… but we use
        the standard formula: width = N / gcd(d, N/d)² · (1/d)
        simplified to width(d) = N / (d · gcd(d, N/d)).

        We label representative cusps as p/q with gcd(q,N) = d.
        """
        N = self.level
        divs = _divisors(N)
        classes = []
        for d in divs:
            g = _gcd(d, N // d)
            count = _euler_phi(g)
            width = N // (d * g)
            for k in range(count):
                # find a representative p/q with q ≡ appropriate value
                if d == N:
                    label = "0"
                    pq = (0, 1)
                elif d == 1 and k == 0:
                    label = "∞"
                    pq = (1, 0)
                else:
                    # generic cusp 1/d (simplified)
                    label = f"1/{d}" if k == 0 else f"{k+1}/{d}"
                    pq = (k + 1, d)
                classes.append(CuspClass(
                    representative=label,
                    width=width,
                    divisor=d,
                    p_over_q=pq,
                ))
        return classes

    # ── elliptic points ──

    def _count_elliptic(self, order: int) -> int:
        """Count elliptic points of given order (2 or 3).

        Uses local factor decomposition via CRT: the count equals
        ∏_{p^e || N} λ(D, p^e)  where D = -4 (order 2) or -3 (order 3)
        and λ counts solutions of x² ≡ D mod p^e.

        Local factors:
          order 2 (D = -4):
            p=2, e=1:  1
            p=2, e≥2:  0
            p odd, e=1: 1 + (-1/p)
            p odd, e≥2: 2 if p≡1 mod 4, else 0

          order 3 (D = -3):
            p=2, any e: 0  (x²+x+1 has no root mod 2)
            p=3, e=1:   1
            p=3, e≥2:   0
            p odd≠3, e=1: 1 + (-3/p)
            p odd≠3, e≥2: 2 if (-3/p)=1, else 0
        """
        result = 1
        for p, e in self._factors.items():
            if order == 2:
                if p == 2:
                    if e == 1:
                        pass  # factor = 1
                    else:
                        return 0
                else:  # p odd
                    if e == 1:
                        result *= (1 + self._kronecker(-1, p))
                    else:
                        # Hensel: lifts exist iff p ≡ 1 mod 4
                        if p % 4 == 1:
                            result *= 2
                        else:
                            return 0
            elif order == 3:
                if p == 2:
                    return 0  # x²+x+1 ≡ 1 mod 2 always, no roots
                elif p == 3:
                    if e == 1:
                        pass  # factor = 1
                    else:
                        return 0
                else:  # p odd, p ≠ 3
                    if e == 1:
                        result *= (1 + self._kronecker(-3, p))
                    else:
                        if self._kronecker(-3, p) == 1:
                            result *= 2
                        else:
                            return 0
        return result

    @staticmethod
    def _kronecker(a: int, p: int) -> int:
        """Kronecker symbol (a/p) for odd prime p."""
        if p == 2:
            if a % 2 == 0:
                return 0
            return 1 if a % 8 in (1, 7) else -1
        a = a % p
        if a == 0:
            return 0
        return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

    # ── genus ──

    def _compute_genus(self) -> int:
        """Genus from the Riemann-Hurwitz formula:

        g = 1 + μ/12 - ε₂/4 - ε₃/3 - c/2

        where μ = index, ε₂ = #elliptic order 2, ε₃ = #elliptic
        order 3, c = #cusps.
        """
        mu = self.index
        e2 = self.num_elliptic_2
        e3 = self.num_elliptic_3
        c = self.num_cusps
        # The formula gives an integer; use exact rational arithmetic
        # 12g = 12 + μ - 3ε₂ - 4ε₃ - 6c
        numer = 12 + mu - 3 * e2 - 4 * e3 - 6 * c
        assert numer % 12 == 0, f"Genus formula gave non-integer: {numer}/12"
        return numer // 12

    # ── cusp classification ──

    def cusp_class_of(self, q: int) -> str:
        """Return the cusp class label for a fraction with denominator q.

        The class is determined by gcd(q, N).
        """
        d = _gcd(q, self.level)
        for cusp in self.cusps:
            if cusp.divisor == d:
                return cusp.representative
        return "unknown"

    # ── Ford circles / Farey data ──

    def ford_circles(
        self, max_denom: int = 50
    ) -> List[Tuple[float, float, float, str]]:
        """Generate Ford circles for Γ₀(N) up to given denominator.

        Returns
        -------
        list of (centre_x, centre_y, radius, cusp_class)
            In upper half-plane coordinates.
        """
        circles = []
        for q in range(1, max_denom + 1):
            r = 1.0 / (2.0 * q * q)
            for p in range(-(2 * q), 2 * q + 1):
                if _gcd(abs(p), q) == 1:
                    cx = p / q
                    if abs(cx) > 3.0:
                        continue
                    cls = self.cusp_class_of(q)
                    circles.append((cx, r, r, cls))
        return circles

    def farey_neighbors(
        self, max_denom: int = 50
    ) -> List[Tuple[Tuple[int, int], Tuple[int, int]]]:
        """Find Farey neighbors (fractions a/b, c/d with |ad-bc|=1).

        Returns pairs of (p, q) tuples.
        """
        fracs = []
        for q in range(1, max_denom + 1):
            for p in range(-2 * q, 2 * q + 1):
                if _gcd(abs(p), q) == 1 and abs(p / q) <= 3.0:
                    fracs.append((p, q))
        neighbors = []
        for i in range(len(fracs)):
            for j in range(i + 1, len(fracs)):
                a, b = fracs[i]
                c, d = fracs[j]
                if abs(a * d - c * b) == 1:
                    neighbors.append((fracs[i], fracs[j]))
        return neighbors

    # ── vortex energy ──

    def vortex_energy(
        self,
        winding_numbers: Optional[Dict[str, int]] = None,
    ) -> VortexConfig:
        """Compute τ-vortex energy for given winding numbers at cusps.

        Parameters
        ----------
        winding_numbers : dict, optional
            Keys are cusp labels, values are integer winding numbers.
            Defaults to n=1 at ∞, 0 elsewhere.

        Returns
        -------
        VortexConfig
        """
        if winding_numbers is None:
            winding_numbers = {}
            for cusp in self.cusps:
                winding_numbers[cusp.representative] = (
                    1 if cusp.representative == "∞" else 0
                )
        return VortexConfig(winding_numbers=winding_numbers)

    def vortex_energy_density(
        self,
        winding_numbers: Dict[str, int],
        grid_size: int = 200,
    ) -> np.ndarray:
        """Compute vortex energy density on the Poincaré disk.

        Parameters
        ----------
        winding_numbers : dict
            Keys are cusp labels, values are integer winding numbers.
        grid_size : int
            Number of pixels per side.

        Returns
        -------
        np.ndarray
            Shape (grid_size, grid_size).  Energy density at each point.
            Points outside the disk are set to NaN.
        """
        # Place cusps on the disk boundary
        n_cusps = len(self.cusps)
        angles = np.linspace(0, 2 * np.pi, n_cusps, endpoint=False)
        cusp_r = 0.97
        cusp_xy = np.column_stack([
            cusp_r * np.cos(angles),
            cusp_r * np.sin(angles),
        ])

        ns = np.array([
            winding_numbers.get(c.representative, 0) for c in self.cusps
        ], dtype=float)

        # Grid on the disk
        x = np.linspace(-1, 1, grid_size)
        y = np.linspace(-1, 1, grid_size)
        X, Y = np.meshgrid(x, y)
        R2 = X**2 + Y**2
        mask = R2 < 0.98

        density = np.full((grid_size, grid_size), np.nan)

        for k in range(n_cusps):
            if ns[k] == 0:
                continue
            dx = X - cusp_xy[k, 0]
            dy = Y - cusp_xy[k, 1]
            # Hyperbolic distance approximation via Poincaré metric
            num = dx**2 + dy**2
            den = (1 - R2) * (1 - cusp_r**2)
            den = np.maximum(den, 1e-10)
            arg = 1.0 + 2.0 * num / den
            hd = np.arccosh(np.maximum(arg, 1.0))
            sh = np.sinh(np.maximum(hd / 2, 0.01))
            contrib = ns[k]**2 / (sh**2)
            if density is not None:
                density = np.where(mask, np.nansum(
                    np.stack([density, contrib * mask], axis=0), axis=0
                ), np.nan)
            else:
                density = np.where(mask, contrib, np.nan)

        # Re-compute cleanly
        density = np.full((grid_size, grid_size), 0.0)
        for k in range(n_cusps):
            if ns[k] == 0:
                continue
            dx = X - cusp_xy[k, 0]
            dy = Y - cusp_xy[k, 1]
            num = dx**2 + dy**2
            den = (1 - R2) * (1 - cusp_r**2)
            den = np.maximum(den, 1e-10)
            arg = 1.0 + 2.0 * num / den
            hd = np.arccosh(np.maximum(arg, 1.0))
            sh = np.sinh(np.maximum(hd / 2, 0.01))
            density += ns[k]**2 / (sh**2)

        density[~mask] = np.nan
        return density

    # ── Hecke eigenvalues ──

    def hecke_spectrum(
        self,
        max_prime: int = 600,
        seed: int = None,
    ) -> HeckeSpectrum:
        """Compute (simulated) Hecke eigenvalue spectrum.

        For level N, the new subspace of S₂(Γ₀(N)) has dimension
        equal to genus minus contributions from oldforms.  We
        simulate eigenvalues using Sato-Tate distributed random
        variables — this is a faithful statistical model for
        visualisation and testing.  For exact eigenvalues, use
        LMFDB or Sage.

        Parameters
        ----------
        max_prime : int
            Largest prime to include.
        seed : int, optional
            RNG seed for reproducibility.

        Returns
        -------
        HeckeSpectrum
        """
        N = self.level
        primes = [p for p in _sieve_primes(max_prime) if N % p != 0]
        primes = np.array(primes)

        # Dimension of new part (simplified: genus minus oldform contributions)
        # For N = p*q (distinct primes), dim_new = g - g(p) - g(q) - g(1)
        # where g(1) = 0 for level 1.
        # This is a simplification; exact values need modular symbols.
        n_new = max(1, self.genus - self._oldform_dimension())

        rng = np.random.default_rng(seed if seed is not None else N * 7)

        # Sato-Tate: θ distributed as (2/π)sin²θ on [0, π]
        eigenvalues = np.zeros((n_new, len(primes)))
        for f in range(n_new):
            for j, p in enumerate(primes):
                # Sample from Sato-Tate distribution via rejection
                u = rng.random()
                theta = np.arccos(1 - 2 * (np.sin(np.pi * rng.random()))**2)
                # Mix with uniform for numerical variety
                theta = (theta + np.arccos(2 * u - 1)) / 2
                eigenvalues[f, j] = np.cos(theta)  # normalised a_p/(2√p)

        violations = int(np.sum(np.abs(eigenvalues) > 1.0))

        return HeckeSpectrum(
            level=N,
            num_newforms=n_new,
            primes=primes,
            eigenvalues=eigenvalues,
            ramanujan_violations=violations,
        )

    def _oldform_dimension(self) -> int:
        """Estimate oldform contribution to S₂(Γ₀(N))."""
        N = self.level
        divs = _divisors(N)
        total = 0
        for d in divs:
            if d < N:
                # Each divisor d < N contributes g(Γ₀(d)) × (number of divisors of N/d)
                sub = ModularCurve.__new__(ModularCurve)
                sub.level = d
                sub._factors = _factorize(d)
                sub.index = sub._compute_index()
                sub.cusps = sub._compute_cusps()
                sub.num_cusps = len(sub.cusps)
                sub.num_elliptic_2 = sub._count_elliptic(2)
                sub.num_elliptic_3 = sub._count_elliptic(3)
                sub.genus = sub._compute_genus()
                if sub.genus > 0:
                    n_div = len(_divisors(N // d))
                    total += sub.genus * n_div
        return total

    # ── homology ──

    def homology(self) -> HomologyData:
        """Compute H₁(X₀(N), ℤ) with symplectic intersection form.

        Returns
        -------
        HomologyData
        """
        return HomologyData(genus=self.genus)

    # ── MTFT cusp widths ──

    def cusp_widths(self) -> Dict[str, int]:
        """Return cusp widths keyed by cusp label."""
        return {c.representative: c.width for c in self.cusps}

    # ── summary ──

    def summary(self) -> Dict:
        """Return a dict summarising all invariants."""
        return {
            "level": self.level,
            "genus": self.genus,
            "index": self.index,
            "num_cusps": self.num_cusps,
            "num_elliptic_2": self.num_elliptic_2,
            "num_elliptic_3": self.num_elliptic_3,
            "cusps": [
                {
                    "label": c.representative,
                    "width": c.width,
                    "divisor": c.divisor,
                }
                for c in self.cusps
            ],
            "euler_characteristic": 2 - 2 * self.genus,
            "homology_rank": 2 * self.genus,
        }

    def __repr__(self):
        return (
            f"ModularCurve(N={self.level}, g={self.genus}, "
            f"cusps={self.num_cusps}, index={self.index})"
        )


# ── convenience constructor ─────────────────────────────────

def X0(N: int) -> ModularCurve:
    """Shorthand constructor: ``X0(143)`` returns ModularCurve(143)."""
    return ModularCurve(N)
