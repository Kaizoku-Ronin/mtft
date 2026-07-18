"""
The 3×3 Jacobian Stiffness Engine on J₀(143)   (Paper 30)
==========================================================

Port of the ``mtft_core.MTFT`` computational engine (March 2026,
consolidated from Papers 1–30) into the mtft package, re-based onto the
independently verified X₀(143) data shipped with v0.6.1.

What it computes
----------------
The lepton-sector stiffness quadratic form on the 3-dimensional space of
Galois orbits (f₁, f₂, f₃) of S₂ⁿᵉʷ(Γ₀(143)):

    M_ij(y) = Σ_n n² e^{−2πyn} (1 − cos 2πn/N) · Tr_i(a_n) Tr_j(a_n)

together with its eigen-decomposition, the sector projections against
the holonomy weights w_n, the skeleton (prime) vs Lambertization
(composite) budget, and the Feigenbaum-ratio diagnostics of Paper 30.

Data provenance (v0.7.0 change)
-------------------------------
The per-orbit traces Tr_i(a_n) are taken directly from the verified
tables ``ORBIT_TRACE_F1/F2/F3`` (n ≤ 50) of ``mtft.x0_143`` — the data
certified by the July 2026 audit — instead of the pre-audit
reconstruction used by ``mtft_core`` (exact values at six primes, Hecke
prime-power recurrences, and a proportional 40/60 split for everything
else, which disagreed with the verified tables at 37 of the first 50
entries).

Truncating the trace sum at n = 50 is exact to ~1e−20 (relative) at the
physical depths y ≈ 0.18 because of the Boltzmann factor e^{−2πyn}; the
truncation error stays below ~1e−6 relative down to y ≈ 0.08.  A
RuntimeWarning is emitted for y < 0.08.

Effect on Paper 30's quoted numbers (y = 0.1812, N = 3):

                    pre-audit (Paper 30)    verified data (this module)
    λ₁              1.52628                 1.531802
    λ₂              6.42705                 6.591154
    λ₃              45.46824                45.466583
    λ₃/λ₂           7.0745                  6.8981   (δ_x^meas = 7.2422)
    λ₂/λ₁           4.2109                  4.3029   (δ_F = 4.6692)

Under the verified data the λ₃/λ₂ ↔ δ_x^meas agreement moves from 2.3%
to 4.8%, while λ₂/λ₁ ↔ δ_F improves from 9.8% to 7.9%.  The couplings
(e-μ ≈ 0.02, μ-τ ≈ 0.72) and the eigenvalue hierarchy are unchanged.

At the canonical confinement depth y_conf = 0.18174:

    λ = (1.515739, 6.507311, 45.102141)

Usage
-----
    from mtft.jacobian import JacobianStiffness

    eng = JacobianStiffness()
    M, evals, evecs = eng.jacobian_matrix(0.18174)
    eng.summary()
"""
from __future__ import annotations

import math
import warnings
from typing import Dict, Tuple

import numpy as np

from mtft.arithmetic import weight_array
from mtft.constants import CriticalDepths, DELTA_X_MEASURED, FEIGENBAUM_DELTA
from mtft.x0_143 import ORBIT_TRACE_F1, ORBIT_TRACE_F2, ORBIT_TRACE_F3

N_VERIFIED = 50           # per-orbit traces are audit-verified for n = 1..50
MIN_RELIABLE_DEPTH = 0.08  # below this the n ≤ 50 truncation is no longer safe

SECTORS = ("electron", "muon", "tau")


def _primes_upto(limit: int) -> list:
    """Eratosthenes sieve."""
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(limit) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i in range(2, limit + 1) if sieve[i]]


class JacobianStiffness:
    """
    Core Jacobian stiffness computation on J₀(143) (Paper 30 engine).

    Parameters
    ----------
    n_max : int
        Cutoff for the weight sums used by :meth:`stiffness`.  The
        trace-weighted quantities always use the verified range
        n ≤ 50 (see module docstring).
    """

    def __init__(self, n_max: int = 2000):
        if n_max < N_VERIFIED:
            raise ValueError(f"n_max must be at least {N_VERIFIED}")
        self.n_max = n_max

        # Holonomy weights w_n, n = 1..n_max  (index 0 ↔ n = 1),
        # sourced from mtft.arithmetic — the package's single weight
        # implementation.
        self._w_full = weight_array(n_max)
        self._w_skel = np.zeros(n_max)
        for p in _primes_upto(n_max):
            self._w_skel[p - 1] = math.log(p) / p
        self._w_lamb = self._w_full - self._w_skel

        # Verified per-orbit traces (n = 1..50) — single source of truth.
        self._traces = np.array(
            [ORBIT_TRACE_F1, ORBIT_TRACE_F2, ORBIT_TRACE_F3], dtype=float
        )
        self._n50 = np.arange(1, N_VERIFIED + 1)

    # ── internals ────────────────────────────────────────────────

    def _weights(self, kind: str) -> np.ndarray:
        try:
            return {"full": self._w_full,
                    "skeleton": self._w_skel,
                    "lambertization": self._w_lamb}[kind]
        except KeyError:
            raise ValueError(
                f"weight must be 'full', 'skeleton' or 'lambertization', got {kind!r}"
            ) from None

    def _base50(self, y: float, N: int) -> np.ndarray:
        """n² e^{−2πyn} (1 − cos 2πn/N) over the verified range n ≤ 50."""
        if y < MIN_RELIABLE_DEPTH:
            warnings.warn(
                f"JacobianStiffness truncates the trace sum at n = {N_VERIFIED}; "
                f"for y < {MIN_RELIABLE_DEPTH} the truncation error is no longer "
                "negligible.",
                RuntimeWarning,
                stacklevel=3,
            )
        n = self._n50
        return n**2 * np.exp(-2 * np.pi * y * n) * (1 - np.cos(2 * np.pi * n / N))

    # ── public API (Paper 30) ────────────────────────────────────

    def stiffness(self, y: float, N: int = 3, weight: str = "full") -> float:
        """μ_N(y) = Σ_{n≤n_max} n² w_n e^{−2πyn} (1 − cos 2πn/N)."""
        n = np.arange(1, self.n_max + 1)
        w = self._weights(weight)
        return float(
            np.sum(n**2 * w * np.exp(-2 * np.pi * y * n) * (1 - np.cos(2 * np.pi * n / N)))
        )

    def jacobian_matrix(
        self, y: float, N: int = 3
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        The 3×3 Jacobian stiffness matrix M_ij(y) on J₀(143).

        Returns ``(M, eigenvalues, eigenvectors)`` with eigenvalues
        ascending (``np.linalg.eigh`` convention).
        """
        base = self._base50(y, N)
        tr = self._traces
        M = np.empty((3, 3))
        for i in range(3):
            for j in range(i, 3):
                v = float(np.sum(base * tr[i] * tr[j]))
                M[i, j] = v
                M[j, i] = v
        evals, evecs = np.linalg.eigh(M)
        return M, evals, evecs

    def sector_projection(
        self, y: float, N: int = 3, weight: str = "full"
    ) -> np.ndarray:
        """Projection vector c_i = ⟨w, f_i⟩_y for (electron, muon, tau)."""
        base = self._base50(y, N)
        w = self._weights(weight)[:N_VERIFIED]
        return np.array([float(np.sum(base * w * t)) for t in self._traces])

    def skeleton_budget(self, y: float, N: int = 3) -> Dict[str, Tuple[float, float]]:
        """
        Per-sector split of the stiffness projection into skeleton
        (prime) and Lambertization (composite) contributions.
        """
        c_bulk = self.sector_projection(y, N, "full")
        c_skel = self.sector_projection(y, N, "skeleton")
        result: Dict[str, Tuple[float, float]] = {}
        for i, name in enumerate(SECTORS):
            if abs(c_bulk[i]) > 1e-15:
                sf = c_skel[i] / c_bulk[i]
                result[name] = (sf, 1 - sf)
            else:
                result[name] = (0.0, 0.0)
        return result

    def lambertization_angle(self, y: float, N: int = 3) -> float:
        """Angle (degrees) between the bulk and skeleton projection directions."""
        cb = self.sector_projection(y, N, "full")
        cs = self.sector_projection(y, N, "skeleton")
        nb, ns = np.linalg.norm(cb), np.linalg.norm(cs)
        if nb < 1e-15 or ns < 1e-15:
            return float("nan")
        cosang = float(np.clip(np.dot(cb / nb, cs / ns), -1.0, 1.0))
        return math.degrees(math.acos(cosang))

    def feigenbaum_ratios(self, y: float, N: int = 3) -> Dict[str, float]:
        """Eigenvalue ratios compared to δ_F and the measured δ_x."""
        _, evals, _ = self.jacobian_matrix(y, N)
        r32 = evals[2] / evals[1] if evals[1] > 1e-10 else float("nan")
        r21 = evals[1] / evals[0] if evals[0] > 1e-10 else float("nan")
        return {
            "lambda_3/lambda_2": r32,
            "delta_x_measured": DELTA_X_MEASURED,
            "discrepancy_x": (
                abs(r32 - DELTA_X_MEASURED) / DELTA_X_MEASURED
                if not math.isnan(r32) else float("nan")
            ),
            "lambda_2/lambda_1": r21,
            "delta_F": FEIGENBAUM_DELTA,
            "discrepancy_F": (
                abs(r21 - FEIGENBAUM_DELTA) / FEIGENBAUM_DELTA
                if not math.isnan(r21) else float("nan")
            ),
        }

    def couplings(self, y: float, N: int = 3) -> Dict[str, float]:
        """Normalized off-diagonal couplings |M_ij| / √(M_ii M_jj)."""
        M, _, _ = self.jacobian_matrix(y, N)
        result: Dict[str, float] = {}
        for i, j, label in [(0, 1, "e-mu"), (0, 2, "e-tau"), (1, 2, "mu-tau")]:
            d = math.sqrt(abs(M[i, i] * M[j, j]))
            result[label] = abs(M[i, j]) / d if d > 0 else float("nan")
        return result

    # ── convenience ──────────────────────────────────────────────

    def summary(self, y: float = CriticalDepths.y_conf) -> None:
        """Print the full Paper 30 diagnostic block at depth y."""
        print(f"Jacobian stiffness engine at y = {y}  (verified X0(143) data)")
        print(f"  mu_3 (full)     = {self.stiffness(y, 3, 'full'):.6f}")
        print(f"  mu_3 (skeleton) = {self.stiffness(y, 3, 'skeleton'):.6f}")

        M, evals, evecs = self.jacobian_matrix(y)
        print(f"\n  Jacobian eigenvalues: {np.round(evals, 6)}")
        for i in range(3):
            v = evecs[:, i]
            print(f"    v_{i+1}: [{v[0]:+.4f}, {v[1]:+.4f}, {v[2]:+.4f}]")

        c = self.couplings(y)
        print(f"\n  Couplings: e-mu={c['e-mu']:.4f}, e-tau={c['e-tau']:.4f}, "
              f"mu-tau={c['mu-tau']:.4f}")

        f = self.feigenbaum_ratios(y)
        print(f"\n  Feigenbaum: l3/l2 = {f['lambda_3/lambda_2']:.4f} "
              f"(delta_x_meas = {f['delta_x_measured']}, "
              f"dev {f['discrepancy_x']*100:.1f}%)")
        print(f"              l2/l1 = {f['lambda_2/lambda_1']:.4f} "
              f"(delta_F = {f['delta_F']:.6f}, dev {f['discrepancy_F']*100:.1f}%)")

        b = self.skeleton_budget(y)
        print("\n  Sector budget (skeleton / Lambertization):")
        for name, (sf, lf) in b.items():
            print(f"    {name:>10s}: {sf*100:.1f}% prime / {lf*100:.1f}% composite")

        print(f"\n  Lambertization angle: {self.lambertization_angle(y):.2f} degrees")


if __name__ == "__main__":
    eng = JacobianStiffness()
    eng.summary()
    print()
    print("Paper 30 reference depth (y = 0.1812):")
    _, ev, _ = eng.jacobian_matrix(0.1812)
    print(f"  eigenvalues = {np.round(ev, 6)}   "
          "(Paper 30 quoted 1.52628, 6.42705, 45.46824 from pre-audit data)")
