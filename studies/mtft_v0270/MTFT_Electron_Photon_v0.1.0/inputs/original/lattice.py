"""
MTFT Lattice Gauge Theory
==========================

Implements the MTFT lattice action (Papers 5, 9, 13, 24):

    S_MTFT[U; y] = S_Wilson[U] + S_arith[U; y]

where:
    S_Wilson = β Σ_{x,μ<ν} (1 − (1/N) Re Tr U_{μν}(x))
    S_arith  = κ Σ_{x⃗} Σ_{n≥1} wₙ e^{−2πyn} (1 − (1/N) Re Tr P(x⃗)ⁿ)

The arithmetic term replaces the single Wilson coupling with an infinite
tower of arithmetically determined couplings {wₙ/g²}.

Key features:
    - Reflection positivity (OS2) proven via Osterwalder-Seiler (Paper 24)
    - μ_N(y) > 0 unconditionally → nonperturbative mass gap
    - Continuum limit via asymptotic freedom scaling
    - SU(N) for any N; specialises to SU(3) for QCD

This module provides data structures and Monte Carlo tools for
lattice simulation, NOT GPU-scale computation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from mtft.arithmetic import weight_array, damped_weights, mass_gap_stiffness
from mtft.constants import PI, FEIGENBAUM_DELTA


# ═══════════════════════════════════════════════════════════════
#  SU(N) group elements
# ═══════════════════════════════════════════════════════════════

def random_su_n(N: int, rng: np.random.Generator | None = None) -> np.ndarray:
    """Sample a random SU(N) matrix via QR decomposition of Gaussian."""
    if rng is None:
        rng = np.random.default_rng()
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / math.sqrt(2)
    Q, R = np.linalg.qr(Z)
    # Fix phases to make det = 1
    d = np.diag(R)
    ph = d / np.abs(d)
    Q = Q @ np.diag(ph.conj())
    det = np.linalg.det(Q)
    Q[0, :] /= det ** (1.0 / N)
    return Q


def su_n_identity(N: int) -> np.ndarray:
    return np.eye(N, dtype=complex)


def su_n_center(N: int, k: int = 1) -> np.ndarray:
    """Z_N center element: ω^k · I where ω = e^{2πi/N}."""
    omega = np.exp(2j * PI * k / N)
    return omega * np.eye(N, dtype=complex)


def su_n_near_identity(N: int, epsilon: float = 0.1,
                       rng: np.random.Generator | None = None) -> np.ndarray:
    """SU(N) matrix near the identity (for HMC proposals)."""
    if rng is None:
        rng = np.random.default_rng()
    # Hermitian generator
    H = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    H = epsilon * (H + H.conj().T) / 2.0
    H -= np.eye(N) * np.trace(H) / N
    U = np.eye(N, dtype=complex)
    # Cayley approximation
    U = np.linalg.solve(np.eye(N) - 0.5j * H, np.eye(N) + 0.5j * H)
    U /= np.linalg.det(U) ** (1.0 / N)
    return U


# ═══════════════════════════════════════════════════════════════
#  Lattice configuration
# ═══════════════════════════════════════════════════════════════

@dataclass
class LatticeConfig:
    """
    4D hypercubic lattice gauge configuration for SU(N).

    The lattice is L⁴ (or L³ × L_t for asymmetric).
    Link variables U_μ(x) ∈ SU(N) live on oriented links.

    Parameters
    ----------
    N : int
        Gauge group rank (SU(N)).
    L : int
        Spatial lattice extent.
    L_t : int or None
        Temporal extent (defaults to L).
    """

    N: int = 3
    L: int = 4
    L_t: Optional[int] = None

    # Internal storage
    links: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        if self.L_t is None:
            self.L_t = self.L
        shape = (self.L_t, self.L, self.L, self.L, 4, self.N, self.N)
        # Initialise to identity (cold start)
        self.links = np.zeros(shape, dtype=complex)
        for t in range(self.L_t):
            for x in range(self.L):
                for y_idx in range(self.L):
                    for z in range(self.L):
                        for mu in range(4):
                            self.links[t, x, y_idx, z, mu] = su_n_identity(self.N)

    @property
    def volume(self) -> int:
        return self.L_t * self.L ** 3

    def hot_start(self, rng: np.random.Generator | None = None):
        """Initialise all links to random SU(N) (hot start)."""
        if rng is None:
            rng = np.random.default_rng()
        for idx in np.ndindex(self.L_t, self.L, self.L, self.L, 4):
            self.links[idx] = random_su_n(self.N, rng)

    def get_link(self, site: tuple, mu: int) -> np.ndarray:
        """U_μ(x) at site = (t, x, y, z)."""
        return self.links[site[0], site[1], site[2], site[3], mu]

    def set_link(self, site: tuple, mu: int, U: np.ndarray):
        self.links[site[0], site[1], site[2], site[3], mu] = U

    def _shift(self, site: tuple, mu: int, direction: int = 1) -> tuple:
        """Shift site by ±1 in direction μ with periodic BC."""
        dims = [self.L_t, self.L, self.L, self.L]
        s = list(site)
        s[mu] = (s[mu] + direction) % dims[mu]
        return tuple(s)


# ═══════════════════════════════════════════════════════════════
#  Plaquette and Polyakov loop
# ═══════════════════════════════════════════════════════════════

def plaquette(cfg: LatticeConfig, site: tuple, mu: int, nu: int) -> np.ndarray:
    """
    Plaquette U_{μν}(x) = U_μ(x) U_ν(x+μ̂) U_μ†(x+ν̂) U_ν†(x).
    """
    U1 = cfg.get_link(site, mu)
    U2 = cfg.get_link(cfg._shift(site, mu), nu)
    U3 = cfg.get_link(cfg._shift(site, nu), mu).conj().T
    U4 = cfg.get_link(site, nu).conj().T
    return U1 @ U2 @ U3 @ U4


def polyakov_loop(cfg: LatticeConfig, spatial_site: tuple) -> np.ndarray:
    """
    Polyakov loop at spatial site (x, y, z):
        P(x⃗) = Π_{t=0}^{L_t−1} U_0(t, x⃗)
    """
    P = su_n_identity(cfg.N)
    for t in range(cfg.L_t):
        P = P @ cfg.get_link((t,) + spatial_site, 0)
    return P


def avg_plaquette(cfg: LatticeConfig) -> float:
    """Average plaquette: ⟨(1/N) Re Tr U_□⟩ over all plaquettes."""
    total = 0.0
    count = 0
    for t in range(cfg.L_t):
        for x in range(cfg.L):
            for y_idx in range(cfg.L):
                for z in range(cfg.L):
                    site = (t, x, y_idx, z)
                    for mu in range(4):
                        for nu in range(mu + 1, 4):
                            P = plaquette(cfg, site, mu, nu)
                            total += np.real(np.trace(P)) / cfg.N
                            count += 1
    return total / count


def avg_polyakov(cfg: LatticeConfig) -> complex:
    """Average Polyakov loop: ⟨(1/N) Tr P(x⃗)⟩ over spatial sites."""
    total = 0.0 + 0j
    count = 0
    for x in range(cfg.L):
        for y_idx in range(cfg.L):
            for z in range(cfg.L):
                P = polyakov_loop(cfg, (x, y_idx, z))
                total += np.trace(P) / cfg.N
                count += 1
    return total / count


# ═══════════════════════════════════════════════════════════════
#  MTFT Lattice Action
# ═══════════════════════════════════════════════════════════════

@dataclass
class MTFTAction:
    """
    The MTFT lattice action:

        S = β Σ (1 − (1/N) Re Tr U_□)
          + κ Σ_{x⃗} Σ_n aₙ(y) (1 − (1/N) Re Tr P(x⃗)ⁿ)

    where aₙ(y) = wₙ e^{−2πyn}.

    Parameters
    ----------
    beta : float
        Wilson coupling β = 2N/g².
    kappa : float
        Arithmetic deformation strength.
    y : float
        Modular depth parameter.
    n_max : int
        Truncation of the arithmetic sum.
    """

    beta: float = 6.0
    kappa: float = 1.0
    y: float = 0.18174      # default: confinement depth
    n_max: int = 50

    def __post_init__(self):
        self._an = None

    @property
    def an(self) -> np.ndarray:
        """Cached damped arithmetic weights."""
        if self._an is None:
            self._an = damped_weights(self.y, self.n_max)
        return self._an

    def wilson_action(self, cfg: LatticeConfig) -> float:
        """S_Wilson = β Σ (1 − (1/N) Re Tr U_□)."""
        total = 0.0
        for t in range(cfg.L_t):
            for x in range(cfg.L):
                for y_idx in range(cfg.L):
                    for z in range(cfg.L):
                        site = (t, x, y_idx, z)
                        for mu in range(4):
                            for nu in range(mu + 1, 4):
                                P = plaquette(cfg, site, mu, nu)
                                total += 1.0 - np.real(np.trace(P)) / cfg.N
        return self.beta * total

    def arithmetic_action(self, cfg: LatticeConfig) -> float:
        """
        S_arith = κ Σ_{x⃗} Σ_n aₙ(y) (1 − (1/N) Re Tr P(x⃗)ⁿ).

        The Polyakov deformation with arithmetic weights.
        """
        an = self.an
        total = 0.0
        for x in range(cfg.L):
            for y_idx in range(cfg.L):
                for z in range(cfg.L):
                    P = polyakov_loop(cfg, (x, y_idx, z))
                    Pn = su_n_identity(cfg.N)
                    for n_idx in range(self.n_max):
                        n = n_idx + 1
                        Pn = Pn @ P
                        tr_Pn = np.real(np.trace(Pn)) / cfg.N
                        total += an[n_idx] * (1.0 - tr_Pn)
        return self.kappa * total

    def total_action(self, cfg: LatticeConfig) -> float:
        """S_MTFT = S_Wilson + S_arith."""
        return self.wilson_action(cfg) + self.arithmetic_action(cfg)

    def mass_gap_bound(self, N: int = 3) -> float:
        """
        Analytic lower bound on the mass gap from arithmetic stiffness:

            m²_gap ≥ κ·μ_N(y) / (K·λ(α))

        Returns μ_N(y) (the arithmetic part).
        """
        return mass_gap_stiffness(self.y, N=N, n_max=self.n_max)


# ═══════════════════════════════════════════════════════════════
#  Single-link Metropolis update
# ═══════════════════════════════════════════════════════════════

def _staple_sum(cfg: LatticeConfig, site: tuple, mu: int,
                beta: float) -> np.ndarray:
    """
    Compute the sum of staples for link U_μ(x).

    The staple for plaquette (μ,ν) is:
        U_ν(x+μ̂) U_μ†(x+ν̂) U_ν†(x)  +  U_ν†(x+μ̂−ν̂) U_μ†(x−ν̂) U_ν(x−ν̂)
    """
    N = cfg.N
    S = np.zeros((N, N), dtype=complex)
    for nu in range(4):
        if nu == mu:
            continue
        # Forward staple
        site_mu = cfg._shift(site, mu)
        site_nu = cfg._shift(site, nu)
        S += (cfg.get_link(site_mu, nu)
              @ cfg.get_link(site_nu, mu).conj().T
              @ cfg.get_link(site, nu).conj().T)
        # Backward staple
        site_mu_nub = cfg._shift(cfg._shift(site, mu), nu, -1)
        site_nub = cfg._shift(site, nu, -1)
        S += (cfg.get_link(site_mu_nub, nu).conj().T
              @ cfg.get_link(site_nub, mu).conj().T
              @ cfg.get_link(site_nub, nu))
    return S


def metropolis_sweep(cfg: LatticeConfig, action: MTFTAction,
                     n_hits: int = 10, epsilon: float = 0.2,
                     rng: np.random.Generator | None = None) -> float:
    """
    One Metropolis sweep over all links.

    For simplicity, only the Wilson term drives the local update
    (the Polyakov term is nonlocal and enters via reweighting or
    global updates in production code).

    Returns acceptance rate.
    """
    if rng is None:
        rng = np.random.default_rng()

    accepted = 0
    total = 0
    N = cfg.N

    for t in range(cfg.L_t):
        for x in range(cfg.L):
            for y_idx in range(cfg.L):
                for z in range(cfg.L):
                    site = (t, x, y_idx, z)
                    for mu in range(4):
                        staple = _staple_sum(cfg, site, mu, action.beta)
                        U_old = cfg.get_link(site, mu)
                        for _ in range(n_hits):
                            R = su_n_near_identity(N, epsilon, rng)
                            U_new = R @ U_old
                            # ΔS = −(β/N) Re Tr((U_new − U_old) · staple†)
                            dS = -(action.beta / N) * np.real(
                                np.trace((U_new - U_old) @ staple.conj().T)
                            )
                            if dS < 0 or rng.random() < math.exp(-dS):
                                U_old = U_new
                                accepted += 1
                            total += 1
                        cfg.set_link(site, mu, U_old)

    return accepted / max(total, 1)


# ═══════════════════════════════════════════════════════════════
#  Observables
# ═══════════════════════════════════════════════════════════════

def wilson_loop(cfg: LatticeConfig, site: tuple,
                R: int, T: int, mu: int = 1, nu: int = 0) -> np.ndarray:
    """
    R×T Wilson loop in the (μ, ν) plane starting at site.

    For extracting the static potential V(R) from ⟨W(R,T)⟩ ~ e^{−V(R)T}.
    """
    N = cfg.N
    # Bottom: R links in μ direction
    W = su_n_identity(N)
    s = site
    for _ in range(R):
        W = W @ cfg.get_link(s, mu)
        s = cfg._shift(s, mu)
    # Right: T links in ν direction
    for _ in range(T):
        W = W @ cfg.get_link(s, nu)
        s = cfg._shift(s, nu)
    # Top: R links in −μ direction
    for _ in range(R):
        s = cfg._shift(s, mu, -1)
        W = W @ cfg.get_link(s, mu).conj().T
    # Left: T links in −ν direction
    for _ in range(T):
        s = cfg._shift(s, nu, -1)
        W = W @ cfg.get_link(s, nu).conj().T
    return W


def creutz_ratio(cfg: LatticeConfig, R: int, T: int) -> float:
    """
    Creutz ratio χ(R,T) = −ln(⟨W(R,T)⟩⟨W(R−1,T−1)⟩ / (⟨W(R,T−1)⟩⟨W(R−1,T)⟩)).

    Converges to the string tension σ as R, T → ∞.
    """
    def avg_W(r, t):
        total = 0.0
        count = 0
        for site in [(0, 0, 0, 0)]:  # simplified: single origin
            W = wilson_loop(cfg, site, r, t)
            total += np.real(np.trace(W)) / cfg.N
            count += 1
        return total / count

    w_rt = avg_W(R, T)
    w_r1t1 = avg_W(R - 1, T - 1)
    w_rt1 = avg_W(R, T - 1)
    w_r1t = avg_W(R - 1, T)

    numer = w_rt * w_r1t1
    denom = w_rt1 * w_r1t
    if abs(denom) < 1e-30 or numer / denom <= 0:
        return float("nan")
    return -math.log(numer / denom)
