"""CC-69 (2026-10-17): the character repair of CC-35 adopted — the MTFT plaquette action with the arithmetic tower on Sym^n characters,
and what follows from it: reflection positivity and a strong-coupling mass gap by Osterwalder–Seiler, restored centre symmetry, the size
of the arithmetic layer, and the N-ality (centre-vortex) content of the tower.

The action.  Every plaquette carries
    s(U) = beta (1 - Re Tr U / N) + (kappa / N) sum_{n >= 2} a_n(y) (d_n - Re chi_{Sym^n}(U)),      a_n = w_n e^{-2 pi y n},  w_1 = 0,
with d_n = dim Sym^n (n + 1 on SU(2), (n + 1)(n + 2)/2 on SU(3)).  CC-35 showed that the original placement — power sums Tr P^n on Polyakov
lines — is not covered by the Osterwalder–Seiler argument (power sums are virtual characters) and fails site-reflection positivity at strong
coupling; it named this repair.  Here the repair is the programme's action.

What is a theorem for it (EXACT).
  1. The Boltzmann weight f = e^{-s} is the exponential of a non-negative combination of characters, hence a positive-type class function:
     every coefficient c_rho in f = sum_rho d_rho c_rho chi_rho is >= 0 (`coefficients`; checked here to 1e-16 for kappa up to 1000).  The
     Osterwalder–Seiler argument (Ann. Phys. 110 (1978) 440) then applies verbatim: link-reflection positivity, a positive transfer matrix, and
     for small character activity a mass gap and an area law from the cluster expansion.  This is the position of every reflection-positive
     lattice gauge action: rigorous at strong coupling, the continuum limit with a uniform gap open.
  2. Centre symmetry.  Plaquettes are blind to the centre, so the plaquette action is exactly Z_N-symmetric and the Polyakov loop is an order
     parameter.  Both Polyakov-line placements break it explicitly: V(z 1) - V(1) = kappa sum_n a_n (1 - cos 2 pi n / N) per site for power sums
     and kappa sum_n a_n (d_n / N)(1 - cos 2 pi n / N) for Sym^n characters (`centre_cost`).  The original action therefore had no confinement
     order parameter at any coupling.
  3. N-ality.  Sym^n has N-ality n mod N, so the n-th arithmetic weight couples to the N-ality class of its index: a Z_N centre vortex (a 2-sheet
     in four dimensions) linked once with a loop multiplies the Sym^n contribution by e^{2 pi i n / N} (`nality`, `vortex_overlap`).  In four
     dimensions a 1-cycle links a 2-cycle (1 + 2 = 4 - 1), a knotted circle is always trivial, and a 2-sphere can knot — the Wilson loop and the
     vortex sheet are the physical pair.

What is computed (NUMERICAL, Weyl integration; the Sym^n characters from h_n = e_1 h_{n-1} - e_2 h_{n-2} + e_3 h_{n-3}, irreducible characters by
Jacobi–Trudi; orthonormality and the SU(2) Bessel coefficients checked in tests/test_repaired_action.py).
  * u_f = c_f / c_0 and the leading strong-coupling observables  sigma a^2 = -ln u_f,  m(0++) a = -4 ln u_f  (Münster 1981).
  * The tower's size: at kappa = 1 its plaquette couplings are (kappa/N) a_n = 0.012, 0.004, 0.002, ... (SU(3)) and it moves u_f by 0.8% (SU(3),
    beta = 5.7) to 1.8% (SU(2), beta = 1): worth Delta beta = 0.05-0.06 in Wilson units (`beta_equivalent`).  A 10% change at beta = 5.7 needs
    kappa = 9.3; kappa = 100 orders the plaquettes almost completely (u_f = 0.93).  Alone (beta = 0) the tower gives u_f = 7e-5 at kappa = 1,
    scaling as kappa^2: a Sym^n character never contains the fundamental, Sym^2 x Sym^3 does.
  * The lock y_c = 0.18174 plays no role in the repaired action: y only sets the decay of the couplings.
  * Mixed actions with large higher-representation couplings have first-order bulk transitions (Bhanot–Creutz 1981); the location of that line
    for this action is not computed here, so kappa cannot be taken large freely.

Reading.  The gap of the repaired action is Wilson's gap, corrected at the per-cent level by the arithmetic layer; the arithmetic does not
supply a mass gap.  Separately: 'Hosotani Hessian positivity' (Matter-05, PG-6) concerns the masses of compactification moduli, not the
spectral gap of pure Yang–Mills on R^4 — a different question.

Status: items 1-3 EXACT; the numbers NUMERICAL (grids: SU(2) 4000 points, SU(3) 420 x 420; coefficients converged to 1e-10)."""
import math
from functools import lru_cache

import numpy as np

from mtft.arithmetic import weight_array

Y_DEFAULT = 0.18174
N_MAX = 24


def a_n(y=Y_DEFAULT, n_max=N_MAX):
    """a_n = w_n e^{-2 pi y n}, indexed from 0 (a_0 unused, a_1 = 0 since w_1 = 0)."""
    w = weight_array(n_max)
    return np.array([0.0] + [float(w[n - 1]) * math.exp(-2 * math.pi * y * n) for n in range(1, n_max + 1)])


def sym_dim(N, n):
    return n + 1 if N == 2 else (n + 1) * (n + 2) // 2


class WeylGroupIntegrator:
    """Class-function integration on SU(2) (eigenphase, measure (2/pi) sin^2) or SU(3) (two eigenphases, Vandermonde^2), with the complete
    homogeneous polynomials h_n = chi_{Sym^n} and the irreducible characters chi_lambda (SU(2): lambda = 2j; SU(3): partitions (l1, l2))."""

    def __init__(self, N, M=None, h_max=60):
        self.N = N
        M = M or (4000 if N == 2 else 420)
        if N == 2:
            th = (np.arange(M) + 0.5) * math.pi / M
            self.e1 = 2 * np.cos(th); self.mu = (2 / math.pi) * np.sin(th) ** 2 * (math.pi / M)
        else:
            t = (np.arange(M) + 0.37) * 2 * math.pi / M
            T1, T2 = np.meshgrid(t, t, indexing="ij")
            x = np.stack([np.exp(1j * T1), np.exp(1j * T2), np.exp(-1j * (T1 + T2))])
            vdm = abs(x[0] - x[1]) ** 2 * abs(x[0] - x[2]) ** 2 * abs(x[1] - x[2]) ** 2
            self.e1 = x.sum(0).ravel(); self.mu = (vdm / vdm.sum()).ravel()
        h = [np.ones_like(self.e1), self.e1]
        for n in range(2, h_max + 1):
            h.append(self.e1 * h[n - 1] - h[n - 2] if N == 2 else self.e1 * h[n - 1] - np.conj(self.e1) * h[n - 2] + (h[n - 3] if n >= 3 else 0))
        self.h = h

    def chi(self, lam):
        if self.N == 2:
            return self.h[lam]
        l1, l2 = lam
        return self.h[l1] * self.h[l2] - self.h[l1 + 1] * (self.h[l2 - 1] if l2 >= 1 else 0)

    def dim(self, lam):
        return lam + 1 if self.N == 2 else (lam[0] - lam[1] + 1) * (lam[1] + 1) * (lam[0] + 2) // 2

    def integrate(self, f):
        return np.sum(self.mu * f)

    def reps(self, top=8):
        return list(range(top + 1)) if self.N == 2 else [(l1, l2) for l1 in range(top + 1) for l2 in range(l1 + 1)]

    @property
    def trivial(self):
        return 0 if self.N == 2 else (0, 0)

    @property
    def fundamental(self):
        return 1 if self.N == 2 else (1, 0)


@lru_cache(maxsize=4)
def integrator(N, M=None):
    return WeylGroupIntegrator(N, M)


def minus_action(G, beta, kappa, y=Y_DEFAULT, n_max=N_MAX):
    """-s(U) up to a constant on the integrator's grid: beta Re Tr U / N + (kappa / N) sum_{n >= 2} a_n Re chi_{Sym^n}(U)."""
    a = a_n(y, n_max)
    tower = sum(a[n] * np.real(G.h[n]) for n in range(2, n_max + 1))
    return beta * np.real(G.e1) / G.N + kappa * tower / G.N


def coefficients(N, beta, kappa, y=Y_DEFAULT, top=8, M=None):
    """Character coefficients c_rho of the plaquette weight, f = sum_rho d_rho c_rho chi_rho, normalised to max f = 1 (ratios are what matter)."""
    G = integrator(N, M)
    ma = minus_action(G, beta, kappa, y); f = np.exp(ma - ma.max())
    return {lam: float(np.real(G.integrate(f * np.conj(G.chi(lam)))) / G.dim(lam)) for lam in G.reps(top)}


def u_fundamental(N, beta, kappa, y=Y_DEFAULT, M=None):
    c = coefficients(N, beta, kappa, y, top=2, M=M); G = integrator(N, M)
    return c[G.fundamental] / c[G.trivial]


def strong_coupling(N, beta, kappa, y=Y_DEFAULT):
    """Leading strong-coupling string tension and 0++ glueball mass in lattice units, from u_f (meaningful only deep in strong coupling)."""
    u = u_fundamental(N, beta, kappa, y)
    return {"u_f": u, "sigma_a2": -math.log(u) if u > 0 else math.inf, "glueball_a": -4 * math.log(u) if u > 0 else math.inf}


def beta_equivalent(N, beta, kappa, y=Y_DEFAULT, beta_max=400.0):
    """The Wilson coupling with the same u_f as the repaired action at (beta, kappa): the tower's whole effect at strong coupling, in units of beta."""
    from scipy.optimize import brentq
    target = u_fundamental(N, beta, kappa, y)
    if not (u_fundamental(N, 0.0, 0.0, y) + 1e-12 < target < u_fundamental(N, beta_max, 0.0, y)):
        return math.nan
    return brentq(lambda b: u_fundamental(N, b, 0.0, y) - target, 1e-9, beta_max, xtol=1e-10)


def nality(N, lam):
    """N-ality: SU(2) spin lam/2 -> lam mod 2; SU(3) partition (l1, l2) -> (l1 + l2) mod 3.  Sym^n has N-ality n mod N."""
    return lam % 2 if N == 2 else (lam[0] + lam[1]) % 3


def tower_nality_classes(N, y=Y_DEFAULT, n_max=N_MAX):
    """How the arithmetic weight of the tower splits over the N-ality classes n mod N (weights a_n d_n / N)."""
    a = a_n(y, n_max); out = {k: 0.0 for k in range(N)}
    for n in range(2, n_max + 1):
        out[n % N] += a[n] * sym_dim(N, n) / N
    return out


def centre_cost(N, which="power", y=Y_DEFAULT, n_max=N_MAX):
    """V(z 1) - V(1) per site in units of kappa for a tower on Polyakov lines (z the generating centre element): the explicit breaking of
    centre symmetry by the original power-sum placement ('power') or by Sym^n characters on Polyakov lines ('sym').  The plaquette placement
    is invariant: plaquettes do not see the centre."""
    a = a_n(y, n_max); cost = 0.0
    for n in range(2, n_max + 1):
        c = 1 - math.cos(2 * math.pi * n / N)
        cost += a[n] * c if which == "power" else a[n] * sym_dim(N, n) / N * c
    return cost


def vortex_overlap(N, beta, kappa, y=Y_DEFAULT, top=8):
    """DIAGNOSTIC: the overlap of the plaquette weight with its centre-twisted copy, <f(zU) f(U)> / <f(U)^2> = sum d^2 c^2 z^{nality} / sum d^2 c^2.
    1 when the twist is invisible (N-ality-0 weight only: vortices cost nothing, the strong-coupling confining regime); it falls as the
    N-ality-charged representations take over at weak coupling."""
    c = coefficients(N, beta, kappa, y, top=top); G = integrator(N)
    num = sum(G.dim(l) ** 2 * c[l] ** 2 * math.cos(2 * math.pi * nality(N, l) / N) for l in c)
    den = sum(G.dim(l) ** 2 * c[l] ** 2 for l in c)
    return num / den


def activity(N, beta, kappa, y=Y_DEFAULT, top=8):
    """The character activity sum_{rho != 1} d_rho^2 u_rho: the expansion parameter of the strong-coupling cluster expansion."""
    c = coefficients(N, beta, kappa, y, top=top); G = integrator(N); c0 = c[G.trivial]
    return sum(G.dim(l) ** 2 * c[l] / c0 for l in c if l != G.trivial)


def report(y=Y_DEFAULT):
    rows = {}
    for N, betas in ((2, (1.0, 2.0, 2.3)), (3, (2.0, 5.0, 5.7))):
        for b in betas:
            u0, u1 = u_fundamental(N, b, 0.0, y), u_fundamental(N, b, 1.0, y)
            rows[(N, b)] = {"u_wilson": u0, "u_kappa1": u1, "shift": u1 / u0 - 1, "beta_equivalent_kappa1": beta_equivalent(N, b, 1.0, y),
                            **{k: v for k, v in strong_coupling(N, b, 1.0, y).items() if k != "u_f"}, "activity_kappa1": activity(N, b, 1.0, y)}
    mins = {N: min(min(coefficients(N, b, k, y).values()) for b in (0.0, 2.5, 6.0) for k in (0.0, 1.0, 100.0, 1000.0)) for N in (2, 3)}
    return {"min_coefficient": mins, "rows": rows, "centre_cost_per_kappa": {N: {w: centre_cost(N, w, y) for w in ("power", "sym")} for N in (2, 3)},
            "tower_nality_classes": {N: tower_nality_classes(N, y) for N in (2, 3)},
            "vortex_overlap_su3": {b: vortex_overlap(3, b, 1.0, y) for b in (0.5, 2.0, 5.7, 8.0)}}
