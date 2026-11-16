"""mtft.maass.trace — the Selberg trace formula on X₀(N), sector by sector, and what it certifies (MA-01, 2026-11-15).

The spectral side is a sum of h(r_j) over Maass cusp forms; the geometric side is arithmetic: areas,
classes of closed geodesics counted through binary quadratic forms, fixed points, mirrors, cusps and
the continuum.  Comparing the two for a list of forms shows whether the list is complete.

Test functions.  h(r) = exp(−(r − r₀)²/2σ²) + exp(−(r + r₀)²/2σ²), with Fourier transform
g(u) = (2σ/√(2π)) cos(r₀u) exp(−σ²u²/2).  A wider h needs shorter geodesics: σ = 0.5 needs lengths
up to about 17, and the package freezes the weights that far (``_data/trace_weights.npz``, levels
1, 11 and 143, recomputed with PARI and equal to the lab's tables).

The formulas (records MP16, MP17, MP19; N squarefree, Γ₀(N) without elliptic points, n = 2^ω(N)
cusps; W_Q the Atkin–Lehner involutions, R the mirror x ↦ −x; ε a sign at each prime, ε_Q its
product over p | Q).  With the pieces

    area        A = (vol/8π) ∫ h(r) r tanh(πr) dr
    geodesics   G = ½ Σ_t w_t g(l_t),       2 cosh(l/2) = t ≥ 3
    glides      G′ = ½ Σ_t w′_t g(l′_t),    2 sinh(l/2) = t ≥ 1
    cusp mirrors  M = (n/2)[g(0)(log 2 + γ) + 4K − 2J] − (g(0)/4) Σ β      (the mirror's cusp ends)
    for Q > 1:  P_Q = ¼ E_Q ∫ h sech(πr) dr + Σ_s w_s g(l_s),   P′_Q = ½ L_Q g(0) + Σ_s w′_s g(l′_s)

the odd and even cusp forms of the sector ε satisfy

    Odd_ε(h)  = [A + G − G′ + M]/n + Σ_{Q>1} ε_Q [P_Q − P′_Q]/(2n)
    Even_ε(h) = [A + G + G′ − M]/n + Σ_{Q>1} ε_Q [P_Q + P′_Q]/(2n)
                + h(0)/2 − g(0) log 2 − Ψ₁(h) − C_ε(h) − δ_ε h(i/2),
    C_ε = C₁ + Σ_p log p [g(0) + Σ_{k≥1} (−ε_p)^k p^{−k/2} g(k log p)],
    C₁  = Ψ_{½}(h) − g(0) log π − 2 Σ_n (Λ(n)/n) g(2 log n) + h(0)/2,   Ψ_a(h) = (1/2π) ∫ h(r) ψ(a + ir) dr,

δ_ε = 1 for the trivial ε (the constant function).  The η-twist (record MP18) has no continuum:
the sum over all its forms is (vol/4π)∫h r tanh(πr) dr + Σ_t w^χ_t g(l_t) − n g(0) log 2.

What it certifies.  For h ≥ 0 every form adds to the spectral side, so geometry − list ≥ 0 and
equals the sum of h over the forms the list lacks.  A missing form at R would add h(R); so no form
is missing wherever h(R) exceeds the residual (``completeness``).  The bound calibrates itself: a
wider h lets in more of the forms above the list's top, but a missing form just below it also
counts for more.

Class numbers.  The weights use L(1, χ_D), which PARI supplied (conditionally on GRH in general).
``l1_dirichlet`` is a second route with no such condition, Dirichlet's finite sum, and
``weight`` rebuilds a weight from it; the gates compare a sample.
"""
from __future__ import annotations

import math
from functools import cache, lru_cache
from pathlib import Path

import mpmath as mp
import numpy as np

from .domain import is_squarefree, prime_factors

EULER = 0.57721566490153286061
C_ELL, C_REFL = 0.25, 0.5          # a fixed point of order 2: ¼∫h sech; a mirror of length L: ½ L g(0)  (MP16's V-checks)
_DATA = Path(__file__).resolve().parent / "_data"
TRUNCATION_ALLOWED = 1e-8      # the largest tail of a geometric sum the formulas accept from the frozen lengths


# ------------------------------------------------------------------------------------------------ test functions
class TestFunction:
    """h(r) = exp(−(r − r₀)²/2σ²) + exp(−(r + r₀)²/2σ²) and its transform g(u)."""

    __test__ = False                                                   # not a pytest class

    def __init__(self, r0: float = 0.0, sigma: float = 0.5):
        if sigma <= 0 or r0 < 0:
            raise ValueError("sigma must be positive and r0 non-negative")
        self.r0, self.sigma = float(r0), float(sigma)
        self._c = 2 * self.sigma / math.sqrt(2 * math.pi)
        self._integrals = None

    def h(self, r):
        r = np.asarray(r, dtype=float)
        s2 = 2 * self.sigma ** 2
        return np.exp(-(r - self.r0) ** 2 / s2) + np.exp(-(r + self.r0) ** 2 / s2)

    def g(self, u):
        u = np.asarray(u, dtype=float)
        return self._c * np.cos(self.r0 * u) * np.exp(-self.sigma ** 2 * u * u / 2)

    @property
    def h_i_half(self) -> float:
        """h(i/2), real: the constant function's term."""
        s2 = 2 * self.sigma ** 2
        return 2 * math.exp((0.25 - self.r0 ** 2) / s2) * math.cos(self.r0 / s2)

    def length_needed(self, tol: float = 1e-13) -> float:
        """The length beyond which |g| is below tol·g(0): the geodesics a geometric side needs."""
        return math.sqrt(2 * math.log(1 / tol)) / self.sigma

    def truncation(self, cut: float) -> float:
        """An estimate of what the classes longer than ``cut`` add to a geometric sum.

        The total weight of the classes of length near l grows like e^{l/2}, so the tail is about
        e^{cut/2} |g(cut)| / (σ² cut − ½) once σ² cut > ½ (infinite otherwise).
        """
        rate = self.sigma ** 2 * cut - 0.5
        if rate <= 0:
            return math.inf
        return math.exp(cut / 2) * self._c * math.exp(-self.sigma ** 2 * cut * cut / 2) / rate

    def gp(self, u):
        """g′(u)."""
        u = np.asarray(u, dtype=float)
        return self._c * (-self.r0 * np.sin(self.r0 * u) - self.sigma ** 2 * u * np.cos(self.r0 * u)) * np.exp(-self.sigma ** 2 * u * u / 2)

    def integrals(self) -> dict:
        """The one-dimensional integrals every formula uses, by composite Gauss–Legendre rules (computed once).

        ``r tanh`` = ∫h r tanh(πr), ``sech`` = ∫h sech(πr), ``order 3`` = ∫h cosh(πr/3)/cosh(πr), all
        over the real line; ``psi a`` = (1/2π)∫h ψ(a + ir); K = −¼∫₀^∞ g′(u) log(sinh(u/2)/2) du and
        J = −∫₀^∞ g′(u) log cosh(u/4) du (the mirror's cusp ends).
        """
        if self._integrals is not None:
            return self._integrals
        r_hi = self.r0 + 12 * self.sigma
        r, w = _gauss(np.linspace(0.0, r_hi, max(2, math.ceil(r_hi / 0.25)) + 1))
        h = self.h(r)
        u_hi = 12 / self.sigma
        edges = np.concatenate([[0.0, 1e-6, 1e-4, 1e-2, 0.1, 0.5], np.linspace(1.0, u_hi, max(2, math.ceil(u_hi / 0.5)) + 1)])
        u, wu = _gauss(edges)
        gp = self.gp(u)
        self._integrals = {"r tanh": 2 * float(np.sum(w * h * r * np.tanh(np.pi * r))),
                "sech": 2 * float(np.sum(w * h / np.cosh(np.pi * r))),
                "order 3": 2 * float(np.sum(w * h * np.cosh(np.pi * r / 3) / np.cosh(np.pi * r))),
                "psi 1": float(np.sum(w * h * digamma(1.0 + 1j * r).real)) / math.pi,
                "psi 1/2": float(np.sum(w * h * digamma(0.5 + 1j * r).real)) / math.pi,
                "K": -0.25 * float(np.sum(wu * gp * np.log(np.sinh(u / 2) / 2))),
                "J": -float(np.sum(wu * gp * np.log(np.cosh(u / 4)))),
                "h(0)": float(self.h(0.0)), "g(0)": float(self.g(0.0))}
        return self._integrals

    def smooth_tail(self, r_top: float, density: float) -> float:
        """∫_{r_top}^∞ density·r·h(r) dr: a generous count of what the forms above r_top add."""
        hi = max(r_top, self.r0) + 12 * self.sigma
        r, w = _gauss(np.linspace(r_top, hi, 25))
        return float(np.sum(w * density * r * self.h(r)))


_GL16 = np.polynomial.legendre.leggauss(16)


def _gauss(edges):
    """Nodes and weights of the 16-point Gauss–Legendre rule on each panel between consecutive edges."""
    edges = np.asarray(edges, dtype=float)
    a, b = edges[:-1, None], edges[1:, None]
    x = ((b - a) / 2 * _GL16[0][None, :] + (a + b) / 2).ravel()
    w = ((b - a) / 2 * _GL16[1][None, :]).ravel()
    return x, w


def digamma(z):
    """ψ(z) for complex arrays with Re z > 0: eight steps of the recurrence, then the asymptotic series (about 1e-15)."""
    z = np.asarray(z, dtype=complex)
    shift = np.zeros_like(z)
    for k in range(8):
        shift += 1.0 / (z + k)
    w = z + 8
    w2 = 1.0 / (w * w)
    series = w2 * (-1 / 12 + w2 * (1 / 120 + w2 * (-1 / 252 + w2 * (1 / 240 + w2 * (-1 / 132 + w2 * (691 / 32760 - w2 / 12))))))
    return np.log(w) - 0.5 / w + series - shift


# ------------------------------------------------------------------------------------------------ arithmetic, second route
def _prime_discriminants(D: int) -> list:
    """A fundamental discriminant as a product of prime discriminants: −4, ±8, and p* = (−1)^((p−1)/2) p."""
    out, m = [], abs(D)
    for p in prime_factors(m):
        if p == 2:
            continue
        out.append(p if p % 4 == 1 else -p)
    rest = D
    for q in out:
        rest //= q
    if rest not in (1, -4, 8, -8):
        raise ValueError(f"{D} is not a fundamental discriminant")
    if rest != 1:
        out.append(rest)
    return out


def kronecker_table(D: int) -> np.ndarray:
    """χ_D(a) = (D/a) for a = 0 … |D| − 1, D a fundamental discriminant."""
    m = abs(D)
    a = np.arange(m)
    chi = np.ones(m, dtype=np.int64)
    for q in _prime_discriminants(D):
        if q == -4:
            tab = np.array([0, 1, 0, -1])
        elif q == 8:
            tab = np.array([0, 1, 0, -1, 0, -1, 0, 1])
        elif q == -8:
            tab = np.array([0, 1, 0, 1, 0, -1, 0, -1])
        else:
            p = abs(q)
            tab = -np.ones(p, dtype=np.int64)
            tab[(np.arange(1, p) ** 2) % p] = 1
            tab[0] = 0
        chi *= tab[a % len(tab)]
    return chi


def kronecker(D: int, p: int) -> int:
    """(D/p) for a prime p."""
    if p == 2:
        return 0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1)
    r = D % p
    if r == 0:
        return 0
    return 1 if pow(r, (p - 1) // 2, p) == 1 else -1


@cache
def l1_dirichlet(D: int) -> float:
    """L(1, χ_D) for a fundamental discriminant D > 0 by Dirichlet's finite sum: −D^{−½} Σ χ_D(a) log sin(πa/D)."""
    if D <= 1:
        raise ValueError("D must be a fundamental discriminant larger than 1")
    chi = kronecker_table(D)
    a = np.arange(1, D)
    return float(-np.sum(chi[1:] * np.log(np.sin(np.pi * a / D))) / math.sqrt(D))


def core_discriminant(Delta: int):
    """Δ = D c₀² with D fundamental: returns (D, c₀)."""
    if Delta <= 0 or Delta % 4 not in (0, 1):
        raise ValueError("Delta must be a positive discriminant")
    m, k, p, n = 1, 1, 2, Delta
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        k *= p ** (e // 2)
        m *= p ** (e % 2)
        p += 1
    m *= n
    return (m, k) if m % 4 == 1 else (4 * m, k // 2)


def weight(Delta: int, N: int, Q: int, twist_trace=None) -> float:
    """W(Δ) of the coset W_Q at level N, from Dirichlet's sum: the second route for a frozen weight.

    W(Δ) = Σ_{f | c₀} f⁻¹ L(1, χ_D) ∏_{p | c₀/f} (1 − χ_D(p)/p) ∏_{p | N/Q} ρ_p, with Δ = D c₀² and
    ρ_p = p + 1 if p | f, else 1 + χ_D(p) (1 when p | c₀/f).  With ``twist_trace`` = t the η-twisted
    weight: the term of index f changes sign when t is even and f is odd.
    """
    D, c0 = core_discriminant(Delta)
    if D == 1:
        raise ValueError("Delta is a square: no hyperbolic class")
    L = l1_dirichlet(D)
    ps = prime_factors(c0)
    tot = 0.0
    for f in range(1, c0 + 1):
        if c0 % f:
            continue
        c = c0 // f
        term = L / f
        for p in ps:
            if c % p == 0:
                term *= 1 - kronecker(D, p) / p
        for p in prime_factors(N // Q):
            term *= (p + 1) if f % p == 0 else (1 + (0 if c % p == 0 else kronecker(D, p)))
        if twist_trace is not None and twist_trace % 2 == 0 and f % 2 == 1:
            term = -term
        tot += term
    return tot


def class_number_negative(D: int) -> int:
    """The number of reduced primitive positive definite forms of discriminant D < 0."""
    cnt, a = 0, 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (b < 0 and a == c) or math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            cnt += 1
        a += 1
    return cnt


def fixed_points(N: int, Q: int) -> float:
    """E_Q: the fixed points of W_Q on X₀(N), counted with the weights of the trace formula."""
    tot = 0.0
    for f in (1, 2):
        if (4 * Q) % (f * f):
            continue
        Dp = -4 * Q // (f * f)
        if Dp % 4 not in (0, 1):
            continue
        w = {-4: 2, -3: 3}.get(Dp, 1)
        rho = 1
        for p in prime_factors(N // Q):
            rho *= (p + 1) if f % p == 0 else 1 + kronecker(Dp, p)
        tot += class_number_negative(Dp) * rho / w
    return tot


def mirror_length(N: int, Q: int) -> float:
    """L_Q = 2√Q W(4Q): the total length of the mirrors of R W_Q, Q > 1."""
    return 2 * math.sqrt(Q) * weight(4 * Q, N, Q)


def mirror_ends(N: int) -> list:
    """The cusp ends of the mirrors of R on X₀(N): two mirrors leave each cusp (x = 0 and x = ½ in its own coordinate).

    Each entry gives the cusp 1/q with its width, the mirror's other end p₂/q₂ with its width, and
    β = log(w₁w₂) + 2 log|p₂q₁ − p₁q₂|.
    """
    if not is_squarefree(N):
        raise ValueError("implemented for squarefree N")
    ps = prime_factors(N)
    divs = sorted({math.prod(ps[i] for i in range(len(ps)) if k >> i & 1) for k in range(2 ** len(ps))})
    out = []
    for d in divs:
        w, p, q, u = N // d, 1, d, 0
        while (2 * q * (1 + u * q)) % N:                               # A = [[1, u], [q, v]], v = 1 + u q, with N | 2 q v
            u += 1
        v = 1 + u * q
        for j in (0, 1):
            num, den = (u, v) if j == 0 else (p * w + 2 * u, q * w + 2 * v)
            if den < 0:
                num, den = -num, -den
            gg = math.gcd(num, den)
            p2, q2 = (num // gg, den // gg) if gg else (1, 0)
            w2 = N // math.gcd(q2, N) if q2 else 1
            cval = abs(p2 * q - p * q2)
            out.append({"cusp": f"1/{d}" if d < N else "inf", "width": w, "mirror": "x = 0" if j == 0 else "x = 1/2",
                        "other end": f"{p2}/{q2}", "other width": w2, "beta": math.log(w * w2) + 2 * math.log(cval)})
    return out


# ------------------------------------------------------------------------------------------------ the frozen weights
@lru_cache(maxsize=1)
def _weights():
    z = np.load(_DATA / "trace_weights.npz")
    return {k: z[k] for k in z.files}


def frozen_levels() -> tuple:
    return (1, 11, 143)


def classes(N: int, Q: int = 1, kind: str = "hyp", twist: bool = False):
    """Frozen classes of a coset: arrays (s, length, weight).  kind is "hyp" (det +Q) or "gli" (det −Q, the mirror's coset)."""
    W = _weights()
    key = f"N{N}_tw_hyp" if twist else f"N{N}_Q{Q}_{kind}"
    if key + "_s" not in W:
        raise ValueError(f"no frozen weights for level {N}, Q = {Q}, kind {kind!r}, twist = {twist}")
    s = W[key + "_s"].astype(float)
    x = s * math.sqrt(Q) / 2
    length = 2 * (np.arccosh(x) if kind == "hyp" else np.arcsinh(x))
    return W[key + "_s"], length, W[key + "_w"]


def length_bound() -> float:
    """The length to which the weights are frozen (17)."""
    return float(_weights()["lmax"])


# ------------------------------------------------------------------------------------------------ the formulas
def lambda_sum(tf: TestFunction) -> float:
    """Σ_n (Λ(n)/n) g(2 log n) over prime powers n."""
    nmax = int(math.exp(tf.length_needed(1e-18) / 2)) + 1
    sieve = np.ones(nmax + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(nmax ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    tot = 0.0
    for p in np.nonzero(sieve)[0]:
        lp, q = math.log(p), int(p)
        while q <= nmax:
            tot += lp / q * float(tf.g(2 * math.log(q)))
            q *= int(p)
    return tot


def continuum_level1(tf: TestFunction, route: str = "series") -> float:
    """C₁(h) = −(1/4π) ∫ h (φ′/φ)(½ + ir) dr for PSL(2, ℤ), by its explicit formula or by quadrature."""
    I = tf.integrals()
    if route == "series":
        return I["psi 1/2"] - I["g(0)"] * math.log(math.pi) - 2 * lambda_sum(tf) + I["h(0)"] / 2
    if route == "quad":
        def theta(r):                                                  # (φ′/φ)(½ + ir), real on the line
            s = mp.mpc(1, 2 * r)
            return -2 * mp.re(mp.digamma(mp.mpc(0.5, r))) + 2 * mp.log(mp.pi) - 4 * mp.re(mp.zeta(s, derivative=1) / mp.zeta(s))
        # the integrand is smooth at r = 0: ζ′/ζ(1 + 2ir) has real part finite there
        r_hi = tf.r0 + 12 * tf.sigma
        r, w = _gauss(np.linspace(0.0, r_hi, max(2, math.ceil(r_hi / 0.5)) + 1))
        return -float(np.sum(w * tf.h(r) * np.array([float(theta(float(x))) for x in r]))) / (2 * math.pi)
    raise ValueError("route is 'series' or 'quad'")


def level_factor(tf: TestFunction, p: int, eps: int) -> float:
    """The prime p's share of the continuum: log p [g(0) + Σ_{k≥1} (−ε_p)^k p^{−k/2} g(k log p)]."""
    L = math.log(p)
    tot, k = float(tf.g(0.0)), 1
    while k * L <= tf.length_needed(1e-18):
        tot += (-eps) ** k * p ** (-k / 2) * float(tf.g(k * L))
        k += 1
    return L * tot


class Level:
    """The trace formulas of X₀(N) for N = 11 or 143 (frozen weights), sector by sector."""

    def __init__(self, N: int = 143):
        if N not in (11, 143):
            raise ValueError("the sector formulas are frozen for N = 11 and 143 (level 1 has its own check: level1_total)")
        self.N = N
        self.primes = prime_factors(N)
        self.n = 2 ** len(self.primes)
        index = 1
        for p in self.primes:
            index *= p + 1
        self.vol = math.pi / 3 * index
        ends = mirror_ends(N)
        self.cusps = len(ends) // 2
        self.sum_beta = 0.5 * sum(e["beta"] for e in ends)
        self.Qs = sorted(math.prod(self.primes[i] for i in range(len(self.primes)) if k >> i & 1) for k in range(1, 2 ** len(self.primes)))
        W = _weights()
        self.E = {Q: float(W[f"N{N}_Q{Q}_E"]) for Q in self.Qs}
        self.L = {Q: float(W[f"N{N}_Q{Q}_L"]) for Q in self.Qs}

    def eps_Q(self, eps, Q: int) -> int:
        return math.prod(e for p, e in zip(self.primes, eps) if Q % p == 0)

    def _sum(self, tf, Q, kind, twist=False, cut=None) -> float:
        _, length, w = classes(self.N, Q, kind, twist)
        if cut is None and tf.truncation(length_bound()) > TRUNCATION_ALLOWED:
            raise ValueError(f"sigma = {tf.sigma} needs geodesics longer than the frozen {length_bound()} (tail about {tf.truncation(length_bound()):.1e}): use sigma >= 0.42")
        if cut is not None:
            keep = length <= cut
            length, w = length[keep], w[keep]
        return float(np.sum(tf.g(length) * w))

    def base_parts(self, tf: TestFunction, cut=None) -> dict:
        """A, G, G′ and M of the module's note (the untwisted coset and its mirror)."""
        I = tf.integrals()
        return {"area": self.vol / (8 * math.pi) * I["r tanh"],
                "geodesics": 0.5 * self._sum(tf, 1, "hyp", cut=cut),
                "glides": 0.5 * self._sum(tf, 1, "gli", cut=cut),
                "cusp mirrors": (self.cusps / 2) * (I["g(0)"] * (math.log(2) + EULER) + 4 * I["K"] - 2 * I["J"]) - I["g(0)"] / 4 * self.sum_beta}

    def coset_parts(self, tf: TestFunction, Q: int, cut=None) -> dict:
        """The pieces of W_Q and of R W_Q, Q > 1: fixed points and geodesics; mirrors and glides."""
        I = tf.integrals()
        return {"fixed points": self.E[Q] * C_ELL * I["sech"], "geodesics": self._sum(tf, Q, "hyp", cut=cut),
                "mirrors": self.L[Q] * C_REFL * I["g(0)"], "glides": self._sum(tf, Q, "gli", cut=cut)}

    def odd(self, tf: TestFunction, eps, cut=None) -> float:
        """Σ h(r_j) over the odd cusp forms with Atkin–Lehner signs eps, from geometry."""
        b = self.base_parts(tf, cut)
        tot = (b["area"] + b["geodesics"] - b["glides"] + b["cusp mirrors"]) / self.n
        for Q in self.Qs:
            c = self.coset_parts(tf, Q, cut)
            tot += self.eps_Q(eps, Q) * (c["fixed points"] + c["geodesics"] - c["mirrors"] - c["glides"]) / (2 * self.n)
        return tot

    def even_pieces(self, tf: TestFunction, eps, cut=None) -> dict:
        """Every piece of Σ h(r_j) over the even cusp forms of the sector; they sum to it."""
        I = tf.integrals()
        b = self.base_parts(tf, cut)
        out = {"area": b["area"] / self.n, "geodesics": b["geodesics"] / self.n, "glides": b["glides"] / self.n,
               "cusp mirrors": -b["cusp mirrors"] / self.n,
               "cusps": I["h(0)"] / 2 - I["g(0)"] * math.log(2) - I["psi 1"]}
        for Q in self.Qs:
            c = self.coset_parts(tf, Q, cut)
            out[f"W{Q}"] = self.eps_Q(eps, Q) * (c["fixed points"] + c["geodesics"] + c["mirrors"] + c["glides"]) / (2 * self.n)
        out["continuum, level 1"] = -continuum_level1(tf)
        for p, e in zip(self.primes, eps):
            out[f"continuum, prime {p}"] = -level_factor(tf, p, e)
        out["constant function"] = -tf.h_i_half if all(e == 1 for e in eps) else 0.0
        return out

    def even(self, tf: TestFunction, eps, cut=None) -> float:
        return sum(self.even_pieces(tf, eps, cut).values())

    def twisted_total(self, tf: TestFunction, cut=None) -> dict:
        """Σ h(r_j) over all forms of the η-twist (level 143): area, classes signed by χ, and the closed cusps."""
        if self.N != 143:
            raise ValueError("the twisted weights are frozen for level 143")
        I = tf.integrals()
        out = {"area": self.vol / (4 * math.pi) * I["r tanh"], "geodesics": self._sum(tf, 1, "hyp", twist=True, cut=cut),
               "cusps": -self.cusps * I["g(0)"] * math.log(2)}
        out["total"] = sum(out.values())
        return out


def level1_total(tf: TestFunction, cut=None) -> dict:
    """PSL(2, ℤ)'s whole formula: every piece, and their sum, which is Σ h(r_j) over its cusp forms.

    The first cusp form has R = 9.5337, so for a test function centred at 0 the sum must vanish:
    a check of every constant and of the class numbers, with no Maass form involved.
    """
    if tf.truncation(length_bound()) > TRUNCATION_ALLOWED:
        raise ValueError(f"sigma = {tf.sigma} needs geodesics longer than the frozen {length_bound()}: use sigma >= 0.42")
    I = tf.integrals()
    _, length, w = classes(1, 1, "hyp")
    if cut is not None:
        length, w = length[length <= cut], w[length <= cut]
    p = {"area": (math.pi / 3) / (4 * math.pi) * I["r tanh"], "geodesics": float(np.sum(tf.g(length) * w)),
         "fixed point of order 2": 0.5 * C_ELL * I["sech"], "fixed point of order 3": I["order 3"] / (3 * math.sqrt(3)),
         "cusp": I["h(0)"] / 2 - I["g(0)"] * math.log(2) - I["psi 1"], "continuum": -continuum_level1(tf),
         "constant function": -tf.h_i_half}
    p["total"] = sum(p.values())
    return p


# ------------------------------------------------------------------------------------------------ what a list is worth
SECTORS = ((1, 1), (1, -1), (-1, 1), (-1, -1))
DENSITY = 3.5            # forms per unit of r·dr in one sector and parity: Weyl's 14 R² over eight


def completeness(residual: float, tf: TestFunction, slack: float = 1e-6):
    """The radius below which no form can be missing, from geometry − list = residual (r₀ = 0 only).

    A missing form at R would add h(R) = 2 exp(−R²/2σ²) to the residual.  Returns None when the
    residual is negative beyond the slack (the list holds something the geometry does not).
    """
    if tf.r0 != 0.0:
        raise ValueError("the positivity radius is implemented for test functions centred at 0")
    if residual < -slack:
        return None
    x = (max(residual, 0.0) + slack) / 2
    return tf.sigma * math.sqrt(2 * math.log(1 / x)) if x < 1 else 0.0


def check_lists(lists: dict, tf: TestFunction, N: int = 143, r_top: float = 2.0, cut=None) -> dict:
    """Compare lists of R by (parity, ε) with the geometric side.

    ``lists[(parity, eps)]`` is the list of R of that sector (parity 0 even, 1 odd; eps a tuple of
    signs, one per prime of N).  For each: the geometric sum, the list's sum, the residual
    (geometry − list), whether it is positive within 1e-6, the positivity radius (r₀ = 0), and
    ``expected``, the smooth count of what forms above ``r_top`` would add.  ``expected`` is an
    average, not a bound: one form just above the list's top can exceed it.
    """
    lev = Level(N)
    expected = tf.smooth_tail(r_top, DENSITY)
    out = {}
    for (parity, eps), Rs in lists.items():
        geo = lev.even(tf, eps, cut) if parity == 0 else lev.odd(tf, eps, cut)
        lst = float(np.sum(tf.h(np.array(Rs, dtype=float)))) if len(Rs) else 0.0
        res = geo - lst
        out[(parity, eps)] = {"geometry": geo, "list": lst, "residual": res, "expected": expected, "positive": res >= -1e-6,
                              "radius": completeness(res, tf) if tf.r0 == 0.0 else None}
    return out


def table_lists(twist: bool = False) -> dict:
    """The frozen table's R by (parity, ε) at level 143 (``mtft.maass.spectrum``)."""
    from . import spectrum
    out = {(par, e): [] for par in (0, 1) for e in SECTORS}
    for f in spectrum.forms(twist=twist):
        out[(f["parity"], f["eps"])].append(f["R"])
    return out


def check_table(tf: TestFunction = None, cut=None) -> dict:
    """The frozen untwisted table against the trace formula: sector by sector, with the positivity radius."""
    return check_lists(table_lists(False), tf or TestFunction(0.0, 0.5), 143, 2.0, cut)


def check_twisted(tf: TestFunction = None, cut=None) -> dict:
    """The frozen twisted table, all sectors together, against the twisted formula."""
    tf = tf or TestFunction(0.0, 0.5)
    geo = Level(143).twisted_total(tf, cut)["total"]
    Rs = np.array([R for v in table_lists(True).values() for R in v])
    lst = float(np.sum(tf.h(Rs)))
    res = geo - lst
    return {"geometry": geo, "list": lst, "residual": res, "expected": tf.smooth_tail(1.8, 8 * DENSITY), "positive": res >= -1e-6,
            "radius": completeness(res, tf) if tf.r0 == 0.0 else None, "forms": int(Rs.size)}


def certificate(sigmas=(0.43, 0.45, 0.5)) -> dict:
    """What the frozen table of Γ₀(143) is worth by the package's own trace formulas.

    Each width σ gives a valid radius for each sector and parity; the certified radius is the
    best of them, and the table's is the smallest over the eight.  Also: every residual is
    positive, and none comes near 2, which an eigenvalue below ¼ would add (h(it) ≥ 2 for r₀ = 0),
    so the formula sees no exceptional eigenvalue.  The class numbers behind the weights are PARI's.
    """
    per, tw_best, residuals = {}, 0.0, []
    for sg in sigmas:
        tf = TestFunction(0.0, sg)
        for key, v in check_table(tf).items():
            residuals.append(v["residual"])
            if v["radius"] is not None and v["radius"] > per.get(key, {"radius": 0.0})["radius"]:
                per[key] = {"radius": v["radius"], "sigma": sg, "residual": v["residual"]}
        tw = check_twisted(tf)
        residuals.append(tw["residual"])
        tw_best = max(tw_best, tw["radius"] or 0.0)
    return {"sectors": per, "untwisted_radius": min(v["radius"] for v in per.values()) if len(per) == 8 else 0.0,
            "twisted_radius": tw_best, "all_positive": min(residuals) >= -1e-6, "largest_residual": max(residuals),
            "no_exceptional_eigenvalue": max(residuals) < 2.0, "sigmas": tuple(sigmas), "lengths_to": length_bound()}
