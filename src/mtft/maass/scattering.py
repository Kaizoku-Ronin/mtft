"""mtft.maass.scattering — the scattering matrix and Eisenstein series of Γ₀(N), N squarefree (EXACT formulas; MA-01, 2026-11-15).

Level 1.  E(z, s) = Σ_{Γ_∞\\PSL(2,ℤ)} Im(γz)^s = y^s + φ(s) y^{1−s} + …, with

    φ(s) = ξ(2s − 1)/ξ(2s) = √π Γ(s − ½) ζ(2s − 1) / (Γ(s) ζ(2s)),     ξ(s) = π^{−s/2} Γ(s/2) ζ(s).

Level N.  Each prime p | N has two cusp types, ∞ and 0, and a cusp of Γ₀(N) is a choice of type at
every prime: at N = 143 the four cusps are ∞ = (∞, ∞), 1/11 = (∞, 0), 1/13 = (0, ∞) and 0 = (0, 0),
the types listed at 11 and then at 13.  The series of a cusp is built from E(dz), d | N, by the rule
f ↦ (p^s f(pz) − f(z))/(p^{2s} − 1) for type ∞ and f ↦ (p^s f(z) − f(pz))/(p^{2s} − 1) for type 0, and

    Φ(s) = φ(s) ⊗_p M_p(s),     M_p(s) = (p^{2s} − 1)⁻¹ [[p − 1, p^s − p^{1−s}], [p^s − p^{1−s}, p − 1]].

Facts gated in the tests (the lab's record MP14, check B1): Φ(s)Φ(1 − s) = I; Φ is unitary on
Re s = ½; every entry has residue 1/area at s = 1 (3/(168π) at N = 143); and Φ(½) = −I, because
φ(½) = −1 and M_p(½) = I.  So every Eisenstein series vanishes identically at the bottom of the
continuum, s = ½.

All of this is established mathematics, reproduced.  Arithmetic is mpmath's, at the working
precision the caller sets (``mp.mp.dps``).
"""
from __future__ import annotations

import itertools

import mpmath as mp

from .domain import is_squarefree, prime_factors


def xi(s):
    """The completed zeta function ξ(s) = π^{−s/2} Γ(s/2) ζ(s)."""
    return mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def phi(s):
    """The scattering function of PSL(2, ℤ): φ(s) = ξ(2s − 1)/ξ(2s)."""
    s = mp.mpmathify(s)
    return mp.sqrt(mp.pi) * mp.gamma(s - mp.mpf(1) / 2) * mp.zeta(2 * s - 1) / (mp.gamma(s) * mp.zeta(2 * s))


def local_matrix(p: int, s):
    """M_p(s): rows and columns in the order (∞, 0) at the prime p."""
    s = mp.mpmathify(s)
    d = mp.mpf(p) ** (2 * s) - 1
    a, b = (p - 1) / d, (mp.mpf(p) ** s - mp.mpf(p) ** (1 - s)) / d
    return mp.matrix([[a, b], [b, a]])


def cusps(N: int = 143) -> list:
    """The cusps of Γ₀(N) in the order of the scattering matrix's rows: a type (∞ or 0) per prime.

    Each is named by a representative 1/d: d is the product of the primes whose type is ∞,
    except that d = N is written "inf" and d = 1 is written "0".
    """
    if not is_squarefree(N):
        raise ValueError("the scattering matrix is implemented for squarefree N")
    ps = prime_factors(N)
    out = []
    for types in itertools.product((1, 0), repeat=len(ps)):          # 1: type ∞ at that prime
        d = 1
        for p, t in zip(ps, types):
            d *= p if t else 1
        out.append("inf" if d == N else ("0" if d == 1 else f"1/{d}"))
    return out


def scattering_matrix(s, N: int = 143):
    """Φ(s) of Γ₀(N): φ(s) times the Kronecker product of M_p(s) over p | N, cusps as in ``cusps(N)``."""
    if not is_squarefree(N):
        raise ValueError("the scattering matrix is implemented for squarefree N")
    s = mp.mpmathify(s)
    if s == mp.mpf(1) / 2:                                 # φ(½) = −1 and M_p(½) = I, by the limit
        return -mp.eye(2 ** len(prime_factors(N)))
    out = mp.matrix([[1]])
    for p in prime_factors(N):
        B = local_matrix(p, s)
        K = mp.matrix(out.rows * 2, out.cols * 2)
        for i in range(out.rows):
            for j in range(out.cols):
                for k in range(2):
                    for l in range(2):
                        K[2 * i + k, 2 * j + l] = out[i, j] * B[k, l]
        out = K
    return phi(s) * out


def _reduce_level1(z):
    """z carried into the standard domain of PSL(2, ℤ)."""
    z = mp.mpc(z)
    for _ in range(10000):
        z = mp.mpc(z.real - mp.floor(z.real + mp.mpf(1) / 2), z.imag)
        if abs(z) >= 1 - mp.mpf(10) ** (-mp.mp.dps + 3):
            return z
        z = -1 / z
    raise RuntimeError("the level-1 reduction did not converge")


def eisenstein_level1(z, s):
    """E(z, s) of PSL(2, ℤ), by its Fourier expansion at the reduced point."""
    s = mp.mpmathify(s)
    w = _reduce_level1(z)
    x, y = w.real, w.imag
    nu = s - mp.mpf(1) / 2
    out = y ** s + phi(s) * y ** (1 - s)
    c = 4 / xi(2 * s)
    n, tot = 1, mp.mpc(0)
    while True:
        k = mp.besselk(nu, 2 * mp.pi * n * y)
        div = mp.fsum(mp.mpf(d) ** (1 - 2 * s) for d in range(1, n + 1) if n % d == 0)
        tot += n ** nu * div * mp.sqrt(y) * k * mp.cos(2 * mp.pi * n * x)
        if n > 3 and abs(k) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
        n += 1
    return out + c * tot


def eisenstein(z, s, N: int = 143) -> dict:
    """The Eisenstein series of every cusp of Γ₀(N) at z: a dictionary keyed as ``cusps(N)``."""
    if not is_squarefree(N):
        raise ValueError("the Eisenstein series are implemented for squarefree N")
    s = mp.mpmathify(s)
    z = mp.mpc(z)
    ps = prime_factors(N)
    base = {}                                                          # E(d z, s), d | N
    for k in range(2 ** len(ps)):
        d = 1
        for i, p in enumerate(ps):
            d *= p if k >> i & 1 else 1
        base[d] = eisenstein_level1(d * z, s)
    den = mp.mpf(1)
    for p in ps:
        den *= mp.mpf(p) ** (2 * s) - 1
    out = {}
    for name, types in zip(cusps(N), itertools.product((1, 0), repeat=len(ps))):
        terms = {1: mp.mpf(1)}                                         # coefficient of E(d z, s)
        for p, t in zip(ps, types):
            ps_ = mp.mpf(p) ** s
            new = {}
            for d, co in terms.items():
                if t:                                                  # type ∞: p^s f(pz) − f(z)
                    new[d * p] = new.get(d * p, 0) + ps_ * co
                    new[d] = new.get(d, 0) - co
                else:                                                  # type 0: p^s f(z) − f(pz)
                    new[d] = new.get(d, 0) + ps_ * co
                    new[d * p] = new.get(d * p, 0) - co
            terms = new
        out[name] = mp.fsum(co * base[d] for d, co in terms.items()) / den
    return out
