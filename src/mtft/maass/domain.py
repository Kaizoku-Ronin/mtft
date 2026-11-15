"""mtft.maass.domain — carrying a point into the Ford domain of Γ₀*(N), with its characters (EXACT rule).

Γ₀*(N), for squarefree N, is Γ₀(N) with all its Atkin–Lehner matrices: the integer matrices
[[Q a, b], [N c, Q d]] of determinant Q, for the divisors Q of N.  A function with Atkin–Lehner
signs ε_p (p | N) satisfies f(Mz) = ε_Q f(z) for M of determinant Q, ε_Q the product of ε_p over
p | Q.  ``reduce_point`` applies, again and again, the element that raises the point most, until
none raises it; it returns the reduced point, the parity of the number of steps that used each
prime (so f(z) = ∏ ε_p^{parity_p} f(z*)), and, for odd N, the η-twist's sign.

The η-twist.  χ = ν_η¹², the multiplier of η¹², is the sign of the permutation that a matrix
makes of the three points of P¹(𝔽₂).  It extends to the integer matrices of odd determinant
(``sign_mod2``), so for odd N to Γ₀*(N).  At N = 143 all four cusp widths (1, 11, 13, 143) are
odd, so χ = −1 round every cusp and the twisted space has no continuous spectrum (record MP18).

``Horocycle`` holds the sample points of Hejhal's method on a closed horocycle y = Y with their
reduced points and characters.
"""
from __future__ import annotations

import itertools
import math

import numpy as np


def prime_factors(N: int) -> list:
    """The primes dividing N, in increasing order."""
    if N < 1:
        raise ValueError("N must be a positive integer")
    out, p, n = [], 2, N
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def is_squarefree(N: int) -> bool:
    return all(N % (p * p) for p in prime_factors(N))


def egcd(a: int, b: int):
    """(g, u, v) with u a + v b = g = gcd(a, b) ≥ 0."""
    if b == 0:
        return (a, 1, 0) if a >= 0 else (-a, -1, 0)
    g, u, v = egcd(b, a % b)
    return g, v, u - (a // b) * v


def sign_mod2(a: int, b: int, c: int, d: int) -> int:
    """χ̃ of an integer matrix of odd determinant: the sign of its permutation of P¹(𝔽₂).

    In GL(2, 𝔽₂) = S₃ the identity and the two elements of trace 1 (the 3-cycles) give +1; the
    three elements of trace 0 and order 2 give −1.
    """
    a, b, c, d = a & 1, b & 1, c & 1, d & 1
    if (a * d - b * c) & 1 == 0:
        raise ValueError("sign_mod2 needs a matrix of odd determinant")
    if (a, b, c, d) == (1, 0, 0, 1):
        return 1
    return 1 if (a + d) & 1 else -1


def reduce_point(x: float, y: float, N: int = 143, use_al: bool = True):
    """Carry z = x + iy into the Ford domain of Γ₀*(N) (of Γ₀(N) if ``use_al`` is False).

    Returns ``(x*, y*, parities, sign)``: ``parities`` has one entry per prime of N (the parity
    of the number of steps whose Q that prime divides), and ``sign`` is the η-twist's sign of the
    whole motion z ↦ z* for odd N (always +1 for even N, where the twist is not defined).
    """
    if y <= 0:
        raise ValueError("the point must lie in the upper half-plane")
    ps = prime_factors(N)
    if use_al and not is_squarefree(N):
        raise ValueError("the Atkin-Lehner reduction is implemented for squarefree N")
    Qs = [1]
    if use_al:
        Qs = sorted({math.prod(s) for k in range(len(ps) + 1) for s in itertools.combinations(ps, k)})
    odd = N % 2 == 1
    par = [0] * len(ps)
    sign = 1
    for _ in range(10000):
        k = math.floor(x + 0.5)
        x -= k
        if odd and k & 1:
            sign = -sign
        best, arg = 1.0 - 1e-13, None
        for Q in Qs:
            Nq = N // Q
            cmax = int(math.sqrt(Q) / (N * y)) + 1
            for c in range(1, cmax + 1):
                fl = math.floor(-N * c * x / Q)
                for d in range(fl - 1, fl + 3):
                    if math.gcd(Q * d, Nq * c) != 1:
                        continue
                    v = ((N * c * x + Q * d) ** 2 + (N * c * y) ** 2) / Q
                    if v < best:
                        best, arg = v, (Q, c, d)
        if arg is None:
            return x, y, tuple(par), sign
        Q, c, d = arg
        _, u, w = egcd(Q * d, (N // Q) * c)              # u Q d + w (N/Q) c = 1
        a, b = u, -w                                       # det [[Q a, b], [N c, Q d]] = Q
        den = (N * c * x + Q * d) ** 2 + (N * c * y) ** 2
        xn = ((Q * a * x + b) * (N * c * x + Q * d) + Q * a * N * c * y * y) / den
        y = Q * y / den
        x = xn
        if odd:
            sign *= sign_mod2(Q * a, b, N * c, Q * d)
        for i, p in enumerate(ps):
            if Q % p == 0:
                par[i] ^= 1
    raise RuntimeError("reduce_point did not converge")


class Horocycle:
    """Sample points x_m = (m − ½)/(2Q), m = 1 … Q, on y = Y, with their reduced points.

    Only half of the horocycle is sampled: parity (cos or sin in x) supplies the other half.
    ``Y`` must lie below the Ford domain, so that every sample point is moved by the reduction.
    """

    def __init__(self, N: int, Y: float, Q: int, use_al: bool = True):
        self.N, self.Y, self.Q, self.use_al = N, Y, Q, use_al
        self.xm = (np.arange(1, Q + 1) - 0.5) / (2 * Q)
        red = [reduce_point(float(x), Y, N, use_al) for x in self.xm]
        self.xs = np.array([r[0] for r in red])
        self.ys = np.array([r[1] for r in red])
        npr = len(prime_factors(N))
        self.par = np.array([list(r[2]) for r in red], dtype=int).reshape(Q, npr)
        self.sign = np.array([r[3] for r in red], dtype=float)
        if not self.ys.min() > Y * 1.0001:
            raise ValueError("the horocycle is not below the Ford domain: a sample point was left where it was")

    def chi(self, eps=(), twist: bool = False) -> np.ndarray:
        """The character of each sample's motion: ∏ ε_p^{parity_p}, times the η-twist's sign if asked."""
        eps = tuple(eps)
        if len(eps) != self.par.shape[1]:
            raise ValueError(f"eps needs one sign for each of the {self.par.shape[1]} primes of N = {self.N}")
        if twist and self.N % 2 == 0:
            raise ValueError("the eta-twist is defined here for odd N only")
        al = np.ones(self.Q) if not eps else np.prod(np.where(self.par == 1, np.array(eps, dtype=float)[None, :], 1.0), axis=1)
        return al * self.sign if twist else al
