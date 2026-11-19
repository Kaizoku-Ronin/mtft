r"""mtft.lines.genus — genus characters, and why the systole of X₀(143) winds round one cusp only (EXACT; KN-1, 2026-11-19).

A hyperbolic A = [[a, b], [c, d]] of trace t has the form Q_A = [c, d − a, −b] of discriminant t² − 4,
and conjugation in SL(2, ℤ) moves Q_A inside its proper class.  For A in Γ₀(N) and δ | N put
A_δ = [[a, δb], [c/δ, d]]: the image of the closed line under τ ↦ δτ, the matrix that
``mtft.lines.rademacher.winding`` evaluates.  For a prime p | c that does not divide d − a,
Dirichlet composition gives

    [c, B, C] = [c/p, B, pC] ∘ [p, B, (c/p)C],

so the class of Q_{A_p} is the class of Q_A times the inverse class of a prime ideal of norm p.
When every genus of the discriminant holds one narrow class, that class is fixed by the genus
characters χ_{d_i}(p), d_i the prime discriminants of the fundamental discriminant.

Trace 4: the systole.  D = 12 = (−3)(−4), h⁺(12) = 2, one class per genus: x² − 3y² (the class of
R²L, Ψ = +1) and −x² + 3y² (the class of RL², Ψ = −1).  For p ≢ 2, 3 the congruence (a − 2)² ≡ 3
(mod p) keeps p off d − a.  So for every trace-4 element of Γ₀(N) and every prime p | N,

    Ψ(A_p) = χ₋₃(p) Ψ(A),       χ₋₃(p) = +1 for p ≡ 1 (mod 3) and −1 for p ≡ 2 (mod 3),

multiplicatively in δ, and the winding of a unit ∏ η(δτ)^{r_δ} along the line is

    φ(A) = (Ψ(A)/24) Σ_δ r_δ χ₋₃(δ).

At 143: (χ₋₃(1), χ₋₃(11), χ₋₃(13), χ₋₃(143)) = (1, −1, 1, −1), so the shortest closed geodesics
wind Ψ = ±1 times round the divisor of u_A = (η(τ)/η(11τ))¹² and not at all round those of φ₁ and
u_B, because 11 ≡ 2 and 13 ≡ 1 (mod 3).  Trace-4 elements exist in Γ₀(p) only when 3 is a square
mod p, p ≡ ±1 (mod 12).

Route 1 is Ψ by Dedekind sums (``mtft.lines.rademacher.Psi``); route 2 is the composition law above.
They agreed on 1193 trace-4 elements at the twenty primes p < 200 with p ≡ ±1 (mod 12), and on all
62 found at 143 (record KN-1; ``check_trace4`` reruns it).  For longer lines Ψ is not a class
character (trace 24 at 143 already shows Ψ = −21 and −9), and no law of this kind is claimed.
"""
from __future__ import annotations

import math

from . import rademacher as _r


def _factor(n: int) -> dict:
    n, out, p = abs(n), {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def jacobi(a: int, n: int) -> int:
    """The Jacobi symbol (a/n) for odd n > 0."""
    if n <= 0 or n % 2 == 0:
        raise ValueError("the Jacobi symbol needs an odd positive n")
    a %= n
    s = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                s = -s
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            s = -s
        a %= n
    return s if n == 1 else 0


def kronecker(a: int, n: int) -> int:
    """The Kronecker symbol (a/n) for n ≥ 1."""
    if n < 1:
        raise ValueError("n must be positive")
    s = 1
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0:
            return 0
        if a % 8 in (3, 5):
            s = -s
    return s * (jacobi(a, n) if n > 1 else 1)


def is_fundamental(D: int) -> bool:
    """Whether D is a fundamental discriminant."""
    if D in (0, 1):
        return False
    if D % 4 == 1:
        return all(e == 1 for e in _factor(D).values())
    if D % 4 == 0:
        m = D // 4
        return m % 4 in (2, 3) and all(e == 1 for e in _factor(m).values())
    return False


def prime_discriminants(D: int) -> tuple:
    """The prime discriminants whose product is the fundamental discriminant D (odd ones first, by size)."""
    if not is_fundamental(D):
        raise ValueError(f"{D} is not a fundamental discriminant")
    odd = [p if p % 4 == 1 else -p for p in sorted(q for q in _factor(D) if q != 2)]
    rest = D
    for d in odd:
        rest //= d
    out = tuple(odd)
    if rest != 1:
        if rest not in (-4, 8, -8):
            raise ArithmeticError("the even part is not a prime discriminant")
        out += (rest,)
    return out


def genus_vector(D: int, n: int) -> tuple:
    """(χ_{d_i}(n)) over the prime discriminants d_i of D, for n coprime to D: the genus of an ideal of norm n."""
    if math.gcd(D, n) != 1:
        raise ValueError("n must be coprime to the discriminant")
    return tuple(kronecker(d, abs(n)) for d in prime_discriminants(D))


def chi_minus3(n: int) -> int:
    """χ₋₃(n) = (−3/n): +1, −1 or 0 as n ≡ 1, 2 or 0 (mod 3)."""
    return kronecker(-3, abs(n))


def trace4_ratio(delta: int) -> int:
    """Ψ(A_δ)/Ψ(A) on the trace-4 elements of Γ₀(N), δ | N coprime to 6: the product of χ₋₃ over the primes of δ."""
    if delta < 1 or math.gcd(delta, 6) != 1:
        raise ValueError("δ must be a positive integer coprime to 6")
    out = 1
    for p, e in _factor(delta).items():
        out *= chi_minus3(p) ** e
    return out


def trace4_windings(psi: int, units: dict | None = None, order=None) -> tuple:
    """The windings round the cusps of a trace-4 line with Ψ = psi: (ψ/24) Σ_δ r_δ χ₋₃(δ) for each unit (default: 143's)."""
    units = _r.UNITS if units is None else units
    order = _r.UNIT_ORDER if order is None else order
    out = []
    for k in order:
        v = sum(e * trace4_ratio(dl) for dl, e in units[k].items())
        if (psi * v) % 24:
            raise ArithmeticError("the predicted winding is not an integer")
        out.append(psi * v // 24)
    return tuple(out)


def trace4_elements(N: int, kmax: int = 4, amax: int = 400) -> list:
    """The trace-4 elements [[a, b], [c, d]] of Γ₀(N) with c = Nk (1 ≤ k ≤ kmax) and |a| ≤ amax."""
    if N < 1 or kmax < 1 or amax < 1:
        raise ValueError("N, kmax and amax must be positive")
    out = []
    for k in range(1, kmax + 1):
        c = N * k
        for a in range(-amax, amax + 1):
            d = 4 - a
            if (a * d - 1) % c == 0:
                out.append(((a, (a * d - 1) // c), (c, d)))
    return out


def check_trace4(N: int, kmax: int = 4, amax: int = 400) -> dict:
    """Route 1 (Dedekind sums) against route 2 (genus characters) on the trace-4 elements of Γ₀(N), N coprime to 6."""
    els = trace4_elements(N, kmax, amax)
    divs = [dl for dl in range(1, N + 1) if N % dl == 0]
    bad = 0
    for g in els:
        (a, b), (c, d) = g
        psi = _r.Psi(g)
        if psi not in (1, -1):
            bad += 1
            continue
        for dl in divs:
            if _r.Psi(((a, dl * b), (c // dl, d))) != trace4_ratio(dl) * psi:
                bad += 1
                break
    return {"N": N, "elements": len(els), "disagreements": bad, "ratios": {dl: trace4_ratio(dl) for dl in divs}}
