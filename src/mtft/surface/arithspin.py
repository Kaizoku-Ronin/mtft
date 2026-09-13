"""mtft.surface.arithspin — arithmetic spin structures, cuspidal group, CM fixed points, three-family purity
(ARITH-SPIN-01 / GAL-01, v0.29.2).  GP-free.

Exact facts for X0(143) established in this module's functions:
- eta-quotient modular units (Newman conditions, Ligozat cusp orders) generate the cuspidal relation lattice
  (for squarefree N all modular units are eta quotients); cuspidal class group Z/420 x Z/10 (order 4200);
  orders of (0)-(oo), (1/11)-(oo), (1/13)-(oo) are 420, 60, 70.
- K ~ 6 D_cusps.  The weight-1 CM form g1 of character chi_{-143} has purely cuspidal divisor (1,6,6,1).
  Cuspidal theta characteristics: S1 = O(3 D) (h^0 = 5, odd), S0 = O(6(0)+6(1/11)) (h^0 = 2, even), and two
  more differing by the cuspidal 2-torsion; all AL-invariant.  h^0 by Serre duality from AL-signed cusp orders
  of the frozen weight-2 eigenbasis.
- u = eta(13t)eta(143t)/(eta(t)eta(11t)) spans H^0(S0) with the constant; u∘W13 = C/u, C = -1/13; u∘W11 = -u.
- W13-fixed points are the four CM points of discriminant -52; u(P) = ±i/sqrt(13), two of each sign;
  hence L3 = S0 (x) O(P1+P2+P3) has h^0 = 3, h^1 = 0 for every triple (Brill–Noether-pure three-family flux),
  W13-symmetric, W11-breaking, defined over Q(i, sqrt 13).
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from math import gcd, lcm
from typing import Dict, List

import mpmath as mp
import numpy as np
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

CUSPS_143 = ("0", "1/11", "1/13", "oo")          # cusp order by denominator class c = 1, 11, 13, 143
AL_ON_CUSPS = {"W11": (1, 0, 3, 2), "W13": (2, 3, 0, 1), "W143": (3, 2, 1, 0)}
AL_MATRICES = {"W11": (66, 5, 143, 11), "W13": (78, 7, 143, 13), "W143": (0, -1, 143, 0)}


def divisors(n): return [d for d in range(1, n + 1) if n % d == 0]


def eta_quotient_lattice(N: int = 143, box: int = 12) -> Dict:
    """Cusp-order vectors of all eta quotients on Gamma_0(N) with exponents in [-box, box] (Newman conditions,
    Ligozat orders), and the Hermite basis of the lattice they generate."""
    D = divisors(N); C = D                                       # cusp representatives 1/c, c | N
    W = [[Fraction(N * gcd(d, c) ** 2, 24 * gcd(c, N // c) * c * d) for d in D] for c in C]
    L = lcm(*[w.denominator for row in W for w in row]); Wi = np.array([[int(w * L) for w in row] for row in W], dtype=np.int64)
    R = np.array(list(itertools.product(range(-box, box + 1), repeat=len(D))), dtype=np.int64)
    ok = R.sum(axis=1) == 0
    ok &= (R @ np.array(D)) % 24 == 0
    ok &= (R @ np.array([N // d for d in D])) % 24 == 0
    # product of d^{r_d} a square: for each prime p | N the total p-adic valuation must be even
    def vp(n, p):
        k = 0
        while n % p == 0: n //= p; k += 1
        return k
    primes = [p for p in range(2, N + 1) if N % p == 0 and all(p % k for k in range(2, int(p ** 0.5) + 1))]
    for p in primes:
        ok &= (R @ np.array([vp(d, p) for d in D])) % 2 == 0
    R = R[ok]
    O = R @ Wi.T
    integral = np.all(O % L == 0, axis=1)
    orders = (O[integral] // L)
    orders = orders[np.any(orders != 0, axis=1)]
    exps = R[integral][np.any((R[integral] @ Wi.T) // L != 0, axis=1)]
    H = hermite_normal_form(Matrix(orders.T.tolist()))
    H = Matrix([[int(x) for x in row] for row in H.tolist()])
    return {"divisors": D, "exponents": exps, "orders": orders, "hermite": H}


def in_lattice(H: Matrix, w) -> bool:
    s = H.solve_least_squares(Matrix(list(w)))
    return (H * s == Matrix(list(w))) and all(x.is_integer for x in s)


def cuspidal_group(H: Matrix) -> Dict:
    """Invariant factors of Z^{#cusps}_0 / Lambda (Smith normal form in the degree-0 basis)."""
    n = H.rows
    Z0 = Matrix([[1 if i == j else (-1 if i == j + 1 else 0) for j in range(n - 1)] for i in range(n)])   # (c_j)-(c_{j+1})
    X = (Z0.T * Z0).inv() * (Z0.T * H)
    S = smith_normal_form(Matrix([[int(x) for x in row] for row in X.tolist()]), domain=ZZ)
    inv = [abs(int(S[i, i])) for i in range(min(S.shape)) if S[i, i] != 0]
    return {"invariant_factors": inv, "order": int(np.prod(inv))}


def cusp_difference_order(H: Matrix, i: int, j: int, nmax: int = 5000):
    w = [0] * H.rows; w[i] = 1; w[j] = -1
    for n in range(1, nmax):
        if in_lattice(H, [n * x for x in w]): return n
    return None


def x0143_theta_characteristics(H: Matrix) -> Dict:
    """The four cuspidal theta characteristics of X0(143) and their AL invariance."""
    K = [1, 11, 11, 1]; Dc = [1, 1, 1, 1]; S0 = [6, 6, 0, 0]
    assert in_lattice(H, [6 * x - k for x, k in zip(Dc, K)]), "K ~ 6 D_cusps expected"
    tors = []
    for s in itertools.product((0, 1), repeat=H.cols):
        v = H * Matrix(s)
        if all(x % 2 == 0 for x in v):
            t = [int(x) // 2 for x in v]
            if not in_lattice(H, t) and not any(in_lattice(H, [a - b for a, b in zip(t, u)]) for u in tors): tors.append(t)
    cands = [S0] + [[a + b for a, b in zip(S0, t)] for t in tors]
    out = []
    for A in cands:
        inv = {name: in_lattice(H, [A[p[k]] - A[k] for k in range(4)]) for name, p in AL_ON_CUSPS.items()}
        out.append({"class": A, "al_invariant": inv})
    return {"K_cuspidal": K, "two_torsion": tors, "characteristics": out}


def effective_representative(H: Matrix, A, box: int = 6):
    if min(A) >= 0: return list(A)
    best = None
    for c in itertools.product(range(-box, box + 1), repeat=H.cols):
        B = [A[i] + int((H * Matrix(c))[i]) for i in range(len(A))]
        if min(B) >= 0 and (best is None or max(B) < max(best)): best = B
    return best


def h0_cuspidal_class(A) -> int:
    """h^0(O(A)) for an effective cuspidal divisor A of degree 12 on X0(143) by Serre duality:
    h^0(O(A)) = h^0(K - A) = #{weight-2 cusp forms with ord_c >= A_c + 1 at every cusp}, with orders at the
    cusps W_Q(oo) read from the AL eigenvalues of the frozen eigenbasis."""
    from .yukawa import load_basis, al_eigenbasis
    forms, labels, dens, vecs = al_eigenbasis(load_basis())
    eps = {"oo": [1] * 13, "1/13": [l[0] for l in labels], "1/11": [l[1] for l in labels], "0": [l[0] * l[1] for l in labels]}
    rows = [[eps[nm][i] * int(forms[i][n]) for i in range(13)] for ci, nm in enumerate(CUSPS_143) for n in range(1, A[ci] + 1)]
    return 13 - (Matrix(rows).rank() if rows else 0)


# ---------------------------------------------------------------- modular units and CM fixed points
def eta(tau, terms: int = 800):
    tau = mp.mpc(tau.real, tau.imag); q = mp.exp(2j * mp.pi * tau); v = mp.exp(2j * mp.pi * tau / 24)   # branch: exp(2 pi i tau/24)
    for n in range(1, terms): v *= (1 - q ** n)
    return v


def eta_quotient(tau, exponents, N: int = 143, terms: int = 800):
    v = mp.mpc(1)
    for d, r in zip(divisors(N), exponents): v *= eta(d * tau, terms) ** r
    return v


UNIT_U = (-1, -1, 1, 1)                                                # u = eta(13t)eta(143t)/(eta(t)eta(11t))


def mobius(M, t): a, b, c, d = M; return (a * t + b) / (c * t + d)


def functional_constants(N: int = 143) -> Dict:
    """C with u∘W13 = C/u and c' with u∘W11 = c' u, at points where both images have Im >~ 0.02."""
    t13 = complex(-13 / 143 + 0.004, 0.028); t11 = complex(-11 / 143 + 0.003, 0.0277)
    C = complex(eta_quotient(mobius(AL_MATRICES["W13"], t13), UNIT_U, N, 1500) * eta_quotient(t13, UNIT_U, N, 1500))
    cp = complex(eta_quotient(mobius(AL_MATRICES["W11"], t11), UNIT_U, N, 1500) / eta_quotient(t11, UNIT_U, N, 1500))
    return {"C_W13": C, "c_W11": cp}


def w13_fixed_points(N: int = 143, amax: int = 40, im_min: float = 0.02) -> List[complex]:
    """Representatives of the W13-fixed CM points (discriminant -52): tau with 143 c tau^2 - 26 a tau - b = 0,
    13 a^2 + 11 b c = -1; only representatives with Im tau >= im_min (reliable q-series evaluation)."""
    pts = []
    for a in range(-amax, amax + 1):
        for c in range(1, amax + 1):
            num = -1 - 13 * a * a
            if num % (11 * c) == 0:
                t = (26 * a + 1j * np.sqrt(52.0)) / (2 * 143 * c); t = complex(t.real - np.round(t.real), t.imag)
                if t.imag >= im_min: pts.append(t)
    return pts


def three_family_purity() -> Dict:
    """u(P)^2 = -1/13 at the W13-fixed points and two of each sign => h^0(S0 (x) O(P1+P2+P3)) = 3, h^1 = 0."""
    vals = [complex(eta_quotient(t, UNIT_U, 143, 1500)) for t in w13_fixed_points()]
    sq = np.array([v * v for v in vals]); signs = sorted(set(int(np.sign(v.imag)) for v in vals))
    return {"u_squared_max_dev_from_-1/13": float(np.max(np.abs(sq + 1 / 13))), "signs_present": signs,
            "h0_L3": 3 if signs == [-1, 1] else None, "h1_L3": 0 if signs == [-1, 1] else None,
            "note": "purity holds iff both signs occur among the four fixed points (u∘W11 = -u pairs them)"}
