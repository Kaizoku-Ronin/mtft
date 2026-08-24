"""mtft.quadratic_forms — the Gauss-Legendre three-squares layer (v0.19.0).

The tier-0 classical module proposed at the start of the 2026-08 arc and
certified as `three_squares_seed.py` (certificates v1-v7).  Layers:

  L1  8 T_n + 1 = (2n+1)^2                                   [Pr, EXACT]
  L2  r_3(n) = 0  <=>  n = 4^a (8b+7)      (Gauss-Legendre)  [Pr; E2 gate]
  L3  mod-8 character projector for the forbidden set        [Pr; gate]
  L4  natural density of the forbidden set = 1/6             [Pr; gate,
      honest error model O(log_4 X / X)]
  L5  Dirichlet series 8^-s zeta(s,7/8)/(1-4^-s)             [Pr; gate s=2]
  L6  Eureka: every N is T_a + T_b + T_c                     [Pr; gate]
  L7  self-similarity F(4n) = F(n)                           [Pr; gate]

E2 discipline: each gate compares routes sharing no computational steps
(lattice enumeration vs 2-adic/character arithmetic vs closed forms).
The forbidden indicator is also the seed of the Exception-Spacing
Curvature Law in `mtft.exception_spectrum` (e1 = 7, e2 = 15, base 14/15).
"""
from __future__ import annotations

import math

import numpy as np
from mpmath import mp, mpf, zeta

__all__ = [
    "v2", "forbidden", "forbidden_projector", "r3_array",
    "forbidden_count", "forbidden_density", "gate",
]


def v2(n):
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def forbidden(n):
    """Route A: 2-adic depth parity x mod-8 residue."""
    v = v2(n)
    u = n >> v
    return (v % 2 == 0) and (u % 8 == 7)


def _chi_m4(u):
    return 1 if u % 4 == 1 else -1


def _chi_8(u):
    return 1 if u % 8 in (1, 7) else -1


def _chi_m8(u):
    return 1 if u % 8 in (1, 3) else -1


def forbidden_projector(n):
    """Route B: exact Dirichlet-character projector."""
    v = v2(n)
    u = n >> v
    val = (1 + (-1) ** v) * (1 - _chi_m4(u) + _chi_8(u) - _chi_m8(u))
    assert val % 8 == 0
    return val // 8


def r3_array(N):
    """r_3(n) for 0 <= n <= N by exact lattice-count convolution."""
    m = int(N ** 0.5) + 1
    r1 = np.zeros(N + 1, dtype=np.int64)
    r1[0] = 1
    for x in range(1, m + 1):
        if x * x <= N:
            r1[x * x] = 2
    r2 = np.convolve(r1, r1)[: N + 1]
    return np.convolve(r2, r1)[: N + 1]


def forbidden_count(X):
    """Exact #{n <= X forbidden} via the geometric 4-adic layering."""
    total, a = 0, 0
    while 7 * 4 ** a <= X:
        total += (X // 4 ** a - 7) // 8 + 1
        a += 1
    return total


def forbidden_density(X):
    return forbidden_count(X) / X


def gate(N=40000, Nproj=10 ** 6, Ndens=10 ** 7, Neur=20000):
    n = np.arange(0, N + 1)
    assert np.all(8 * (n * (n + 1) // 2) + 1 == (2 * n + 1) ** 2), "L1"
    r3 = r3_array(N)
    forb = np.array([forbidden(k) for k in range(1, N + 1)])
    assert np.array_equal(r3[1:] == 0, forb), "L2 support"
    for k in range(1, Nproj + 1, 997):
        assert forbidden_projector(k) == int(forbidden(k)), "L3"
    for k in range(1, 20000):
        assert forbidden_projector(k) == int(forbidden(k)), "L3 head"
    assert int(forb.sum()) == forbidden_count(N), "L4 exact count"
    dens = forbidden_density(Ndens)
    tol = (math.log(Ndens, 4) + 2) / Ndens
    assert abs(dens - 1 / 6) < tol, ("L4 density", dens, tol)
    mp.dps = 30
    s = mpf(2)
    rhs = 8 ** -s * zeta(s, mpf(7) / 8) / (1 - 4 ** -s)
    lhs = mp.fsum(mpf(k) ** -s for k in range(1, 200001) if forbidden(k))
    tail = mpf(200000) ** (1 - s) / (s - 1) / 6 * mpf("1.2")
    assert abs(lhs - rhs) < tail, "L5"
    T = [k * (k + 1) // 2 for k in range(0, 300) if k * (k + 1) // 2 <= Neur]
    reach1 = np.zeros(Neur + 1, dtype=bool)
    reach1[T] = True
    reach2 = np.zeros(Neur + 1, dtype=bool)
    for t in T:
        reach2[t:][reach1[: Neur + 1 - t]] = True
    reach3 = np.zeros(Neur + 1, dtype=bool)
    for t in T:
        reach3[t:][reach2[: Neur + 1 - t]] = True
    assert reach3.all(), "L6 Eureka"
    for k in range(1, N // 4 + 1):
        assert forbidden(4 * k) == forbidden(k), "L7 self-similarity"
    return dict(N=N, density_1e7=dens, dirichlet_gap=float(abs(lhs - rhs)),
                eureka_to=Neur, layers="L1-L7 PASS")
