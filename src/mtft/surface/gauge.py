"""mtft.surface.gauge — gauge theory ON the modular surface (reading A).

Layers and tags
---------------
EXACT / STANDARD  (theorems of gauge theory on a closed Riemann surface,
                   evaluated at the exact invariants of X0(N)):
  * U(1) flux quantisation: (1/2pi) int F = m in Z; on the closed
    curvature -1 surface of area A = 4 pi (g-1) the minimal action in sector
    m is S_m = 2 pi^2 m^2 / (e^2 A); the flux sum is a theta series.
  * flat U(1) connections: H^1(X, R)/H^1(X, Z) = U(1)^{2g}; on the punctured
    curve b_1 = 2g + n_cusp - 1.
  * 2D Yang-Mills (Migdal/Witten): Z = sum_R (dim R)^{2-2g} exp(-g^2 A C_2(R)/2)
    for U(1), SU(2), SU(3), SU(N); with n boundary holonomies U_r the
    character-insertion form sum_R (dim R)^{2-2g-n} e^{...} prod chi_R(U_r).
  * Riemann-Roch on a line bundle of degree d: h^0 - h^1 = d + 1 - g.
  * Aharony-Seiberg-Tachikawa line-operator census for su(N): Lagrangian
    subgroups of (Z/N)^2 with the Dirac pairing; sigma(N) of them, grouped by
    the gauge group SU(N)/Z_k (k | N) with k discrete theta variants; the
    2 pi theta shift acts by the Witten effect (e, m) -> (e + m, m) and its
    orbits have size k / gcd(k, N/k).
    THEOREM (gated here by two independent routes): for squarefree N the
    cyclic Lagrangian subgroups are exactly P^1(Z/N), the Witten shift is the
    Manin T action, and the theta-orbit sizes are the cusp widths of
    Gamma_0(N).  Cusps of X0(N) are theta-orbit / line-operator sectors —
    not "stable particles" (the printable-kit wording is retired).

PHYSICS_OVERLAY (interpretation, not derived): reading X0(N) as a coupling
space of a 4D theory (reading B), BPS charge overlays, electroweak skeletons,
the M4 x X0(N) lift.  Nothing in this module asserts those.

Maass eigenmodes of the hyperbolic Laplacian are spectral dynamics on the
surface, not electromagnetic waves; the Whitney spectra of :mod:`hodge` are
not Maass eigenvalues (different metric in the same conformal class).
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import gcd, pi
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from .manin import Invariants, canonical_factory, divisors, invariants, orbits, projective_line


# ------------------------------------------------------------- U(1) Maxwell
def closed_area(inv: Invariants) -> float:
    if inv.genus < 2:
        raise ValueError("closed constant-curvature -1 metric needs g >= 2")
    return 4 * pi * (inv.genus - 1)


def flux_action(inv: Invariants, m: int, e2: float = 1.0) -> float:
    """S_m = 2 pi^2 m^2 / (e^2 A), the minimal Maxwell action in flux sector m (EXACT formula)."""
    return 2 * pi ** 2 * m * m / (e2 * closed_area(inv))


def flux_partition_sum(inv: Invariants, e2: float = 1.0, mmax: int = 200) -> float:
    """sum_m exp(-S_m): a theta series theta_3(0, q) with q = exp(-pi/(2 e^2 (g-1)))."""
    return float(sum(np.exp(-flux_action(inv, m, e2)) for m in range(-mmax, mmax + 1)))


def flat_connection_torus_dimension(inv: Invariants, punctured: bool = False) -> int:
    return inv.b1_open if punctured else inv.b1_compact


# ---------------------------------------------------------- 2D Yang-Mills
@dataclass(frozen=True)
class Irrep:
    label: str
    dim: int
    casimir: float           # C_2(R) with the normalisation stated in ``group``


def irreps(group: str, cutoff: int) -> List[Irrep]:
    """Irreps with C_2 below ``cutoff`` (SU(2): j(j+1); SU(3): (p^2+q^2+pq+3p+3q)/3; U(1): n^2)."""
    out: List[Irrep] = []
    if group == "U(1)":
        n = 0
        while n * n <= cutoff:
            out.append(Irrep(f"n={n}", 1, n * n))
            if n:
                out.append(Irrep(f"n=-{n}", 1, n * n))
            n += 1
    elif group == "SU(2)":
        two_j = 0
        while (two_j / 2) * (two_j / 2 + 1) <= cutoff:
            j = two_j / 2
            out.append(Irrep(f"j={j:g}", two_j + 1, j * (j + 1)))
            two_j += 1
    elif group == "SU(3)":
        for p in range(0, 60):
            for q in range(0, 60):
                c2 = (p * p + q * q + p * q + 3 * p + 3 * q) / 3
                if c2 <= cutoff:
                    out.append(Irrep(f"({p},{q})", (p + 1) * (q + 1) * (p + q + 2) // 2, c2))
    else:
        raise ValueError("group must be U(1), SU(2) or SU(3)")
    return sorted(out, key=lambda r: (r.casimir, r.dim, r.label))


def ym_partition_function(inv: Invariants, group: str = "SU(2)", g_ym2: float = 1.0,
                          area: float | None = None, cutoff: float = 60.0) -> Dict:
    """Z = sum_R (dim R)^{2-2g} exp(-g_YM^2 A C_2(R)/2) on the closed surface (EXACT series, truncated)."""
    A = closed_area(inv) if area is None else area
    chi = 2 - 2 * inv.genus
    terms = []
    for R in irreps(group, cutoff):
        w = R.dim ** chi * np.exp(-g_ym2 * A * R.casimir / 2)
        terms.append((R.label, R.dim, R.casimir, float(w)))
    Z = sum(t[3] for t in terms)
    tail = terms[-1][3] / Z if terms else 0.0
    return {"class": "EXACT_SERIES", "group": group, "genus": inv.genus, "area": A, "g_ym2": g_ym2,
            "Z": float(Z), "terms": terms, "last_term_relative": float(tail),
            "converged": tail < 1e-12}


def su2_character(j: float, theta: float) -> float:
    """chi_j(U) for U with eigenvalues e^{+-i theta/2}."""
    if abs(np.sin(theta / 2)) < 1e-12:
        return 2 * j + 1
    return float(np.sin((2 * j + 1) * theta / 2) / np.sin(theta / 2))


def ym_partition_function_with_cusps(inv: Invariants, holonomies: Sequence[float], g_ym2: float = 1.0,
                                     area: float | None = None, cutoff: float = 60.0) -> Dict:
    """SU(2) with n boundary circles (truncated cusps) carrying holonomy angles theta_r:
    Z = sum_j (2j+1)^{2-2g-n} exp(-g^2 A j(j+1)/2) prod_r chi_j(theta_r).  EXACT series."""
    A = closed_area(inv) if area is None else area
    n = len(holonomies)
    chi = 2 - 2 * inv.genus - n
    terms = []
    for R in irreps("SU(2)", cutoff):
        j = (R.dim - 1) / 2
        w = R.dim ** chi * np.exp(-g_ym2 * A * R.casimir / 2) * np.prod([su2_character(j, t) for t in holonomies])
        terms.append((R.label, float(w)))
    Z = sum(t[1] for t in terms)
    return {"class": "EXACT_SERIES", "group": "SU(2)", "cusps": n, "Z": float(Z), "terms": terms}


def same_genus_control(N1: int, N2: int, g_ym2: float = 0.01) -> Dict:
    """Closed 2D YM cannot separate same-genus levels; cusp holonomies can."""
    i1, i2 = invariants(N1), invariants(N2)
    if i1.genus != i2.genus:
        raise ValueError("control requires equal genus")
    z1 = ym_partition_function(i1, "SU(2)", g_ym2)["Z"]
    z2 = ym_partition_function(i2, "SU(2)", g_ym2)["Z"]
    demo = pi / 3
    c1 = ym_partition_function_with_cusps(i1, [demo] * i1.cusps, g_ym2)["Z"]
    c2 = ym_partition_function_with_cusps(i2, [demo] * i2.cusps, g_ym2)["Z"]
    rel = abs(c1 - c2) / max(abs(c1), abs(c2), 1e-300)
    return {"closed_equal": abs(z1 - z2) < 1e-12 * max(abs(z1), 1), "closed": (z1, z2),
            "cusps": (i1.cusps, i2.cusps), "with_cusp_holonomy": (c1, c2),
            "cusp_relative_difference": rel,
            "note": "nontrivial reps carry (dim R)^{2-2g}: at g=13 the separation is O(2^-24) relative"}


# ------------------------------------------------------------ Riemann-Roch
def line_bundle_index(inv: Invariants, degree: int) -> int:
    """h^0(L) - h^1(L) = deg L + 1 - g (EXACT)."""
    return degree + 1 - inv.genus


# ------------------------------------------- AST line-operator census
def lagrangian_subgroups(N: int) -> List[frozenset]:
    """All order-N subgroups of (Z/N)^2 isotropic for <(e,m),(e',m')> = e m' - m e' mod N.
    Brute force (N <= ~60 comfortably): route B for the census, independent of Manin symbols."""
    pts = [(e, m) for e in range(N) for m in range(N)]

    def span(gens):
        S = {(0, 0)}
        frontier = [(0, 0)]
        while frontier:
            x = frontier.pop()
            for g_ in gens:
                y = ((x[0] + g_[0]) % N, (x[1] + g_[1]) % N)
                if y not in S:
                    S.add(y)
                    frontier.append(y)
        return frozenset(S)

    found = set()
    for a in pts:
        for b in pts:
            if (a[0] * b[1] - a[1] * b[0]) % N:
                continue
            H = span([a, b])
            if len(H) == N and all((x[0] * y[1] - x[1] * y[0]) % N == 0 for x in H for y in H):
                found.add(H)
    return sorted(found, key=lambda H: sorted(H))


def gauge_group_of(H: frozenset, N: int) -> int:
    """k such that the electric sublattice H cap (Z/N x 0) = k Z/N: the group SU(N)/Z_k."""
    electric = sorted(e for e, m in H if m == 0)
    k = N
    for e in electric:
        k = gcd(k, e)
    return k if k else N


def witten_shift(H: frozenset, N: int) -> frozenset:
    return frozenset(((e + m) % N, m) for e, m in H)


def line_operator_census(N: int) -> Dict:
    """AST census and the theorem gate against Gamma_0(N) cusp data (two independent routes)."""
    inv = invariants(N)
    Ls = lagrangian_subgroups(N)
    by_k: Dict[int, List[frozenset]] = {}
    for H in Ls:
        by_k.setdefault(gauge_group_of(H, N), []).append(H)
    # theta orbits under the Witten shift
    seen, orbit_sizes = set(), []
    for H in Ls:
        if H in seen:
            continue
        size, cur = 0, H
        while cur not in seen:
            seen.add(cur)
            size += 1
            cur = witten_shift(cur, N)
        orbit_sizes.append(size)
    cyclic = [H for H in Ls if any(len({(i * e % N, i * m % N) for i in range(N)}) == N for e, m in H)]
    gates = [
        {"name": "count_is_sigma_N", "class": "EXACT", "status": "PASS" if len(Ls) == inv.sigma else "FAIL",
         "evidence": f"{len(Ls)} vs sigma={inv.sigma}"},
        {"name": "variants_per_gauge_group_k", "class": "EXACT",
         "status": "PASS" if all(len(by_k.get(k, [])) == k for k in divisors(N)) else "FAIL",
         "evidence": {k: len(v) for k, v in sorted(by_k.items())}},
        {"name": "cyclic_lagrangians_equal_P1_cardinality", "class": "EXACT",
         "status": "PASS" if len(cyclic) == inv.index else "FAIL", "evidence": f"{len(cyclic)} vs psi={inv.index}"},
    ]
    if inv.squarefree:
        gates.append({"name": "theta_orbit_sizes_equal_cusp_widths", "class": "EXACT",
                      "status": "PASS" if sorted(orbit_sizes) == list(inv.cusp_widths) else "FAIL",
                      "evidence": f"orbits {sorted(orbit_sizes)} vs widths {list(inv.cusp_widths)}"})
    return {"N": N, "sigma": inv.sigma, "lagrangian_count": len(Ls), "by_gauge_group": {k: len(v) for k, v in sorted(by_k.items())},
            "theta_orbit_sizes": sorted(orbit_sizes), "cusp_widths": list(inv.cusp_widths), "gates": gates,
            "dictionary": {f"SU({N})/Z_{k}": f"cusp class d={N // k}, width {(k // gcd(k, N // k))}" for k in divisors(N)},
            "class": "EXACT (theorem for squarefree N; PHYSICS_OVERLAY only in the naming)"}


def manin_theta_orbits(N: int) -> List[int]:
    """Route A for the same census on P^1(Z/N): sizes of Manin T-orbits (cusp widths)."""
    canon = canonical_factory(N)
    T = lambda x: canon(x[0], x[0] + x[1])
    return sorted(len(o) for o in orbits(projective_line(N, canon), T))
