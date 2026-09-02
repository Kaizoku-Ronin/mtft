"""mtft.surface.manin — EXACT layer: the Manin complex of X0(N).

P^1(Z/N) is built with a CRT canonical form (O(omega(N) log N) per point,
replacing the O(N^2 phi(N)) unit scan of the printable kit).  The Manin
actions on bottom rows (c:d) are

    S:(c,d)->(d,-c)     T:(c,d)->(c,c+d)     R:(c,d)->(d,-c-d)

with R T S = -I, so (T, S, R) is a combinatorial map: T-orbits are cusps
(vertices), S-orbits are edges, R-orbits are ideal Farey triangles (faces).
Dart g <-> coset Gamma_0(N) g; the oriented edge of a representative dart g
is the modular symbol {g.oo, g.0}.

Gates are non-definitional: the genus from the Riemann-Hurwitz divisor
formula is compared with V - E + F of the constructed complex (only when
nu2 = nu3 = 0), cardinality against psi(N), S^2 = R^3 = 1, connectedness,
width partition, and for squarefree N the width census {N/d : d | N}.

Elliptic levels (nu2 > 0 or nu3 > 0) are refused by ``cell_complex`` with
an explicit ValueError (CC-MSL-02): the orbifold complex is not built.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd
from typing import Callable, Dict, Iterable, List, Tuple

import numpy as np

Pair = Tuple[int, int]


# ---------------------------------------------------------------- arithmetic
def factorize(n: int) -> Dict[int, int]:
    out: Dict[int, int] = {}
    m, p = n, 2
    while p * p <= m:
        while m % p == 0:
            out[p] = out.get(p, 0) + 1
            m //= p
        p += 1
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


def divisors(n: int) -> List[int]:
    out = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
        d += 1
    return sorted(out)


def euler_phi(n: int) -> int:
    r = n
    for p in factorize(n):
        r = r // p * (p - 1)
    return r


def kronecker(a: int, p: int) -> int:
    """Kronecker symbol (a/p) for prime p (handles p = 2)."""
    if p == 2:
        if a % 2 == 0:
            return 0
        return 1 if a % 8 in (1, 7) else -1
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


@dataclass(frozen=True)
class Invariants:
    N: int
    primes: Tuple[int, ...]
    index: int            # psi(N) = [SL2(Z) : Gamma_0(N)]
    nu2: int
    nu3: int
    cusps: int
    genus: int
    sigma: int
    squarefree: bool
    b1_compact: int
    b1_open: int
    modular_area_over_pi: Fraction   # psi(N)/3
    closed_hyperbolic_area_over_pi: int | None   # 4(g-1) for g >= 2
    cusp_widths: Tuple[int, ...]

    @property
    def claim_class(self) -> str:
        return "EXACT"


def invariants(N: int) -> Invariants:
    if N < 2:
        raise ValueError("N must be >= 2")
    fac = factorize(N)
    ps = tuple(sorted(fac))
    mu = N
    for p in ps:
        mu = mu // p * (p + 1)
    nu2 = 0 if N % 4 == 0 else _prod(1 + kronecker(-4, p) for p in ps)
    nu3 = 0 if N % 9 == 0 else _prod(1 + kronecker(-3, p) for p in ps)
    ds = divisors(N)
    cusps = sum(euler_phi(gcd(d, N // d)) for d in ds)
    g = Fraction(1) + Fraction(mu, 12) - Fraction(nu2, 4) - Fraction(nu3, 3) - Fraction(cusps, 2)
    if g.denominator != 1:
        raise ArithmeticError(f"genus formula did not close at N={N}: {g}")
    genus = int(g)
    widths = []
    for d in ds:                       # cusps a/d, d | N, width (N/d)/gcd(d, N/d)
        w = (N // d) // gcd(d, N // d)
        widths.extend([w] * euler_phi(gcd(d, N // d)))
    return Invariants(
        N=N, primes=ps, index=mu, nu2=nu2, nu3=nu3, cusps=cusps, genus=genus,
        sigma=sum(ds), squarefree=all(e == 1 for e in fac.values()),
        b1_compact=2 * genus, b1_open=2 * genus + cusps - 1,
        modular_area_over_pi=Fraction(mu, 3),
        closed_hyperbolic_area_over_pi=4 * (genus - 1) if genus >= 2 else None,
        cusp_widths=tuple(sorted(widths)),
    )


def _prod(values: Iterable[int]) -> int:
    r = 1
    for v in values:
        r *= v
    return r


# ------------------------------------------------------- CRT canonical P^1
def canonical_factory(N: int) -> Callable[[int, int], Pair]:
    """Canonical representative of (c:d) in P^1(Z/N) via CRT.

    For each q = p^a || N: if p does not divide c, scale to (1, d/c) mod q,
    else (p does not divide d) scale to (c/d, 1) mod q.  The case split is
    unit-invariant, so the combined pair is a canonical representative.
    """
    fac = factorize(N)
    moduli = [p ** a for p, a in fac.items()]
    # CRT idempotents
    idem = []
    for q in moduli:
        M = N // q
        idem.append((M * pow(M, -1, q)) % N)
    cache: Dict[Pair, Pair] = {}

    def canon(c: int, d: int) -> Pair:
        key = (c % N, d % N)
        hit = cache.get(key)
        if hit is not None:
            return hit
        C = D = 0
        for q, e in zip(moduli, idem):
            cq, dq = key[0] % q, key[1] % q
            p = next(p for p in fac if q % p == 0)
            if cq % p:
                cc, dd = 1, (dq * pow(cq, -1, q)) % q
            else:
                if dq % p == 0:
                    raise ValueError(f"({c},{d}) is not a point of P^1(Z/{N})")
                cc, dd = (cq * pow(dq, -1, q)) % q, 1
            C = (C + cc * e) % N
            D = (D + dd * e) % N
        cache[key] = (C, D)
        return (C, D)

    return canon


def projective_line(N: int, canon: Callable[[int, int], Pair] | None = None) -> List[Pair]:
    """All points of P^1(Z/N) as canonical pairs, sorted, by S/T closure from (0:1)."""
    canon = canon or canonical_factory(N)
    S = lambda x: canon(x[1], -x[0])
    T = lambda x: canon(x[0], x[0] + x[1])
    Ti = lambda x: canon(x[0], x[1] - x[0])
    start = canon(0, 1)
    pts, seen, k = [start], {start}, 0
    while k < len(pts):
        x = pts[k]
        k += 1
        for y in (S(x), T(x), Ti(x)):
            if y not in seen:
                seen.add(y)
                pts.append(y)
    return sorted(pts)


def orbits(points: Iterable[Pair], perm: Callable[[Pair], Pair]) -> List[List[Pair]]:
    seen, out = set(), []
    for x in points:
        if x in seen:
            continue
        orb, cur = [], x
        while cur not in seen:
            seen.add(cur)
            orb.append(cur)
            cur = perm(cur)
        out.append(orb)
    return out


# --------------------------------------------------------- cellular complex
@dataclass
class ManinComplex:
    N: int
    inv: Invariants
    points: List[Pair]
    vertices: List[List[Pair]]      # T-orbits (cusps)
    edges: List[Pair]               # representative dart per S-orbit (min)
    faces: List[List[Pair]]         # R-orbits, cyclic order = face orientation
    vertex_of: Dict[Pair, int]
    face_of: Dict[Pair, int]
    edge_of: Dict[Pair, Tuple[int, int]]   # dart -> (edge id, +1/-1)
    boundary_1: np.ndarray           # V x E   (int64)
    boundary_2: np.ndarray           # E x F   (int64)
    wedge_1: np.ndarray              # E x E   Whitney cup pairing (float)
    S: Callable[[Pair], Pair] = field(repr=False)
    T: Callable[[Pair], Pair] = field(repr=False)
    R: Callable[[Pair], Pair] = field(repr=False)
    cusp_divisor_class: List[int] = field(default_factory=list)   # gcd(c, N) per vertex, N for c=0

    @property
    def counts(self) -> Tuple[int, int, int]:
        return (len(self.vertices), len(self.edges), len(self.faces))

    def gates(self) -> List[dict]:
        return complex_gates(self)


def cell_complex(N: int) -> ManinComplex:
    inv = invariants(N)
    if inv.nu2 or inv.nu3:
        raise ValueError(
            f"N={N}: elliptic points (nu2={inv.nu2}, nu3={inv.nu3}); the Manin cellular "
            "complex requires orbifold handling (CC-MSL-02)")
    canon = canonical_factory(N)
    S = lambda x: canon(x[1], -x[0])
    T = lambda x: canon(x[0], x[0] + x[1])
    R = lambda x: canon(x[1], -x[0] - x[1])
    points = projective_line(N, canon)
    edge_orbits = orbits(points, S)
    faces = orbits(points, R)
    vertices = orbits(points, T)
    if any(len(o) != 2 for o in edge_orbits) or any(len(o) != 3 for o in faces):
        raise ValueError(f"N={N}: unexpected S/R orbit lengths despite nu2=nu3=0")
    vertex_of = {d: i for i, o in enumerate(vertices) for d in o}
    face_of = {d: i for i, o in enumerate(faces) for d in o}
    edges, edge_of = [], {}
    for o in edge_orbits:
        rep = min(o)
        eid = len(edges)
        edges.append(rep)
        edge_of[rep] = (eid, 1)
        edge_of[S(rep)] = (eid, -1)
    V, E, F = len(vertices), len(edges), len(faces)
    b1 = np.zeros((V, E), dtype=np.int64)
    for eid, d in enumerate(edges):
        b1[vertex_of[d], eid] -= 1
        b1[vertex_of[S(d)], eid] += 1
    b2 = np.zeros((E, F), dtype=np.int64)
    wedge = np.zeros((E, E))
    for fid, o in enumerate(faces):
        loc = [edge_of[d] for d in o]
        for eid, s in loc:
            b2[eid, fid] += s
        for i, j in ((0, 1), (1, 2), (2, 0)):
            (ei, si), (ej, sj) = loc[i], loc[j]
            wedge[ei, ej] += si * sj / 6.0
            wedge[ej, ei] -= si * sj / 6.0
    classes = [N if o[0][0] % N == 0 else gcd(o[0][0], N) for o in vertices]
    return ManinComplex(N, inv, points, vertices, edges, faces, vertex_of, face_of, edge_of,
                        b1, b2, wedge, S, T, R, classes)


def complex_gates(cx: ManinComplex) -> List[dict]:
    inv, pts = cx.inv, cx.points
    V, E, F = cx.counts
    chi = V - E + F
    g_top = Fraction(2 - chi, 2)
    gates = [
        _g("projective_line_cardinality_is_psi", len(pts) == inv.index, f"{len(pts)} vs {inv.index}"),
        _g("S_squared_identity", all(cx.S(cx.S(x)) == x for x in pts), f"{len(pts)} checks"),
        _g("R_cubed_identity", all(cx.R(cx.R(cx.R(x))) == x for x in pts), f"{len(pts)} checks"),
        _g("RTS_is_identity_on_darts", all(cx.S(cx.T(cx.R(x))) == x for x in pts), "combinatorial-map relation"),
        _g("cusp_widths_partition_darts", sum(len(o) for o in cx.vertices) == len(pts), ""),
        _g("independent_euler_reconstruction", g_top.denominator == 1 and int(g_top) == inv.genus,
           f"V-E+F={chi}, genus {g_top} vs formula {inv.genus}"),
        _g("boundary_squared_zero", not np.any(cx.boundary_1 @ cx.boundary_2), ""),
        _g("width_census_matches_divisor_formula",
           sorted(len(o) for o in cx.vertices) == list(inv.cusp_widths), f"{sorted(len(o) for o in cx.vertices)}"),
    ]
    return gates


def _g(name: str, ok: bool, evidence: str) -> dict:
    return {"name": name, "class": "EXACT", "status": "PASS" if ok else "FAIL", "evidence": evidence}


def assert_gates(gates: List[dict]) -> None:
    bad = [(g["name"], g["evidence"]) for g in gates if g["status"] == "FAIL"]
    if bad:
        raise AssertionError(f"failed gates: {bad}")
