r"""mtft.lines.walls — walls between twists: the wall distance, the systole band and one-way bottleneck costs.

EXACT within a radius; WL-01, 2026-11-20.  The shortest closed geodesic in each class mod 2, the exact wall ball and the wall
distance between twists.

A wall in 4D across which the twist changes is a vortex wrapping a curve on the closed curve X₀(143) (FF-17).  A curve of
class Y changes ζ by the twist dual to Y, so the cheapest wall from ζ₁ to ζ₂ has length

    N(y) = least total length of closed geodesics whose classes add to y mod 2,   y = C(ζ₁ ⊕ ζ₂),

with C the cup form of ``twists``.  d_W(ζ₁, ζ₂) = N(C(ζ₁ ⊕ ζ₂)) is a translation-invariant metric on the twists (FF-20).

Layers:

  level1_classes        primitive hyperbolic classes of PSL(2, ℤ) up to a trace, as canonical words in R and L (FF-14)
  primitive_classes     their lifts to X₀(N): every primitive oriented closed geodesic up to a length, with its class (FF-14)
  shortest_by_class     the shortest geodesic in each nonzero class mod 2
  wall_ball             every class within a radius, with its least wall and the geodesics that make it (a Dijkstra on the
                        classes reached).  Exact: a wall of length ≤ r uses geodesics ≤ r only, all of them enumerated
  three_sums            the second route: every sum of at most three distinct classes, enumerated directly
  wall_distance         d_W between two twists within a radius (None means: farther than the radius, certified)
  systoles, systole_band  the four systoles of X₀(143) (2 arccosh 2), and how many of them a twist makes antiperiodic.
                        FF-20 found the twist energy fixed to within a few percent by this count (five bands that do not
                        overlap over all 2²⁶ twists); that energy is not in the package, so the bands' depths are the record's
  thin_wall_exponent    Coleman's B = 27π²σ⁴/(2ε³) for a release ε > 0; infinite otherwise
  bottleneck            the least, over routes, of a route's hardest step: c(x, z) ≤ max(c(x, y), c(y, z)), a
                        quasi-ultrametric (not symmetric when the step costs are one-way)

Records: FF-14 (the enumeration), FF-17 (walls), FF-20 (the ball sizes 1, 5, 7, 33, 289 and 1,020 at radii 2.5, 2.75,
5, 5.5, 8 and 10; the bands; the bottleneck), in the MTFT Python project under claude/experiments/.  Wave 8 compared
``wall_ball(16)`` with FF-20's own programme class by class: the same 120,827 classes, N equal to the record's float32.
"""
from __future__ import annotations

import heapq
import math
from functools import lru_cache

from . import darts as _d
from . import twists as _tw

SYSTOLE = 2 * math.acosh(2)


# ------------------------------------------------------------------ level-1 classes and their lifts (FF-14)
def _block(a: int, b: int):
    return ((1 + a * b, a), (b, 1))


def _mul(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))


def _canonical(blocks) -> bool:
    """The block sequence is its own smallest rotation and is not periodic (a primitive cyclic word)."""
    for s in range(1, len(blocks)):
        rot = blocks[s:] + blocks[:s]
        if rot <= blocks:
            return False
    return True


def level1_classes(tmax: int) -> list:
    """Primitive hyperbolic classes of PSL(2, ℤ) with trace ≤ tmax, as (blocks of R^a L^b, the matrix)."""
    out = []

    def tr(M):
        return M[0][0] + M[1][1]

    def dfs(blocks, P):
        if blocks and _canonical(blocks):
            out.append((tuple(blocks), P))
        a0 = blocks[0][0] if blocks else 1
        a = a0
        while True:
            b = blocks[0][1] if (blocks and a == a0) else 1
            if tr(_mul(P, _block(a, b))) > tmax:
                if a == a0 and b > 1:
                    a += 1
                    continue
                break
            while True:
                Pb = _mul(P, _block(a, b))
                if tr(Pb) > tmax:
                    break
                blocks.append((a, b))
                dfs(blocks, Pb)
                blocks.pop()
                b += 1
            a += 1

    dfs([], ((1, 0), (0, 1)))
    return out


def lifts(A, lmax: float, N: int = 143) -> list:
    """The primitive closed geodesics of X₀(N) over the level-1 class A with length ≤ lmax: (length, ell, cls)."""
    G = _tw.graph(N)
    DA = G.darts
    LA = _d.length(A)
    perm = DA.permutation(A)
    wA = _d.word(A)
    kmax = int(lmax / LA + 1e-9)
    seen, res = [False] * DA.size, []
    for i in range(DA.size):
        if seen[i]:
            continue
        o, j = [], i
        while not seen[j]:
            seen[j] = True
            o.append(j)
            j = perm[j]
        ell = len(o)
        if ell > kmax:
            continue
        x = min(o)
        vec = [0] * len(G.edges)
        y = x
        for _ in range(ell):
            y = G.trace_path(wA, y, vec)
        if y != x:
            raise ArithmeticError("a cycle did not close")
        res.append((ell * LA, ell, tuple(G.reduce(vec))))
    return res


@lru_cache(maxsize=8)
def primitive_classes(lmax: float, N: int = 143) -> tuple:
    """Every primitive oriented closed geodesic of X₀(N) with length ≤ lmax: (length, ell, cls), by length.

    A geodesic and its reverse are both listed (they lie over a level-1 class and its inverse); mod 2 their classes agree.
    At N = 143 the eight shortest, of length 2 arccosh 2, are the four systoles in both orientations."""
    tmax = math.floor(2 * math.cosh(lmax / 2) + 1e-9)
    out = []
    for blocks, A in level1_classes(tmax):
        if _d.length(A) <= lmax + 1e-12:
            out.extend(lifts([[A[0][0], A[0][1]], [A[1][0], A[1][1]]], lmax, N))
    out.sort(key=lambda r: r[0])
    return tuple(out)


@lru_cache(maxsize=8)
def shortest_by_class(lmax: float, N: int = 143) -> tuple:
    """(Y, length) for the shortest geodesic in each nonzero class mod 2 with length ≤ lmax, by length."""
    seen = {}
    for L, ell, cls in primitive_classes(lmax, N):
        y = _tw.class_bits(cls, N)
        if y and y not in seen:
            seen[y] = L
    return tuple(sorted(seen.items(), key=lambda kv: kv[1]))


# ------------------------------------------------------------------ the wall ball, its second route, the distance
@lru_cache(maxsize=8)
def wall_ball(radius: float, N: int = 143) -> dict:
    """{Y: (N(Y), the lengths of the geodesics of a least wall)} for every class with N(Y) ≤ radius (Y = 0 included).

    N(Y) is exact.  The geodesics are one least wall.  Lengths over one level-1 class are commensurable (ℓ L_A), so a class
    can have least walls with different numbers of geodesics: within 16, 64 classes do, all through such sums (for
    example 5L = 2L + 3L over the systole's level-1 class, L = 2 arccosh 2).  Which one is returned depends on the order
    of the search."""
    gens = shortest_by_class(radius, N)
    best = {0: (0.0, -1, 0)}
    heap = [(0.0, 0)]
    done = set()
    while heap:
        c, y = heapq.heappop(heap)
        if y in done:
            continue
        done.add(y)
        rem = radius - c
        for i, (gy, gl) in enumerate(gens):
            if gl > rem + 1e-12:
                break
            z, nc = y ^ gy, c + gl
            if z not in done and (z not in best or nc < best[z][0] - 1e-12):
                best[z] = (nc, i, y)
                heapq.heappush(heap, (nc, z))
    out = {}
    for y in done:
        pieces, x = [], y
        while x:
            c, i, p = best[x]
            pieces.append(gens[i][1])
            x = p
        out[y] = (best[y][0], tuple(sorted(pieces)))
    return out


def ball_sizes(radii, radius: float | None = None, N: int = 143) -> dict:
    """|B(r)|: twists within wall distance r of a given twist (the same from every twist)."""
    ball = wall_ball(max(radii) if radius is None else radius, N)
    return {r: sum(1 for c, _ in ball.values() if c <= r + 1e-12) for r in radii}


def three_sums(radius: float, N: int = 143) -> dict:
    """{Y: least total length of at most three distinct classes adding to Y, if ≤ radius} (the second route; Y ≠ 0)."""
    gens = shortest_by_class(radius, N)
    R = {}

    def put(y, v):
        if y and v < R.get(y, math.inf):
            R[y] = v

    n = len(gens)
    for i in range(n):
        yi, li = gens[i]
        if li > radius:
            break
        put(yi, li)
        for j in range(i + 1, n):
            lj = li + gens[j][1]
            if lj > radius:
                break
            yj = yi ^ gens[j][0]
            put(yj, lj)
            for k in range(j + 1, n):
                lk = lj + gens[k][1]
                if lk > radius:
                    break
                put(yj ^ gens[k][0], lk)
    return R


def compare_routes(radius: float, N: int = 143) -> dict:
    """The ball against the direct sums: the sums never beat the ball, and agree wherever the ball's wall has ≤ 3 geodesics."""
    ball, R = wall_ball(radius, N), three_sums(radius, N)
    beats = sum(1 for y, v in R.items() if v < ball.get(y, (math.inf,))[0] - 1e-9)
    small = [(y, c) for y, (c, p) in ball.items() if 1 <= len(p) <= 3]
    disagree = sum(1 for y, c in small if abs(R.get(y, math.inf) - c) > 1e-9)
    return {"classes in the ball": len(ball), "walls of 1-3 geodesics": len(small), "sums that beat the ball": beats,
            "walls of 1-3 geodesics where the routes disagree": disagree}


def wall_distance(z1: int, z2: int, radius: float, N: int = 143):
    """d_W(ζ₁, ζ₂) if it is ≤ radius, else None (then the distance is certified to exceed the radius)."""
    hit = wall_ball(radius, N).get(_tw.dual_class(z1 ^ z2, N))
    return None if hit is None else hit[0]


def wall(z1: int, z2: int, radius: float, N: int = 143):
    """The geodesic lengths of a least wall from ζ₁ to ζ₂ within the radius, else None."""
    hit = wall_ball(radius, N).get(_tw.dual_class(z1 ^ z2, N))
    return None if hit is None else hit[1]


# ------------------------------------------------------------------ systoles and bands
@lru_cache(maxsize=4)
def systoles(N: int = 143) -> tuple:
    """The classes (26 bits) of the systoles of X₀(N): geodesics of length 2 arccosh 2 at N = 143."""
    gens = shortest_by_class(SYSTOLE + 0.5, N)
    lmin = gens[0][1]
    return tuple(y for y, l in gens if abs(l - lmin) < 1e-9)


def systole_band(z: int, N: int = 143) -> int:
    """How many systoles the twist makes antiperiodic (0 to 4 at N = 143)."""
    return sum(_tw.parity(z & y) for y in systoles(N))


def check_systoles(N: int = 143) -> dict:
    """The systoles' classes, their rank mod 2, their pairwise intersections, and the sizes of the bands they make."""
    S = systoles(N)
    pair = [[_tw.intersection(a, b, N) for b in S] for a in S]
    r = _tw.rank_f2(S)
    from math import comb
    return {"systole classes": len(S), "rank mod 2": r, "pairwise intersections": pair,
            "band sizes over all twists": {k: comb(r, k) * 2 ** (_tw.dim(N) - r) for k in range(r + 1)}}


# ------------------------------------------------------------------ one-way costs
def thin_wall_exponent(sigma: float, eps: float) -> float:
    """Coleman's thin-wall B = 27π²σ⁴/(2ε³) for a release ε > 0; infinite uphill and between degenerate vacua."""
    return math.inf if eps <= 0 else 27 * math.pi ** 2 * sigma ** 4 / (2 * eps ** 3)


def bottleneck(cost):
    """Least over routes of the hardest step, from a square table of one-step costs (math.inf for no step).

    The result c satisfies c[x][z] ≤ max(c[x][y], c[y][z]): a route costs its hardest step, because the time a chain of
    tunnellings takes is ruled by its largest exponent."""
    n = len(cost)
    c = [[0.0 if i == j else cost[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            cik = c[i][k]
            if cik == math.inf:
                continue
            for j in range(n):
                c[i][j] = min(c[i][j], max(cik, c[k][j]))
    return c
