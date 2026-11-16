r"""mtft.lines.darts — the darts of P¹(ℤ/N), their graph, and the class of a group element (EXACT; LN-01, 2026-11-16).

Γ = Γ₀(N) has ψ(N) right cosets in SL(2, ℤ), and the coset of [[a, b], [c, d]] is the point
(c : d) of P¹(ℤ/N); a matrix g moves it by (c : d) ↦ (c : d) g.  At N = 143 there are 168.  The
points are called *darts* because three permutations make a graph of them, with

    S = [[0, −1], [1, 0]],   T = [[1, 1], [0, 1]],   R = S⁻¹T = [[0, 1], [−1, −1]]   (T = S R, R³ = 1):

    vertex = orbit of R (three darts),   edge = orbit of S (two darts),   face = orbit of T.

A face is a cusp of the curve and its length is the cusp's width.  The graph is the dual of the
Farey tessellation folded by Γ; it is a deformation retract of the punctured curve Γ\ℍ, so its
cycles are the first homology of the punctured curve, a free group of rank b₁ = E − V + 1.  At
N = 143: 56 vertices, 84 edges (one a loop), 4 faces of lengths 1, 11, 13, 143, b₁ = 29, genus 13.

The class of an element.  A matrix is a word in S and T (``word``).  Read from left to right the
word moves a dart, every letter S crossing the dart's edge and every letter T crossing the edge
and turning at the vertex.  ``Graph.cls`` adds the signed crossings; an element of Γ brings the
dart (0 : 1) back, and its 29 crossings of the edges outside a spanning tree are its class in
H₁(Γ\ℍ, ℤ) = ℤ²⁹.  The class is additive and unchanged by conjugation in Γ.  Each edge outside the
tree gives one generator of Γ (``Graph.generators``), whose class is a unit vector.

The construction (order of the points, breadth-first tree from the vertex of (0 : 1)) is the lab
record MP12's, so the generators are the same matrices in the same order; ``mtft.lines.table``
checks that against the frozen list.  Levels with elliptic points (S or R fixing a dart) have no
such graph here: ``Graph`` refuses them.  ``Darts`` works at every level.
"""
from __future__ import annotations

import math
from functools import lru_cache

I2 = ((1, 0), (0, 1))
S = ((0, -1), (1, 0))
T = ((1, 1), (0, 1))
TI = ((1, -1), (0, 1))
R = ((0, 1), (-1, -1))            # R = S^-1 T, R^3 = 1, T = S R


def mm(A, B):
    """The product of two 2 x 2 matrices given as tuples of rows."""
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))


def inv(A):
    """The inverse of a matrix of determinant 1."""
    return ((A[1][1], -A[0][1]), (-A[1][0], A[0][0]))


def neg(A):
    return ((-A[0][0], -A[0][1]), (-A[1][0], -A[1][1]))


def mpow(A, n: int):
    """A to the power n (n may be negative for a matrix of determinant 1)."""
    if n < 0:
        A, n = inv(A), -n
    out, B = I2, A
    while n:
        if n & 1:
            out = mm(out, B)
        B = mm(B, B)
        n >>= 1
    return out


R2 = mm(R, R)


def as_matrix(g):
    """A matrix as a tuple of two tuples of Python integers; refuses determinants other than 1."""
    (a, b), (c, d) = g
    a, b, c, d = int(a), int(b), int(c), int(d)
    if a * d - b * c != 1:
        raise ValueError("the matrix must have integer entries and determinant 1")
    return ((a, b), (c, d))


def trace(A) -> int:
    return A[0][0] + A[1][1]


def is_hyperbolic(A) -> bool:
    return abs(trace(A)) > 2


def length(A) -> float:
    """The length 2 arccosh(|tr A| / 2) of the closed geodesic of a hyperbolic matrix."""
    t = abs(trace(A))
    if t <= 2:
        raise ValueError("the matrix is not hyperbolic")
    return 2.0 * math.acosh(t / 2.0)


def word(g) -> list:
    """g = T^q0 S T^q1 S ... T^qn up to sign, as tokens ('T', q) and ('S',)."""
    (a, b), (c, d) = g
    out = []
    while c != 0:
        q = a // c
        a, b = a - q * c, b - q * d
        out.append(("T", q))
        out.append(("S",))
        a, b, c, d = c, d, -a, -b
    out.append(("T", b * d))
    return out


def unword(tokens):
    """The matrix a list of tokens spells (the inverse of ``word`` up to sign)."""
    g = I2
    for tok in tokens:
        g = mm(g, S if tok[0] == "S" else mpow(T, tok[1]))
    return g


def chi6(g) -> int:
    """The character of SL(2, ℤ) to ℤ/6 with χ(T) = 1 and χ(S) = 3, read off the word."""
    return sum(tok[1] if tok[0] == "T" else 3 for tok in word(g)) % 6


def orbits(perm) -> list:
    """The cycles of a permutation given as a list, each from its smallest element."""
    seen, out = [False] * len(perm), []
    for i in range(len(perm)):
        if not seen[i]:
            o, j = [], i
            while not seen[j]:
                seen[j] = True
                o.append(j)
                j = perm[j]
            out.append(o)
    return out


class Darts:
    """P¹(ℤ/N): the right cosets of Γ₀(N) in SL(2, ℤ), with the action of matrices on them.

    ``points`` lists the canonical pairs (c, d) in increasing order (the smallest pair of each
    class under the units of ℤ/N), ``index`` inverts it and ``base`` is the dart (0 : 1), the
    coset of Γ₀(N) itself.
    """

    def __init__(self, N: int = 143):
        N = int(N)
        if N < 1:
            raise ValueError("the level must be a positive integer")
        self.N = N
        units = [u for u in range(N) if math.gcd(u, N) == 1]
        canon = {}
        for c in range(N):
            for d in range(N):
                if (c, d) in canon or math.gcd(math.gcd(c, d), N) != 1:
                    continue
                for u in units:                       # (c, d) is the smallest of its class
                    canon[((u * c) % N, (u * d) % N)] = (c, d)
        self._canon = canon
        self.points = sorted(set(canon.values()))
        self.index = {p: i for i, p in enumerate(self.points)}
        self.size = len(self.points)
        self.base = self.index[canon[(0, 1 % N)]]
        self.pS = self.permutation(S)
        self.pT = self.permutation(T)
        self.pR = self.permutation(R)
        self.pR2 = self.permutation(R2)

    def dart(self, c: int, d: int) -> int:
        """The index of the point (c : d)."""
        N = self.N
        key = (c % N, d % N)
        if key not in self._canon:
            raise ValueError(f"({c} : {d}) is not a point of the projective line modulo {N}")
        return self.index[self._canon[key]]

    def act(self, i: int, g) -> int:
        """The dart (c : d) g."""
        c, d = self.points[i]
        return self.dart(c * g[0][0] + d * g[1][0], c * g[0][1] + d * g[1][1])

    def permutation(self, g) -> list:
        return [self.act(i, g) for i in range(self.size)]

    def cycles(self, g) -> list:
        """The cycles of g on the darts, sorted by length and then by smallest dart."""
        return sorted(orbits(self.permutation(g)), key=lambda o: (len(o), min(o)))

    def cycle_type(self, g) -> list:
        """The lengths of the cycles of g, in increasing order."""
        return sorted(len(o) for o in orbits(self.permutation(g)))

    def in_group(self, g) -> bool:
        """Whether g lies in Γ₀(N)."""
        return g[1][0] % self.N == 0

    def elliptic(self) -> tuple:
        """The numbers of darts fixed by S and by R (elliptic points of orders 2 and 3)."""
        return (sum(1 for i, j in enumerate(self.pS) if i == j), sum(1 for i, j in enumerate(self.pR) if i == j))


class Graph:
    """The graph on the darts of P¹(ℤ/N) and the class map of Γ₀(N), for levels without elliptic points."""

    def __init__(self, N: int = 143):
        D = darts(N)
        if D.elliptic() != (0, 0):
            raise ValueError(f"level {N} has elliptic points: the graph on its darts is not built here")
        self.darts, self.N = D, D.N
        pS, pR, pR2, pT = D.pS, D.pR, D.pR2, D.pT
        self.vertices, self.edges, self.faces = orbits(pR), orbits(pS), orbits(pT)
        self.vof = {x: k for k, o in enumerate(self.vertices) for x in o}
        self.eof = {x: k for k, o in enumerate(self.edges) for x in o}
        self.sgn = {x: (1 if x == min(o) else -1) for o in self.edges for x in o}
        self.fof = {x: k for k, o in enumerate(self.faces) for x in o}
        self.fpos = {x: j for o in self.faces for j, x in enumerate(o)}
        self.loops = [k for k, o in enumerate(self.edges) if self.vof[o[0]] == self.vof[o[1]]]
        x0 = D.base
        tree, seen, queue = set(), {self.vof[x0]}, [self.vof[x0]]       # breadth-first spanning tree
        while queue:
            v = queue.pop(0)
            for x in self.vertices[v]:
                w = self.vof[pS[x]]
                if w not in seen:
                    seen.add(w)
                    tree.add(self.eof[x])
                    queue.append(w)
        if len(seen) != len(self.vertices):
            raise ValueError("the graph is not connected")
        self.tree = tree
        self.nontree = [k for k in range(len(self.edges)) if k not in tree]
        self._ntpos = {k: j for j, k in enumerate(self.nontree)}
        self.b1 = len(self.nontree)
        nv, ne, nf = len(self.vertices), len(self.edges), len(self.faces)
        self.genus = (2 - (nv - ne + nf)) // 2
        rep = {x0: I2}                                                  # coset representatives along the tree
        queue = [x0]
        while queue:
            x = queue.pop(0)
            for y, g in ((pR[x], R), (pR2[x], R2)) + (((pS[x], S),) if self.eof[x] in tree else ()):
                if y not in rep:
                    rep[y] = mm(rep[x], g)
                    queue.append(y)
        self.rep = rep
        gens = []
        for k in self.nontree:
            tail = self.tail(k)
            gens.append(mm(mm(rep[tail], S), inv(rep[pS[tail]])))
        self.generators = gens

    # ------------------------------------------------------------------ the pieces of the graph
    def tail(self, k: int) -> int:
        """The dart of edge k that carries the sign +1."""
        o = self.edges[k]
        return o[0] if self.sgn[o[0]] == 1 else o[1]

    def summary(self) -> dict:
        return {"N": self.N, "darts": self.darts.size, "V": len(self.vertices), "E": len(self.edges), "F": len(self.faces),
                "loops": len(self.loops), "faces": sorted(len(o) for o in self.faces), "b1": self.b1, "genus": self.genus}

    def cusp_of_face(self, f: int) -> int:
        """The denominator gcd(c, N) of the cusp a face stands for (N for the cusp ∞, 1 for the cusp 0)."""
        c, _ = self.darts.points[self.faces[f][0]]
        return math.gcd(c, self.N)

    def cusp_loop(self, f: int):
        """The element of Γ₀(N) that goes once round the cusp of face f (a conjugate of T^width)."""
        o = self.faces[f]
        x = min(o)
        return mm(mm(self.rep[x], mpow(T, len(o))), inv(self.rep[x]))

    def face_cycle(self, f: int) -> tuple:
        """The class of the boundary of face f, straight from its darts."""
        vec = [0] * len(self.edges)
        for y in self.faces[f]:
            vec[self.eof[y]] += self.sgn[y]
        return tuple(vec[k] for k in self.nontree)

    # ------------------------------------------------------------------ paths and classes
    def trace_path(self, tokens, x: int, vec: list) -> int:
        """Move the dart x along the tokens, adding the signed edge crossings to vec; return the last dart."""
        D, eof, sgn = self.darts, self.eof, self.sgn
        for tok in tokens:
            if tok[0] == "S":
                vec[eof[x]] += sgn[x]
                x = D.pS[x]
                continue
            k = tok[1]
            if k == 0:
                continue
            o = self.faces[self.fof[x]]
            w, p = len(o), self.fpos[x]
            full, rem = divmod(abs(k), w)
            if full:
                turns = full if k > 0 else -full
                for y in o:
                    vec[eof[y]] += turns * sgn[y]
            if k > 0:
                for j in range(rem):
                    y = o[(p + j) % w]
                    vec[eof[y]] += sgn[y]
            else:
                for j in range(1, rem + 1):
                    y = o[(p - j) % w]
                    vec[eof[y]] -= sgn[y]
            x = o[(p + k) % w]
        return x

    def crossings(self, g, start: int | None = None) -> list:
        """The signed crossings of all the edges by the path of g from a dart (the base dart by default)."""
        vec = [0] * len(self.edges)
        x = self.darts.base if start is None else start
        end = self.trace_path(word(g), x, vec)
        if end != x:
            raise ValueError("not a closed path: the element does not fix the starting dart")
        return vec

    def cls(self, g, start: int | None = None) -> tuple:
        """The class of g in ℤ^b₁: its signed crossings of the edges outside the spanning tree."""
        vec = self.crossings(g, start)
        return tuple(vec[k] for k in self.nontree)

    def reduce(self, vec) -> tuple:
        """The b₁ coordinates of a full vector of edge crossings."""
        return tuple(vec[k] for k in self.nontree)

    def is_cycle(self, vec) -> bool:
        """Whether a vector of edge crossings is a cycle: at every vertex as many exits as entries."""
        bal = [0] * len(self.vertices)
        pS = self.darts.pS
        for k in range(len(self.edges)):
            tail = self.tail(k)
            bal[self.vof[tail]] -= vec[k]
            bal[self.vof[pS[tail]]] += vec[k]
        return all(v == 0 for v in bal)

    def element(self, c) -> tuple:
        """An element of Γ₀(N) with the class c: the product of the generators to those powers."""
        if len(c) != self.b1:
            raise ValueError(f"a class has {self.b1} coordinates")
        g = I2
        for e, n in enumerate(c):
            if n:
                g = mm(g, mpow(self.generators[e], int(n)))
        return g

    def chi_form(self) -> list:
        """χ on the generators: χ(γ) = Σ cls_e χ_e modulo 6 for every γ in Γ₀(N)."""
        return [chi6(g) for g in self.generators]

    def random_element(self, rng, length: int = 6):
        """A random word of the given length in the generators and their inverses (``rng``: a ``random.Random``)."""
        g = I2
        for _ in range(length):
            h = rng.choice(self.generators)
            g = mm(g, h if rng.random() < 0.5 else inv(h))
        return g


@lru_cache(maxsize=8)
def darts(N: int = 143) -> Darts:
    """The darts of level N (cached)."""
    return Darts(N)


@lru_cache(maxsize=8)
def graph(N: int = 143) -> Graph:
    """The graph of level N (cached); raises ValueError at levels with elliptic points."""
    return Graph(N)
