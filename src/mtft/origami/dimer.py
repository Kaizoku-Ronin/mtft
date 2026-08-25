"""mtft.origami.dimer — weighted planar bipartite graphs and their dimer ensembles.

General machinery (not tied to any one instance):

    DimerGraph          almost-perfect matchings, boundary measurement,
                        Plucker vector, Grassmann-Plucker certification,
                        coarse-graining by sufficient statistics
    ensemble_conservation   the gate that catches a dropped micro-state

Reference: P. Galashin, "Amplituhedra and Origami, I: Tree Level",
arXiv:2410.09574v2.  Boundary-measurement convention is (1.1)/Section 1.1:
  d(a) = { i : u_i black and used by a } U { i : u_i white and unused }.
"""
from __future__ import annotations

from itertools import combinations

import sympy as sp

__all__ = ["DimerGraph", "ensemble_conservation"]


class DimerGraph:
    """A weighted bipartite graph in a disk with degree-1 boundary vertices.

    Parameters
    ----------
    color : dict  vertex -> 'w' | 'b'
    edges : list of (u, v, weight)      weight may be a sympy expression
    boundary : list of vertices, in cyclic order (these are the u_i)

    All vertices not in ``boundary`` are interior.  An almost-perfect matching
    (APM) covers every interior vertex exactly once and every boundary vertex
    at most once.
    """

    def __init__(self, color, edges, boundary):
        self.color = dict(color)
        self.edges = [(u, v, sp.sympify(w)) for (u, v, w) in edges]
        self.boundary = list(boundary)
        self.interior = [v for v in self.color if v not in set(self.boundary)]
        self._inc = {v: [] for v in self.color}
        for k, (u, v, _w) in enumerate(self.edges):
            self._inc[u].append(k)
            self._inc[v].append(k)
        self._apms = None

    # ------------------------------------------------------------------
    @property
    def n(self):
        return len(self.boundary)

    def apms(self):
        """All APMs as a list of (frozenset of edge indices, weight)."""
        if self._apms is not None:
            return self._apms
        out = []
        interior = self.interior

        def rec(idx, covered, chosen, wt):
            while idx < len(interior) and interior[idx] in covered:
                idx += 1
            if idx == len(interior):
                out.append((frozenset(chosen), sp.expand(wt)))
                return
            v = interior[idx]
            for k in self._inc[v]:
                a, b, w = self.edges[k]
                o = b if a == v else a
                if o in covered:
                    continue
                rec(idx + 1, covered | {v, o}, chosen + [k], wt * w)

        rec(0, frozenset(), [], sp.Integer(1))
        self._apms = out
        return out

    def partition_function(self):
        return sp.expand(sum(w for _e, w in self.apms()))

    # ------------------------------------------------------------------
    def boundary_set(self, chosen):
        used = set()
        for k in chosen:
            a, b, _w = self.edges[k]
            for v in (a, b):
                if v in set(self.boundary):
                    used.add(self.boundary.index(v))
        return frozenset(
            i for i in range(self.n)
            if (self.color[self.boundary[i]] == "b" and i in used)
            or (self.color[self.boundary[i]] == "w" and i not in used)
        )

    def boundary_measurement(self):
        """dict  frozenset(I) -> Delta_I(Gamma, wt).  Keys are 0-based."""
        D = {}
        for chosen, w in self.apms():
            I = self.boundary_set(chosen)
            D[I] = sp.expand(D.get(I, 0) + w)
        return D

    def k(self):
        """Type (k, n); raises if the APMs do not all share a boundary size."""
        sizes = {len(I) for I in self.boundary_measurement()}
        if len(sizes) != 1:
            raise ValueError(f"not of a single type: boundary sizes {sizes}")
        return sizes.pop()

    def pluckers(self):
        """Full Plucker vector, zero-filled over all C(n, k) subsets."""
        D = self.boundary_measurement()
        kk = self.k()
        return {frozenset(I): D.get(frozenset(I), sp.Integer(0))
                for I in combinations(range(self.n), kk)}

    def verify_grassmann_plucker(self):
        """Certify that the Plucker vector is that of an actual point of Gr(k,n).

        Checks all three-term relations
            D_abc D_ade - D_abd D_ace + D_abe D_acd = 0
        (valid for k = 3).  Returns (n_checked, list_of_violations).
        """
        P = self.pluckers()
        if self.k() != 3:
            raise NotImplementedError("three-term check implemented for k=3")

        def d(*idx):
            return P[frozenset(idx)]

        bad, cnt = [], 0
        for a in range(self.n):
            rest = [x for x in range(self.n) if x != a]
            for b, c, e, f in combinations(rest, 4):
                cnt += 1
                val = sp.expand(d(a, b, c) * d(a, e, f)
                                - d(a, b, e) * d(a, c, f)
                                + d(a, b, f) * d(a, c, e))
                if val != 0:
                    bad.append(((a, b, c, e, f), val))
        return cnt, bad

    def cyclic_symmetry(self, shift):
        """True if the Plucker vector is invariant under I -> I + shift."""
        P = self.pluckers()
        return all(P[I] == P[frozenset((i + shift) % self.n for i in I)]
                   for I in P)

    # ------------------------------------------------------------------
    def sufficient_classes(self, varying):
        """Coarse-grain by occupation of the edges in ``varying``.

        ``varying`` is a list of edge indices.  Returns dict
        stat-tuple -> class coefficient, with the varying weights DIVIDED OUT,
        so that  Z = sum_k coeff_k * prod(w_e ** k_e).
        """
        cls = {}
        for chosen, w in self.apms():
            stat = tuple(1 if e in chosen else 0 for e in varying)
            coeff = w
            for e, bit in zip(varying, stat):
                if bit:
                    coeff = sp.cancel(coeff / self.edges[e][2])
            cls[stat] = sp.expand(cls.get(stat, 0) + coeff)
        return cls

    def class_partition_function(self, varying):
        """Z rebuilt from the coarse-grained classes (must equal Z exactly)."""
        cls = self.sufficient_classes(varying)
        tot = sp.Integer(0)
        for stat, coeff in cls.items():
            term = coeff
            for e, bit in zip(varying, stat):
                if bit:
                    term *= self.edges[e][2]
            tot += term
        return sp.expand(tot)


def ensemble_conservation(graph: DimerGraph, varying):
    """ENSEMBLE CONSERVATION gate.

    Aggregating micro-states into equal-sufficient-statistic classes must
    preserve the partition function EXACTLY:

        sum_{APM} wt(a)  ==  sum_{classes} coeff_k * theta^k

    This is the gate that catches a dropped matching (see CORRECTION RECORD:
    the (2,4) section-B class (0,0) was 1+q+s instead of (1+q)(1+s), the
    deficit being exactly the qs matching).  Raises AssertionError on failure;
    returns the class dictionary on success.
    """
    micro = graph.partition_function()
    coarse = graph.class_partition_function(varying)
    assert sp.expand(micro - coarse) == 0, (
        f"ensemble conservation FAILED: micro {micro} != coarse {coarse}, "
        f"deficit {sp.expand(micro - coarse)}"
    )
    return graph.sufficient_classes(varying)
