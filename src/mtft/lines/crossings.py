r"""mtft.lines.crossings — every crossing of two closed geodesics on X₀(N), exactly, and the linking numbers (LN-01, 2026-11-16).

No closed line of M = Γ₀(N)\PSL(2, ℝ) over the motions used here bounds by itself, so none has a
linking number with another.  A line γ̂ together with its reverse γ̂ʳ (the same geodesic with its
tangent reversed) does bound: K_γ = γ̂ ∪ γ̂ʳ has class (Ψ − Ψ; cls − cls) = 0.  For two such unions
(lab record MP13, derived before its run from the half-turn annulus over γ, in the orientation
dx ∧ dy ∧ dθ)

    Lk(K_γ, K_δ) = −I(γ, δ) ≤ 0,

I being the number of crossings of the two geodesics on the curve: pairs of branches through a
common point.  For a union L whose 26 periods vanish, Lk(K_γ, L) = −½ N(γ, L).  The formula is the
lab's derivation [derived]; the counts here are EXACT.

Counting.  Each closed geodesic of the modular surface is cut along the fundamental domain F into
exact arcs (``mtft.lines.cutting``).  At level N a closed line over the cycle of a dart x puts the
arc of frame m into the sheet x·m, one of the ψ(N) copies of F that tile the curve.  Two arcs in
one sheet cross where their semicircles do, at a point with rational x and y².  A crossing on the
boundary of F lies in two sheets (six at a corner), so one half-open convention counts it once:
the interior, the left edge and the left half of the bottom arc belong to the sheet; at a lift of
i the lower-numbered of the two sheets counts; at a lift of ρ all the branches through the point
are gathered in one frame and every pair is counted.  The mirrored convention (right edge, right
half, higher sheet, the corner ρ + 1) must give the same numbers: ``Crossings.check_conventions``.

``iota`` is the signed count, +1 when δ crosses γ from its right to its left; it is the
intersection form of the closed curve on the two classes, so it sees only the 26 periods, and most
of the linking is geometry that homology does not see.  Levels with elliptic points are refused.

Provenance: the Modular Plane lab, records MP13 (test T-G9) and MP14 (test T-G10), 12 November 2026.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from itertools import combinations

from . import closed as _c
from . import cutting as _k
from . import darts as _d

_INCLUDED = {"default": ("int", "left", "bot-"), "mirrored": ("int", "right", "bot+")}
_CORNER = {"default": "rho", "mirrored": "rho1"}


def _sgn(v) -> int:
    return (v > 0) - (v < 0)


class Crossings:
    """The closed geodesics over some motions at level N, with all their crossings.

    ``motions`` maps names to primitive hyperbolic matrices (the lab's three by default).  ``lines``
    lists the closed lines, motion by motion, each motion's by length and then by smallest dart;
    ``geodesics`` the distinct unoriented geodesics they trace, each oriented like its first line.
    """

    def __init__(self, N: int = 143, motions: dict | None = None):
        self.N = int(N)
        self.motions = {k: _d.as_matrix(A) for k, A in (motions if motions is not None else _c.MOTIONS).items()}
        G = _d.graph(self.N)                                            # refuses levels with elliptic points
        self.darts = D = G.darts
        self.cuts = {k: _k.walk(A) for k, A in self.motions.items()}
        self._perm = {(k, j): D.permutation(ar.m) for k, cut in self.cuts.items() for j, ar in enumerate(cut.arcs)}
        self.lines = [l for k, A in self.motions.items() for l in _c.lifts(A, self.N, k)]
        self.geodesics, self.geo_of_line, self.direction = [], [], []
        seen = {}
        for n, l in enumerate(self.lines):
            arcs = self.cuts[l.motion].arcs
            key = frozenset((self._perm[(l.motion, j)][x], _k.unoriented(ar.form)) for x in l.cycle for j, ar in enumerate(arcs))
            if key not in seen:
                seen[key] = len(self.geodesics)
                self.geodesics.append({"id": len(self.geodesics), "motion": l.motion, "ell": l.ell, "rep": n, "lines": [n]})
                self.direction.append(1)
            else:
                self.geodesics[seen[key]]["lines"].append(n)
                self.direction.append(-1)                               # a second line on a geodesic runs the other way
            self.geo_of_line.append(seen[key])
        self._line_at = {(l.motion, x): n for n, l in enumerate(self.lines) for x in l.cycle}
        self._xt = {}
        allarc = [(k, j, ar) for k, cut in self.cuts.items() for j, ar in enumerate(cut.arcs)]
        for (k1, j1, a1), (k2, j2, a2) in combinations(allarc, 2):
            p = _k.cross(a1.form, a2.form)
            if p is not None:
                sg = a1.s * a2.s * _sgn(a2.center - a1.center)
                self._xt[(k1, j1, k2, j2)] = (_k.kind(p), sg)
                self._xt[(k2, j2, k1, j1)] = (_k.kind(p), -sg)
        self._cache = {}

    # ------------------------------------------------------------------ the count
    def count(self, convention: str = "default") -> dict:
        """All crossings in one half-open convention: {'I', 'iota', 'places'}.

        I[g][h] counts the crossings of geodesics g and h (I[g][g] the self-crossings of g);
        iota[g][h] is the signed count; places counts crossings by where they sit on the tile
        ('int', 'left' or 'right', 'bot-' or 'bot+', 'i', 'rho' or 'rho1').
        """
        if convention in self._cache:
            return self._cache[convention]
        inc, corner = _INCLUDED[convention], _CORNER[convention]
        D, ng = self.darts, len(self.geodesics)
        Im = [[0] * ng for _ in range(ng)]
        Io = [[0] * ng for _ in range(ng)]
        places = {}
        by_sheet, inst = {}, []
        for g, G in enumerate(self.geodesics):
            l = self.lines[G["rep"]]
            for xp in l.cycle:
                for j in range(len(self.cuts[l.motion].arcs)):
                    it = (g, l.motion, j, self._perm[(l.motion, j)][xp])
                    inst.append(it)
                    by_sheet.setdefault(it[3], []).append(it)

        def add(g, h, sg, loc):
            places[loc] = places.get(loc, 0) + 1
            if g == h:
                Im[g][g] += 1
            else:
                Im[g][h] += 1
                Im[h][g] += 1
                Io[g][h] += sg
                Io[h][g] -= sg

        for s, lst in by_sheet.items():
            sS = D.pS[s]
            for p_, q_ in combinations(lst, 2):
                hit = self._xt.get((p_[1], p_[2], q_[1], q_[2]))
                if hit is None:
                    continue
                loc, sg = hit
                if loc == "i":
                    if (s < sS) if convention == "default" else (s > sS):
                        add(p_[0], q_[0], sg, loc)
                elif loc in inc:
                    add(p_[0], q_[0], sg, loc)
        groups = {}                                                     # the corner rule: all branches through a lift of rho
        for g, k, j, s in inst:
            ar = self.cuts[k].arcs[j]
            if corner not in (_k.kind(ar.entry), _k.kind(ar.exit)):
                continue
            Nm = ar.N
            if corner == "rho1":                                        # rho + 1 in the sheet s is rho in the sheet s T
                s, Nm = D.pT[s], _k.conj(_d.T, Nm)
            orb = [s, D.pR[s], D.pR2[s]]
            u = min(orb)
            r = orb.index(u)                                            # u = s R^r: the arc seen in u's frame is conj(R^r, N)
            if r:
                Nm = _k.conj(_d.R if r == 1 else _d.R2, Nm)
            a, b, c = _k.form(Nm)
            if a - Fr(b, 2) + c != 0:
                raise ArithmeticError("a branch gathered at a corner does not pass through it")
            groups.setdefault(u, []).append((g, 1 if a > 0 else -1, Fr(-b, 2 * a)))
        for lst in groups.values():
            for (g, sa, ca), (h, sb, cb) in combinations(lst, 2):
                add(g, h, sa * sb * _sgn(cb - ca), corner)
        out = {"I": Im, "iota": Io, "places": places, "sheets": len(by_sheet), "instances": len(inst), "corner_points": len(groups)}
        self._cache[convention] = out
        return out

    @property
    def I(self) -> list:
        return self.count()["I"]

    @property
    def iota(self) -> list:
        return self.count()["iota"]

    def check_conventions(self) -> bool:
        """The two half-open conventions give the same counts and the same signs."""
        a, b = self.count("default"), self.count("mirrored")
        return a["I"] == b["I"] and a["iota"] == b["iota"]

    # ------------------------------------------------------------------ what follows from the count
    def self_crossings(self) -> list:
        I = self.I
        return [I[g][g] for g in range(len(I))]

    def simple(self) -> list:
        """The geodesics that do not cross themselves."""
        return [g for g, n in enumerate(self.self_crossings()) if n == 0]

    def linking(self) -> list:
        """Lk(K_g, K_h) = −I(g, h) for g ≠ h; None on the diagonal."""
        I = self.I
        n = len(I)
        return [[(-I[g][h] if g != h else None) for h in range(n)] for g in range(n)]

    def family(self, motion: str) -> list:
        return [g for g, G in enumerate(self.geodesics) if G["motion"] == motion]

    def total(self, motion_a: str, motion_b: str) -> int:
        """The crossings between all the geodesics over one motion and all those over another."""
        I = self.I
        return sum(I[a][b] for a in self.family(motion_a) for b in self.family(motion_b))

    def label(self, g: int) -> str:
        G = self.geodesics[g]
        return f"{G['motion']}{G['ell']}(" + "-".join(str(n + 1) for n in G["lines"]) + ")"

    def geodesic_at(self, motion: str, dart) -> tuple:
        """(geodesic, direction) of the closed line of a motion through a dart (c, d): direction −1 if it runs against the geodesic's own."""
        n = self._line_at[(motion, self.darts.dart(*dart))]
        return self.geo_of_line[n], self.direction[n]

    def mirror(self) -> list:
        """For each geodesic: (its mirror image, the orientation factor), or (None, 0) if the image is not among them."""
        D, out = self.darts, []
        for G in self.geodesics:
            l = self.lines[G["rep"]]
            rel = _k.mirror_relation(l.A)
            if rel is None:
                out.append((None, 0))
                continue
            e, h = rel
            c, d = l.dart
            n = self._line_at[(l.motion, D.act(D.dart(c, -d), h))]
            out.append((self.geo_of_line[n], e * self.direction[n]))
        return out

    def reversed_by_mirror(self) -> list:
        """The geodesics that are their own mirror images with the direction reversed."""
        return [g for g, (img, e) in enumerate(self.mirror()) if img == g and e == -1]

    def mirror_check(self) -> bool:
        """The mirror permutes the geodesics, keeps every crossing count and reverses every signed count."""
        mir = self.mirror()
        if any(img is None for img, _ in mir):
            return False
        I, Io, n = self.I, self.iota, len(mir)
        img, eps = [m[0] for m in mir], [m[1] for m in mir]
        return sorted(img) == list(range(n)) and all(I[img[g]][img[g]] == I[g][g] for g in range(n)) and \
            all(I[img[g]][img[h]] == I[g][h] and eps[g] * eps[h] * Io[img[g]][img[h]] == -Io[g][h]
                for g in range(n) for h in range(n) if g != h)

    def summary(self) -> dict:
        c = self.count()
        n = len(self.geodesics)
        pairs = n * (n - 1) // 2
        return {"N": self.N, "lines": len(self.lines), "geodesics": n,
                "crossings_between": sum(c["I"][g][h] for g in range(n) for h in range(g + 1, n)),
                "self_crossings": sum(c["I"][g][g] for g in range(n)),
                "unlinked_pairs": sum(1 for g in range(n) for h in range(g + 1, n) if c["I"][g][h] == 0), "pairs": pairs,
                "signed_total": sum(abs(c["iota"][g][h]) for g in range(n) for h in range(g + 1, n)),
                "simple": len(self.simple()), "places": dict(sorted(c["places"].items()))}


def crossings_143() -> Crossings:
    """The lab's three motions at level 143: 28 lines on 24 geodesics."""
    return _cached(143, tuple(_c.MOTIONS.items()))


def crossings_143_extended() -> Crossings:
    """The six motions of the records MP13 and MP14 at level 143: 44 lines on 38 geodesics."""
    return _cached(143, tuple({**_c.MOTIONS, **_c.MORE_MOTIONS}.items()))


_STORE: dict = {}


def _cached(N: int, items: tuple) -> Crossings:
    key = (N, items)
    if key not in _STORE:
        _STORE[key] = Crossings(N, dict(items))
    return _STORE[key]


def crossings(N: int = 143, motions: dict | None = None) -> Crossings:
    """The crossings at level N over the given motions (cached by level and motions)."""
    m = motions if motions is not None else _c.MOTIONS
    return _cached(int(N), tuple((k, _d.as_matrix(A)) for k, A in m.items()))
