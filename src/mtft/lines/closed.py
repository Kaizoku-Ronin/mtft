r"""mtft.lines.closed — the closed lines of Γ₀(N)\PSL(2, ℝ) over a motion, and their integers (EXACT; LN-01, 2026-11-16).

The 3-fold M = Γ₀(N)\PSL(2, ℝ) is the unit tangent bundle of the curve, and the geodesic flow on
it is the frozen frame of the lab's instrument 12.7.  A hyperbolic matrix A of SL(2, ℤ) (a
*motion*) gives one closed orbit of the flow at level 1.  At level N it permutes the darts of
P¹(ℤ/N), and every cycle of length ℓ through the coset of g gives one closed line of M over it,

    γ = g A^ℓ g⁻¹ ∈ Γ₀(N),      of length ℓ · 2 arccosh(|tr A| / 2).

Its class in H₁(M) ⊗ ℚ = ℚ^(b₁ + 1) is (Ψ; cls): Rademacher's Ψ along the fibre and the class of γ
in H₁ of the punctured curve.  ``lifts`` builds the lines of one motion; each ``Line`` carries ℓ,
its darts, γ, Ψ, cls and (at level 143) the three windings.  Facts the tests gate:

    Ψ(γ) = ℓ Ψ(A);      Ψ ≡ χ(cls) modulo 6;      Σ over the lifts of A of cls = 0;
    Σ over the lifts of A of Ψ = ψ(N) Ψ(A)        (ψ(143) = 168).

So no union of all the lifts bounds unless Ψ(A) = 0, and no single line of the 28 below bounds;
a line together with its reverse always does (``mtft.lines.crossings``).

``MOTIONS`` holds the lab's three: the cat map [[2, 1], [1, 1]] (cycles 7⁴ 35⁴, Ψ = 0), the motion
of Sphereland 10, [[84, −47], [143, −80]] (1⁴ 5⁴ 6⁴ 30⁴, Ψ = −1; trace 4, so its four lines of length
1 are systoles of X₀(143)), and the Pell unit [[12, 1], [143, 12]] (1, 11, 13, 143, Ψ = −21).
``lines_143`` returns their 28 lines in the order of the lab record MP12.

Provenance: the Modular Plane lab, record MP12 (test T-G4, 12 November 2026).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from . import darts as _d
from . import rademacher as _r

MOTIONS = {"cat": ((2, 1), (1, 1)), "s10": ((84, -47), (143, -80)), "pell": ((12, 1), (143, 12))}
"""The three motions of the lab record MP12."""

MORE_MOTIONS = {"uv3": ((4, 1), (3, 1)), "u2v2": ((5, 2), (2, 1)), "uv5": ((6, 1), (5, 1))}
"""The three motions added in the record MP14: UV³, U²V², UV⁵ (U = T, V = [[1, 0], [1, 1]])."""


@dataclass(frozen=True)
class Line:
    """One closed line of Γ₀(N)\\PSL(2, ℝ): a cycle of a motion on the darts, with its integers."""
    motion: str
    A: tuple
    N: int
    ell: int
    dart: tuple            # the smallest dart (c, d) of the cycle
    cycle: tuple           # the darts' indices, from the smallest, in the order A moves them
    gamma: tuple           # g A^ell g⁻¹, an element of Γ₀(N)
    psi: int
    cls: tuple
    chi: int
    windings: tuple | None

    @property
    def length(self) -> float:
        return self.ell * _d.length(self.A)

    @property
    def trace(self) -> int:
        return abs(_d.trace(self.gamma))

    def describe(self) -> dict:
        return {"motion": self.motion, "ell": self.ell, "dart": list(self.dart), "psi": self.psi,
                "phi": None if self.windings is None else list(self.windings), "cls": list(self.cls)}


def _windings_or_none(g, N: int):
    return _r.windings(g) if N == 143 else None


@lru_cache(maxsize=64)
def _lifts(A, N: int, name: str) -> tuple:
    G = _d.graph(N)
    D = G.darts
    wA = _d.word(A)
    psiA = _r.Psi(A)
    out = []
    for o in D.cycles(A):
        x, ell = o[0], len(o)                     # o starts from its smallest dart
        vec = [0] * len(G.edges)
        y = x
        for _ in range(ell):
            y = G.trace_path(wA, y, vec)
        if y != x:
            raise ArithmeticError("the cycle did not close")
        gam = _d.mm(_d.mm(G.rep[x], _d.mpow(A, ell)), _d.inv(G.rep[x]))
        if not D.in_group(gam):
            raise ArithmeticError("the lift is not in the group")
        c = G.reduce(vec)
        psi = _r.Psi(gam)
        if psi != ell * psiA:
            raise ArithmeticError("Ψ of a lift is not ℓ times Ψ of the motion")
        out.append(Line(name, A, N, ell, D.points[x], tuple(o), gam, psi, c, _d.chi6(gam), _windings_or_none(gam, N)))
    return tuple(out)


def lifts(A, N: int = 143, name: str | None = None) -> tuple:
    """The closed lines over the hyperbolic matrix A at level N, sorted by length and then by smallest dart."""
    A = _d.as_matrix(A)
    if not _d.is_hyperbolic(A):
        raise ValueError("a motion must be hyperbolic")
    return _lifts(A, int(N), name if name is not None else f"{A[0][0]},{A[0][1]};{A[1][0]},{A[1][1]}")


def lines_143() -> tuple:
    """The 28 closed lines of the lab record MP12, in its order (numbered from 1 there)."""
    return tuple(l for name, A in MOTIONS.items() for l in lifts(A, 143, name))


def invariants(gamma, N: int = 143) -> dict:
    """Ψ, the class, χ and (at level 143) the windings of an element of Γ₀(N) given directly."""
    gamma = _d.as_matrix(gamma)
    G = _d.graph(N)
    if not G.darts.in_group(gamma):
        raise ValueError(f"the matrix is not in the group of level {N}")
    out = {"trace": abs(_d.trace(gamma)), "psi": _r.Psi(gamma), "cls": G.cls(gamma), "chi": _d.chi6(gamma),
           "windings": _windings_or_none(gamma, N)}
    out["length"] = _d.length(gamma) if _d.is_hyperbolic(gamma) else 0.0
    return out


def transfer(A, N: int = 143) -> dict:
    """The sums over all the lifts of a motion: the classes add up to zero, the fibre numbers to ψ(N) Ψ(A)."""
    L = lifts(A, N)
    b1 = len(L[0].cls)
    out = {"lifts": len(L), "cycle_type": sorted(l.ell for l in L), "psi": _r.Psi(_d.as_matrix(A)),
           "sum_cls": tuple(sum(l.cls[j] for l in L) for j in range(b1)), "sum_psi": sum(l.psi for l in L),
           "darts": _d.darts(N).size}
    if N == 143:
        out["sum_windings"] = tuple(sum(l.windings[i] for l in L) for i in range(3))
    return out


def partner(line: Line) -> tuple:
    """The cycle of darts that S carries the line's cycle to (for a motion S conjugates to its inverse, the reverse line)."""
    pS = _d.darts(line.N).pS
    return tuple(sorted(pS[x] for x in line.cycle))
