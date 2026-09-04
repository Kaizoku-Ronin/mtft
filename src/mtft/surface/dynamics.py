"""mtft.surface.dynamics — linear Hamiltonian dynamics on H_1(X0(143), R) with genericity controls.

Port of the Modular Surface Laboratory Waves 11-13 core (Sol, 2026-09-02):
Hamiltonian flows x' = F x with F = J A, A self-adjoint for the Hodge metric
G = E J, and the greedy real Lie closure with the absolute two-tier residual
gate.  Audit (Claude, 2026-09-03): full closure dim 351 = dim sp(26, R) is the
generic outcome for any pair of Hamiltonians with inter-block coupling
(Kuranishi genericity), so it carries no arithmetic information by itself.
``closure_controls`` therefore ships three mandatory controls:

  random        two random G-self-adjoint Hamiltonians      -> 351 (generic)
  hecke_commut  J P_block for the Hecke block projectors     ->   4 (abelian)
  block_diag    two random block-diagonal Hamiltonians       -> 127 = dim sp(2)+sp(4)+sp(8)+sp(12)

The informative quantity is the distance from the 127-dimensional block
subalgebra, reported as ``off_block_fraction`` of each generator.  A closure
dimension that random generators also reach is not evidence (falsifier F8).

Frames: ``stage_from_periods`` uses the mtft.periods symplectic frame (as the
lab did); ``stage_from_frozen`` uses the surface cycle frame.  Do not mix.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Sequence

import numpy as np

SP26 = 13 * 27


class AmbiguousClosure(RuntimeError):
    pass


@dataclass
class Stage:
    E: np.ndarray            # symplectic (intersection) form on cycles
    J: np.ndarray            # complex structure, J^2 = -I
    G: np.ndarray            # Hodge metric E J (or -E J), SPD
    projectors: Dict[str, np.ndarray]
    frame: str

    def orthonormal(self) -> "Stage":
        """Conjugate by G^{1/2} so the Hodge metric is the identity; Frobenius norms then measure
        Hodge-metric size and the closure residual gate is frame-independent."""
        w, U = np.linalg.eigh((self.G + self.G.T) / 2)
        R = U @ np.diag(np.sqrt(w)) @ U.T
        Ri = U @ np.diag(1 / np.sqrt(w)) @ U.T
        conj = lambda A: R @ A @ Ri
        return Stage(Ri @ self.E @ Ri, conj(self.J), np.eye(len(w)),
                     {k: conj(P) for k, P in self.projectors.items()}, self.frame + "|G-orthonormal")

    def gates(self) -> Dict[str, float]:
        n = self.E.shape[0]
        rel = lambda a, b: float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))
        return {"J_squared": rel(self.J @ self.J, -np.eye(n)),
                "G_symmetric": rel(self.G, self.G.T),
                "G_min_eigenvalue": float(np.linalg.eigvalsh((self.G + self.G.T) / 2)[0]),
                "projectors_resolve_identity": rel(sum(self.projectors.values()), np.eye(n)) if self.projectors else 0.0}


def hecke_block_projectors(T: np.ndarray, G: np.ndarray, charpoly_factors: Sequence) -> Dict[str, np.ndarray]:
    """G-orthogonal projectors onto the generalized eigenspaces of the integer Hecke matrix T,
    one per irreducible factor (sympy Polys) of its characteristic polynomial."""
    import sympy as sp
    x = sp.Symbol("x")
    Tf = T.astype(float)
    out = {}
    for label, f in charpoly_factors:
        others = np.eye(T.shape[0])
        for lab2, f2 in charpoly_factors:
            if lab2 != label:
                coeffs = [float(c) for c in sp.Poly(f2, x).all_coeffs()]
                M = np.zeros_like(Tf)
                for c in coeffs:
                    M = M @ Tf + c * np.eye(T.shape[0])
                others = others @ M
        # image of 'others' is the block; G-orthogonal projector onto it
        U, s, _ = np.linalg.svd(others)
        r = int(np.sum(s > 1e-8 * s[0]))
        B = U[:, :r]
        out[label] = B @ np.linalg.solve(B.T @ G @ B, B.T @ G)
    return out


def stage_from_frozen() -> Stage:
    from .frozen import x0143
    from .hodge_structure import _rel  # noqa: F401
    import sympy as sp
    d = x0143()
    E = d["intersection_cycles"].astype(float)
    J = d["J_true"]
    G = E @ J
    G = (G + G.T) / 2
    if np.linalg.eigvalsh(G)[0] < 0:
        G = -G
    x = sp.Symbol("x")
    factors = [("ell", x**2), ("ghost", (x + 2)**4), ("q4", (x**4 - 3*x**3 - x**2 + 5*x + 1)**2),
               ("q6", (x**6 - 10*x**4 + 2*x**3 + 24*x**2 - 7*x - 12)**2)]
    proj = hecke_block_projectors(d["T2"], G, factors)
    return Stage(E, J, G, proj, "surface_cycle_frame")


def stage_from_periods(dps: int = 30) -> Stage:
    from mtft.periods.involutions import transported_intersection
    from mtft.periods.physics import hodge_metric_hecke, hodge_structure_hecke
    G = np.asarray(hodge_metric_hecke(dps), dtype=float)
    G = (G + G.T) / 2
    J = np.asarray(hodge_structure_hecke(dps), dtype=float)
    E = np.asarray(transported_intersection(), dtype=float)
    return Stage(E, J, G, {}, "periods_symplectic_frame")


def hamiltonian(stage: Stage, energy: np.ndarray) -> np.ndarray:
    """F = J A with A = G^-1 S the G-self-adjoint operator of a symmetric energy matrix S."""
    S = (energy + energy.T) / 2
    return stage.J @ np.linalg.solve(stage.G, S)


def flow_gates(stage: Stage, F: np.ndarray, times: Iterable[float] = (0.5, 1.0, 2.0)) -> Dict[str, float]:
    from scipy.linalg import expm
    E = stage.E
    rel = lambda a, b: float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))
    out = {"hamiltonian_residual": float(np.linalg.norm(F.T @ E + E @ F) / np.linalg.norm(F.T @ E))}
    for t in times:
        U = expm(t * F)
        out[f"symplectic_t{t}"] = rel(U.T @ E @ U, E)
        out[f"reversible_t{t}"] = rel(expm(-t * F) @ U, np.eye(E.shape[0]))
    return out


def symplectic_projector(E: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    return lambda B: (B - np.linalg.solve(E, B.T @ E)) / 2


def lie_closure(generators: Iterable[np.ndarray], E: np.ndarray, tol_hi: float = 1e-6, tol_lo: float = 1e-9,
                target_dim: int = SP26, projector=None) -> Dict:
    """Greedy real Lie closure with an absolute two-tier residual gate (after Sol, Wave 11)."""
    n = E.shape[0]
    projector = projector or symplectic_projector(E)
    basis: List[np.ndarray] = []
    rows = np.empty((0, n * n))
    depths: List[int] = []
    acc, rej = [], []

    def attempt(C: np.ndarray, depth: int) -> None:
        nonlocal rows
        w = C.reshape(-1).astype(float).copy()
        if np.linalg.norm(w) < 1e-13:
            rej.append(float(np.linalg.norm(w)))
            return
        for _ in range(5):
            w = projector(w.reshape(n, n)).reshape(-1)
            if len(rows):
                w -= rows.T @ (rows @ w)
        r = float(np.linalg.norm(w))
        if tol_lo < r <= tol_hi:
            raise AmbiguousClosure(f"residual {r:.3e} in ({tol_lo:g}, {tol_hi:g}]")
        if r <= tol_lo:
            rej.append(r)
            return
        w /= r
        basis.append(w.reshape(n, n))
        rows = np.vstack((rows, w))
        depths.append(depth)
        acc.append(r)

    for F in generators:
        attempt(np.asarray(F), 0)
    i = 0
    while i < len(basis) and len(basis) < target_dim:
        for j in range(i):
            attempt(basis[i] @ basis[j] - basis[j] @ basis[i], max(depths[i], depths[j]) + 1)
            if len(basis) >= target_dim:
                break
        i += 1
    growth = [sum(x <= d for x in depths) for d in sorted(set(depths))]
    return {"dimension": len(basis), "target": target_dim, "growth": growth,
            "min_accepted": min(acc) if acc else None, "max_rejected": max(rej) if rej else 0.0,
            "separation": (min(acc) / max(rej)) if acc and rej and max(rej) else None}


def off_block_fraction(stage: Stage, F: np.ndarray) -> float:
    """Norm fraction of F outside the Hecke block-diagonal subalgebra (0 for block-diagonal)."""
    P = stage.projectors
    if not P:
        raise ValueError("stage has no Hecke projectors")
    diag = sum(Pk @ F @ Pk for Pk in P.values())
    return float(np.linalg.norm(F - diag) / np.linalg.norm(F))


def _closure_dim(gens, E, bands=((1e-6, 1e-9), (1e-5, 1e-7), (1e-4, 1e-6))):
    """Closure dimension with the first unambiguous band; records the band used."""
    for hi, lo in bands:
        try:
            r = lie_closure(gens, E, tol_hi=hi, tol_lo=lo)
            return {"dimension": r["dimension"], "band": (hi, lo), "separation": r["separation"]}
        except AmbiguousClosure as exc:
            last = str(exc)
    return {"dimension": None, "band": None, "verdict": f"AMBIGUOUS: {last}"}


def closure_controls(stage: Stage, seed: int = 143, orthonormal: bool = True) -> Dict:
    if orthonormal and not np.allclose(stage.G, np.eye(stage.G.shape[0])):
        stage = stage.orthonormal()
    rng = np.random.default_rng(seed)
    n = stage.E.shape[0]
    rand = lambda: hamiltonian(stage, rng.standard_normal((n, n)))
    out = {"random_pair": _closure_dim([rand(), rand()], stage.E)}
    if stage.projectors:
        P = stage.projectors
        out["hecke_commuting"] = _closure_dim([stage.J @ Pk for Pk in P.values()], stage.E)

        def block():
            A = np.linalg.solve(stage.G, (lambda S: (S + S.T) / 2)(rng.standard_normal((n, n))))
            return stage.J @ sum(Pk @ A @ Pk for Pk in P.values())
        dims = [int(round(np.trace(Pk))) for Pk in P.values()]
        out["block_diagonal"] = _closure_dim([block(), block()], stage.E)
        out["block_diagonal_expected"] = sum(d * (d + 1) // 2 for d in dims)
    out["expected"] = {"random_pair": SP26, "hecke_commuting": len(stage.projectors) if stage.projectors else None}
    return out
