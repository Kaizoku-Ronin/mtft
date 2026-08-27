"""Quadratic-Hamiltonian layer on the X_0(143) Hodge stage.

The graph couplings V are graph-self-adjoint but not Hodge-self-adjoint;
only the Hodge-Hermitian part A = (V + V^{dagger_H})/2 enters the scalar
quadratic energy H_V(x) = x^T G A x / 2.  This module implements the
Hodge adjoint, the Hamiltonian J-split A = A_+ + A_-, the pairing-strength
parameter rho, the 13 Williamson/symplectic frequencies of the flow J A,
and the oldspace / Hecke-block routing of the pairing channel A_-.

Classes: operator identities CERTIFIED(tol) on the certified period stage;
the degree channel is the EXACT null control (V_degree = 3I, A_- = 0);
physical (Hz/mass/SM) interpretation remains PHENO and is not made here.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from .physics import (hodge_structure_hecke, hodge_metric_hecke,
                      graph_coupling, metric_hs_norm)
from .involutions import oldspace_projector, hecke_blocks

POTENTIALS = ("degree", "width", "distance")


@lru_cache(maxsize=4)
def _stage(dps: int = 60):
    G = hodge_metric_hecke(dps)
    J = hodge_structure_hecke(dps)
    w, Q = np.linalg.eigh((G + G.T) / 2)
    Gh = Q @ np.diag(np.sqrt(w)) @ Q.T
    Ghi = Q @ np.diag(1 / np.sqrt(w)) @ Q.T
    return G, J, Gh, Ghi


def hodge_adjoint(V, dps: int = 60):
    """V^{dagger_H} = G^{-1} V^T G."""
    G, *_ = _stage(dps)
    return np.linalg.solve(G, np.asarray(V, float).T @ G)


def hermitian_split(V, dps: int = 60):
    """(A, K): Hodge-Hermitian and Hodge-skew parts of V."""
    adj = hodge_adjoint(V, dps)
    V = np.asarray(V, float)
    return 0.5 * (V + adj), 0.5 * (V - adj)


def hamiltonian_split(A, dps: int = 60):
    """(A_+, A_-): J-linear and J-antilinear parts of the Hamiltonian A."""
    _, J, *_ = _stage(dps)
    A = np.asarray(A, float)
    return 0.5 * (A - J @ A @ J), 0.5 * (A + J @ A @ J)


def channel_report(potential: str = "width", dps: int = 60):
    """Full Hamiltonianized channel report for one graph potential."""
    G, J, Gh, Ghi = _stage(dps)
    V, _, _ = graph_coupling(potential)
    A, K = hermitian_split(V, dps)
    Ap, Am = hamiltonian_split(A, dps)
    hs = lambda X: metric_hs_norm(X, G)
    nV, nA, nAm = hs(V), hs(A), hs(Am)
    # inertia of the pairing form via the G-orthonormal frame
    S = Gh @ Am @ Ghi
    ev = np.linalg.eigvalsh((S + S.T) / 2)
    inertia = (int((ev > 1e-10 * max(1, abs(ev).max())).sum()),
               int((ev < -1e-10 * max(1, abs(ev).max())).sum()))
    return {
        "potential": potential,
        "hermiticity_defect_rel": hs(K) / nV if nV else 0.0,
        "hamiltonian_antilinear_fraction": nAm / nA if nA else 0.0,
        "hamiltonian_power_fraction": (nAm / nA) ** 2 if nA else 0.0,
        "pairing_inertia": inertia,
        "epistemic": "CERTIFIED(tol) diagnostic; degree channel EXACT zero",
    }


def pairing_stability(potential: str = "width", dps: int = 60) -> float:
    """rho = ||A_+^{-1/2} A_- A_+^{-1/2}||_2 in the G-orthonormal frame.

    rho < 1 certifies the full quadratic Hamiltonian stays positive.
    Anchors: width 0.1234286299, distance 0.4248813827, degree 0.
    """
    _, _, Gh, Ghi = _stage(dps)
    V, _, _ = graph_coupling(potential)
    A, _ = hermitian_split(V, dps)
    Ap, Am = hamiltonian_split(A, dps)
    Sp = Gh @ Ap @ Ghi
    Sm = Gh @ Am @ Ghi
    Sp = (Sp + Sp.T) / 2
    w, Q = np.linalg.eigh(Sp)
    if w.min() <= 0:
        raise ArithmeticError("A_+ not positive definite")
    Ri = Q @ np.diag(1 / np.sqrt(w)) @ Q.T
    C = Ri @ ((Sm + Sm.T) / 2) @ Ri
    return float(np.linalg.norm(C, 2))


def symplectic_frequencies(potential: str = "width", dps: int = 60):
    """The 13 positive Williamson frequencies of the flow F = J A.

    spec(F) = {±i w_k}; degree gives w = 3 (x13) exactly.
    """
    _, J, *_ = _stage(dps)
    V, _, _ = graph_coupling(potential)
    A, _ = hermitian_split(V, dps)
    ev = np.linalg.eigvals(J @ A)
    assert abs(ev.real).max() < 1e-8 * max(1.0, abs(ev.imag).max())
    return np.sort(np.abs(ev.imag))[1::2]


def _projector(Bcols, G):
    B = np.array(Bcols, float)
    return B @ np.linalg.solve(B.T @ G @ B, B.T @ G)


def oldspace_routing(potential: str = "width", dps: int = 60):
    """Power split of A_- into new/old sectors via the exact P_old.

    Honest negative on record: the ghost carries only ~12.3% (width) and
    ~14.8% (distance) of the pairing power; new->new dominates.
    """
    G, *_ = _stage(dps)
    V, _, _ = graph_coupling(potential)
    A, _ = hermitian_split(V, dps)
    _, Am = hamiltonian_split(A, dps)
    Po = np.array([[float(x) for x in r] for r in oldspace_projector()])
    Pn = np.eye(26) - Po
    hs2 = lambda X: metric_hs_norm(X, G) ** 2
    tot = hs2(Am)
    parts = {"new_new": hs2(Pn @ Am @ Pn), "old_old": hs2(Po @ Am @ Po),
             "new_old": hs2(Pn @ Am @ Po), "old_new": hs2(Po @ Am @ Pn)}
    out = {k: v / tot for k, v in parts.items()}
    out["old_involved"] = 1.0 - out["new_new"]
    out["closure_residual"] = abs(sum(parts.values()) / tot - 1.0)
    return out


def hecke_block_routing(potential: str = "width", dps: int = 60):
    """4x4 Hecke-block power matrix of A_- over (ell, ghost, q4, q6).

    Width is mostly intra-block (~75.4%); distance is majority inter-block
    (~54.5%) — the two nontrivial potentials route differently through the
    arithmetic sectors.
    """
    G, *_ = _stage(dps)
    V, _, _ = graph_coupling(potential)
    A, _ = hermitian_split(V, dps)
    _, Am = hamiltonian_split(A, dps)
    blocks = hecke_blocks()
    names = ("ell", "ghost", "q4", "q6")
    P = {b: _projector([[float(v[i]) for v in blocks[b]] for i in range(26)], G)
         for b in names}
    hs2 = lambda X: metric_hs_norm(X, G) ** 2
    tot = hs2(Am)
    M = {(a, b): hs2(P[a] @ Am @ P[b]) / tot for a in names for b in names}
    diag = sum(M[(a, a)] for a in names)
    return {"matrix": M, "intra_block": diag, "inter_block": 1.0 - diag,
            "closure_residual": abs(sum(M.values()) - 1.0)}


__all__ = ["POTENTIALS", "hodge_adjoint", "hermitian_split",
           "hamiltonian_split", "channel_report", "pairing_stability",
           "symplectic_frequencies", "oldspace_routing",
           "hecke_block_routing"]
