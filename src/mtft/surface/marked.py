"""mtft.surface.marked — marked two-factor geometry: decomposition, degeneracy maps, readouts (EXACT).

Astra's PLANCK-02 witness: an additive synthetic Hamiltonian h_syn = 60 Z⊗I + 48 I⊗Z has the same
spectrum {−108,−12,12,108} as the connected arithmetic interaction, so no spectrum-only action can
tell them apart.  Marked observables can (MARKED-03).  This module keeps the factor embeddings:

  decomposition    h = τ(h) I + h_A + h_B + C(h)   (partial traces in the G-orthonormal frame)
  η(h)             ‖C(h)‖² / ‖h − τ(h) I‖²         (1 for the arithmetic interaction, 0 for h_syn)
  degeneracy maps  J_{q,i} = I⊗e_i (adjoin q), J_{p,j} = e_j⊗I; gates J♯J = (r+1)I, J₀♯J₁ = a_r I
  transfers        K_ij = J_i♯ H J_j; diagonal transfers are blind to α, K₀₁ − K₁₀ = 320 α R reads it
  two objectives   least leakage Tr V_i(H_α) = 409600/81 + 7200 α²  (minimum α = 0)
                   dist²(H_α, Herm ℝ[U_p,U_q]) = 57600/41 (α − 1)²  (zero α = 1)
so the choice of principle is explicit.  Values quoted are for 11a1 with (p, q) = (13, 17).
"""
from __future__ import annotations

from typing import Dict

import sympy as sp

from .oldsector import connected_part


def _orthonormal(Gp: sp.Matrix, Gq: sp.Matrix):
    Lf = sp.kronecker_product(Gp.cholesky(), Gq.cholesky())
    return Lf, Lf.inv()


def decomposition(h: sp.Matrix, Gp: sp.Matrix, Gq: sp.Matrix) -> Dict:
    if any(isinstance(x, (float, complex)) or getattr(x, "is_Float", False) or (x.is_number and not x.is_real) for x in h):
        raise TypeError("decomposition: exact real (rational/symbolic) entries required; complex/float inputs are unsupported")
    Lf, Li = _orthonormal(Gp, Gq)
    Ho = Lf.T * h * Li.T
    I2 = sp.eye(2)
    tr1 = sp.Matrix(2, 2, lambda i, j: sum(Ho[2 * k + i, 2 * k + j] for k in range(2)))
    tr2 = sp.Matrix(2, 2, lambda i, j: sum(Ho[2 * i + k, 2 * j + k] for k in range(2)))
    tau = Ho.trace() / 4
    hA = sp.kronecker_product(tr2 / 2 - tau * I2, I2)
    hB = sp.kronecker_product(I2, tr1 / 2 - tau * I2)
    C = Ho - tau * sp.eye(4) - hA - hB
    norm2 = lambda M: (M.T * M).trace() / 4          # normalised trace tau = Tr/4 (Astra's convention)
    cen = norm2(Ho - tau * sp.eye(4))
    return {"tau": tau, "local_p_norm2": norm2(hA), "local_q_norm2": norm2(hB), "connected_norm2": norm2(C),
            "centered_norm2": cen, "eta": (norm2(C) / cen) if cen != 0 else None,
            "orthogonal_sum_identity": sp.simplify(cen - norm2(hA) - norm2(hB) - norm2(C)) == 0,
            "connected_orthonormal_frame": C}


def degeneracy_maps(sector: Dict) -> Dict:
    """Maps from the two-dimensional p-sector (adjoin q) and q-sector (adjoin p) into the 4-dim sector,
    with metric adjoints J♯ = G_src⁻¹ Jᵀ G_tgt and the two map gates."""
    I2 = sp.eye(2)
    e = [sp.Matrix([[1], [0]]), sp.Matrix([[0], [1]])]
    G = sector["G"]
    Gp, Gq = sector["local_p"]["G"], sector["local_q"]["G"]
    Jq = [sp.kronecker_product(I2, e[i]) for i in range(2)]      # 143-sector -> 2431 along q
    Jp = [sp.kronecker_product(e[j], I2) for j in range(2)]      # 187-sector -> 2431 along p
    adj = lambda J, Gs: Gs.inv() * J.T * G
    q, p = sector["q"], sector["p"]
    aq, ap = sector["local_q"]["B"][0, 0], sector["local_p"]["B"][0, 0]
    gates = {
        "Jq_degree": all((adj(Jq[i], Gp) * Jq[i] - (q + 1) * I2).is_zero_matrix for i in range(2)),
        "Jq_cross_is_a_q": (adj(Jq[0], Gp) * Jq[1] - aq * I2).is_zero_matrix,
        "Jp_degree": all((adj(Jp[j], Gq) * Jp[j] - (p + 1) * I2).is_zero_matrix for j in range(2)),
        "Jp_cross_is_a_p": (adj(Jp[0], Gq) * Jp[1] - ap * I2).is_zero_matrix,
    }
    return {"J_q": Jq, "J_q_adj": [adj(J, Gp) for J in Jq], "J_p": Jp, "J_p_adj": [adj(J, Gq) for J in Jp], "gates": gates}


def readout_coefficient(sector: Dict):
    """c with K₀₁ − K₁₀ = c α R along the q-maps: c = (q+1)² − a_q²  (320 at q = 17, a_q = −2).
    CC-24: v0.28.0 hardcoded 320, so the reversed (17,13) sector read 9α/16 (found by Astra)."""
    q = sector["q"]
    aq = sector["local_q"]["B"][0, 0]
    return sp.Integer((q + 1) ** 2 - aq ** 2)


def _require_exact(M: sp.Matrix, name: str = "input") -> None:
    if any(not isinstance(x, (sp.Rational, sp.Integer)) and not (x.is_number and x.is_rational) for x in M):
        raise TypeError(f"{name}: exact rational entries required (no floats/complex); symbols allowed only for family parameters")


def cross_map_readout(sector: Dict, H: sp.Matrix, maps: Dict) -> Dict:
    """K_ij = J_i♯ H J_j along the q-maps; α = Tr(R♯(K₀₁−K₁₀)) / (c Tr(R♯R)) with c = readout_coefficient."""
    K = [[maps["J_q_adj"][i] * H * maps["J_q"][j] for j in range(2)] for i in range(2)]
    Gp = sector["local_p"]["G"]
    R = sector["local_p"]["R"]
    Rs = Gp.inv() * R.T * Gp
    return {"K": K, "diagonal_equal": (K[0][0] - K[1][1]).is_zero_matrix,
            "antisym_transfer": K[0][1] - K[1][0], "R_pairing": (Rs * (K[0][1] - K[1][0])).trace(), "R_norm2": (Rs * R).trace()}


def least_leakage(sector: Dict, H: sp.Matrix, maps: Dict, branch: int = 0):
    """Tr V_i(H) with E_i(H) = J_i♯ H J_i / (q+1) and V_i = E_i(H²) − E_i(H)² (positive metric variance)."""
    d = sector["q"] + 1
    J, Js = maps["J_q"][branch], maps["J_q_adj"][branch]
    E = lambda X: Js * X * J / d
    return (E(H * H) - E(H) * E(H)).trace()


def hermitian_hecke_distance2(sector: Dict, H: sp.Matrix):
    """Squared metric Hilbert–Schmidt distance from H to the real span of Herm(Q[U_p,U_q])."""
    herm = sector["herm"]
    G = sector["G"]
    Lf, Li = _orthonormal(sector["local_p"]["G"], sector["local_q"]["G"])
    on = lambda M: Lf.T * M * Li.T
    span = [on(herm(M)).reshape(16, 1) for M in (sp.eye(4), sector["U_p"], sector["U_q"], sector["U_p"] * sector["U_q"])]
    A = sp.Matrix.hstack(*span)
    v = on(H).reshape(16, 1)
    c = (A.T * A).inv() * A.T * v
    r = v - A * c
    return (r.T * r)[0, 0] / 4                        # tau-norm


def two_objectives(sector: Dict, alpha=None) -> Dict:
    a = sp.Symbol("alpha") if alpha is None else alpha
    H = sector["H_prod"] + a * sector["skew_correction"]
    maps = degeneracy_maps(sector)
    return {"least_leakage_trace": sp.expand(least_leakage(sector, H, maps)),
            "herm_hecke_distance2": sp.factor(hermitian_hecke_distance2(sector, H)),
            "readout_alpha": sp.simplify(cross_map_readout(sector, H, maps)["R_pairing"] / (readout_coefficient(sector) * cross_map_readout(sector, H, maps)["R_norm2"]))}
