"""mtft.surface.oldsector — exact two-prime oldform sectors and the Hermitian-part selection rule.

For a weight-2 newform f of level N0 and primes p, q not dividing N0, the level-N0 pq oldform
sector spanned by f(z), p f(pz), q f(qz), pq f(pqz) is a tensor product of local 2-dimensional
sectors.  On the local basis (f(z), p f(pz)) the operator U_p and the Petersson Gram matrix are

    B_p = [[a_p, p], [-1, 0]],      G_p = [[p+1, a_p], [a_p, p+1]]

(Schulze-Pillot–Yenirce factorisation), with metric adjoint B_p† = G_p^-1 B_pᵀ G_p, Hermitian
part A_p = (B_p + B_p†)/2 and skew part R_p = (B_p − B_p†)/2.  U_p = B_p ⊗ I, U_q = I ⊗ B_q,
G = G_p ⊗ G_q; all cross-prime commutators vanish.

Selection rule tested here (EXACT, within the four-word family H_α = A_p⊗A_q + α R_p⊗R_q):
Herm(Q[U_p, U_q]) is 4-dimensional and contains exactly ONE interaction direction, α = 1, i.e.
H_pq = Herm(U_p U_q) = A_p⊗A_q + R_p⊗R_q.  The product of Hermitian parts (α = 0) is not the
Hermitian part of any Hecke element.  The interaction term R_p⊗R_q is a product of the skew
(non-normal) parts, which vanish for good primes: in this class, coupling lives only at the
level's primes.  Status: a stated principle (Hamiltonians are Hermitian parts of Hecke algebra
elements), canonical and basis-independent, not derived from the curve.  Verified against
Astra's composite-correspondence study (2026-09-10): connected spectrum {−108, −12, 12, 108}
and raw charpoly (x²−92x−1584)(x²+100x−1584) for f = 11a1, p, q = 13, 17.
"""
from __future__ import annotations

from typing import Dict, Tuple

import sympy as sp


def _is_prime(n: int) -> bool:
    return n > 1 and all(n % k for k in range(2, int(n ** 0.5) + 1))


def local_block(p: int, a_p: int) -> Dict[str, sp.Matrix]:
    if not _is_prime(p):
        raise ValueError(f"p={p} is not prime")
    if a_p * a_p > 4 * p:
        raise ValueError(f"a_p={a_p} violates the Hasse bound |a_p| <= 2 sqrt(p)")
    B = sp.Matrix([[a_p, p], [-1, 0]])
    G = sp.Matrix([[p + 1, a_p], [a_p, p + 1]])
    Bd = G.inv() * B.T * G
    return {"B": B, "G": G, "B_adj": Bd, "A": (B + Bd) / 2, "R": (B - Bd) / 2}


def two_prime_sector(p: int, q: int, a_p: int, a_q: int) -> Dict:
    if p == q:
        raise ValueError("the two added primes must be distinct")
    Lp, Lq = local_block(p, a_p), local_block(q, a_q)
    I2 = sp.eye(2)
    K = sp.kronecker_product
    G = K(Lp["G"], Lq["G"])
    Up, Uq = K(Lp["B"], I2), K(I2, Lq["B"])
    herm = lambda M: (M + G.inv() * M.T * G) / 2
    H_add = K(Lp["A"], I2) + K(I2, Lq["A"])
    H_prod = K(Lp["A"], Lq["A"])
    H_pq = herm(Up * Uq)
    return {"p": p, "q": q, "local_p": Lp, "local_q": Lq, "G": G, "U_p": Up, "U_q": Uq, "herm": herm,
            "H_add": H_add, "H_prod": H_prod, "H_pq": H_pq,
            "skew_correction": K(Lp["R"], Lq["R"]),      # Herm(U_p U_q) - Herm(U_p) Herm(U_q)
            "interaction_term": K(Lp["R"], Lq["R"]),     # compatibility alias for skew_correction
            "cross_commutators_zero": all((X * Y - Y * X).is_zero_matrix for X in (Up, herm(Up) * 2 - Up) for Y in (Uq, herm(Uq) * 2 - Uq))}


def connected_part(H: sp.Matrix, Gp: sp.Matrix, Gq: sp.Matrix) -> sp.Matrix:
    """Remove scalar and single-factor terms by partial traces in the G-orthonormal frame
    (returned in orthonormal coordinates; the spectrum is frame-independent)."""
    Lp, Lq = Gp.cholesky(), Gq.cholesky()
    Lf = sp.kronecker_product(Lp, Lq)
    Ho = Lf.T * H * Lf.inv().T
    I2 = sp.eye(2)
    tr1 = sp.Matrix(2, 2, lambda i, j: sum(Ho[2 * k + i, 2 * k + j] for k in range(2)))
    tr2 = sp.Matrix(2, 2, lambda i, j: sum(Ho[2 * i + k, 2 * j + k] for k in range(2)))
    return Ho - sp.kronecker_product(I2, tr1) / 2 - sp.kronecker_product(tr2, I2) / 2 + Ho.trace() * sp.eye(4) / 4


def schmidt_rank(C: sp.Matrix) -> int:
    R = sp.Matrix(4, 4, lambda I_, J_: C[2 * (I_ // 2) + (J_ // 2), 2 * (I_ % 2) + (J_ % 2)])
    return int(R.rank())


def hermitian_hecke_selection(sector: Dict) -> Dict:
    """Which α in H_α = A_p⊗A_q + α R_p⊗R_q lies in Herm(Q[U_p, U_q])?  Exactly α = 1."""
    herm = sector["herm"]
    Up, Uq = sector["U_p"], sector["U_q"]
    span = [herm(M) for M in (sp.eye(4), Up, Uq, Up * Uq)]
    cols = sp.Matrix.hstack(*[m.reshape(16, 1) for m in span])
    a = sp.Symbol("alpha")
    Ha = sector["H_prod"] + a * sector["interaction_term"]
    cs = sp.symbols("c0:4")
    eqs = list(cols * sp.Matrix(cs) - Ha.reshape(16, 1))
    solution = sp.solve(eqs, list(cs) + [a], dict=True)
    if not solution:
        status, alpha = "no_member_in_family", None
    elif a in solution[0] and solution[0][a].free_symbols == set():
        status, alpha = "unique", solution[0][a]
    elif a in solution[0]:
        status, alpha = "not_unique", None
    else:                                             # alpha absent: family member for every alpha
        status, alpha = "not_identifiable", None
    return {"dim_Herm_Hecke": int(cols.rank()), "status": status, "alpha_selected": alpha,
            "H_selected": sector["H_pq"] if status == "unique" else None}


def good_prime_replacement_control(sector: Dict, a_ell: int) -> Dict:
    """Replace the second local operator by the scalar a_ell I (a good prime T_ell acts as the
    scalar a_ell on the fixed newform oldclass): the skew correction and the whole connected
    term vanish.  This is the representation-specific 'bad-prime only' statement."""
    I2 = sp.eye(2)
    Lp = sector["local_p"]
    Up = sp.kronecker_product(Lp["B"], I2)
    Uq = a_ell * sp.eye(4)
    G = sector["G"]
    herm = lambda M: (M + G.inv() * M.T * G) / 2
    Hpq = herm(Up * Uq)
    skew = Hpq - herm(Up) * herm(Uq)
    C = connected_part(Hpq, Lp["G"], sector["local_q"]["G"])
    return {"skew_correction_zero": skew.is_zero_matrix, "connected_term_zero": C.is_zero_matrix}
