"""R2C-01 (Astra, 2026-09-18; verified independently 2026-09-19): 6D chirality assignments for M1's ten oriented bifundamentals.

Convention: n_L - n_R (4D) = eps_ij d_ij with d_ij = m_i - m_j.  Selection rule (Spin(6, C) = SL_4(C)): the ordinary scalar
bilinear needs OPPOSITE 6D chiralities (4 (x) 4^* contains 1, 4 (x) 4 does not); an odd Clifford (vector) insertion needs EQUAL
chiralities.  Result: preserving the six family indices forces the six family chiralities to +1, while the four elementary-
scalar Yukawa contractions need eps_cL = -eps_ca = -eps_cb and eps_Ld = -eps_ad = -eps_bd — no assignment does both
(1024 / 16 family-preserving / 64 scalar-compatible / 0 both).  With all family chiralities equal, the VECTOR insertion is
allowed: a Higgs that is an internal gauge component (the flux-recombination picture) is compatible with the three families;
an elementary 6D scalar Higgs is not.  Exact, conditional on the declared class (one Weyl field per sector, original flux)."""
import itertools, functools, numpy as np
from fractions import Fraction

SECTORS = ("cL", "ca", "cb", "cd", "La", "Lb", "Ld", "ab", "ad", "bd")
DIMS = {"cL": 6, "ca": 3, "cb": 3, "cd": 3, "La": 2, "Lb": 2, "Ld": 2, "ab": 1, "ad": 1, "bd": 1}
M1_DEGREES = {"c": 0, "L": -3, "a": 3, "b": 3, "d": 0}
FAMILY_TARGET = {"cL": 3, "ca": -3, "cb": -3, "Ld": -3, "ad": 3, "bd": 3}          # Q, u^c, d^c, L, nu^c, e^c
SCALAR_TRIANGLES = (("cL", "ca"), ("cL", "cb"), ("Ld", "ad"), ("Ld", "bd"))         # Q u^c H_u, Q d^c H_d, L nu^c H_u, L e^c H_d

def sector_degrees(m=M1_DEGREES): return {s: m[s[0]] - m[s[1]] for s in SECTORS}

def clifford_6d():
    s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex); I = np.eye(2)
    kron = lambda *ms: functools.reduce(np.kron, ms)
    g = [kron(s1, I, I), kron(s2, I, I), kron(s3, s1, I), kron(s3, s2, I), kron(s3, s3, s1), kron(s3, s3, s2)]
    Gs = 1j * g[0] @ g[1] @ g[2] @ g[3] @ g[4] @ g[5]; C = g[1] @ g[3] @ g[5]
    return {"gamma": g, "Gamma_star": Gs, "C": C, "clifford_ok": all(np.allclose(g[a] @ g[b] + g[b] @ g[a], 2 * np.eye(8) * (a == b)) for a in range(6) for b in range(6)),
            "C_ok": all(np.allclose(C @ g[a] @ np.linalg.inv(C), -g[a].T) for a in range(6))}

def bilinear_selection(eps1, eps2, kind="scalar"):
    """Rank of the algebraic contraction between Weyl projectors: scalar P_{-e1} P_{e2}, transpose P_{e1}^T C P_{e2}, vector P_{-e1} gamma P_{e2}."""
    cl = clifford_6d(); Gs, C, g = cl["Gamma_star"], cl["C"], cl["gamma"]; P = lambda e: (np.eye(8) + e * Gs) / 2
    M = {"scalar": P(-eps1) @ P(eps2), "transpose": P(eps1).T @ C @ P(eps2), "vector": P(-eps1) @ g[4] @ P(eps2)}[kind]
    r = int(np.linalg.matrix_rank(M, tol=1e-12)); return {"kind": kind, "rank": r, "allowed": r > 0, "rule": "opposite chiralities" if kind != "vector" else "equal chiralities"}

def signed_index(eps, degree): return eps * degree

def enumerate_chiralities(m=M1_DEGREES):
    d = sector_degrees(m); counts = dict(total=0, family=0, full_ledger=0, scalar=0, both=0, family_and_LaLb_equal=0); witness_both = None; family_rows = []
    for eps in itertools.product((1, -1), repeat=len(SECTORS)):
        E = dict(zip(SECTORS, eps)); idx = {s: E[s] * d[s] for s in SECTORS}; counts["total"] += 1
        fam = all(idx[s] == FAMILY_TARGET[s] for s in FAMILY_TARGET); full = all(idx[s] == d[s] for s in SECTORS)
        scalar = all(E[a] == -E[b] for a, b in SCALAR_TRIANGLES)
        counts["family"] += fam; counts["full_ledger"] += full; counts["scalar"] += scalar; counts["both"] += fam and scalar
        if fam: counts["family_and_LaLb_equal"] += (E["La"] == E["Lb"]); family_rows.append(E)
        if fam and scalar and witness_both is None: witness_both = E
    conflict = {"family_requires": "eps_cL = eps_ca = eps_cb = +1 (and eps_Ld = eps_ad = eps_bd = +1)", "scalar_requires": "eps_cL = -eps_ca = -eps_cb (and eps_Ld = -eps_ad = -eps_bd)", "witness_pair": ("cL", "ca")}
    return {"counts": counts, "family_rows": family_rows, "no_go": counts["both"] == 0, "conflict": conflict, "vector_higgs_allowed_with_family_signs": bilinear_selection(1, 1, "vector")["allowed"]}

def gravitational_p2_coefficient(E):
    """[p2] I8 = -n_grav/1440 with n_grav = sum eps_ij dim(R_ij) the signed 6D representation dimension (irreducible term)."""
    n = sum(E[s] * DIMS[s] for s in SECTORS); return {"n_grav": n, "p2_coefficient": Fraction(-n, 1440)}

def m1_scalar_higgs_no_go():
    r = enumerate_chiralities(); table = [(E["cd"], E["ab"], gravitational_p2_coefficient(E)) for E in r["family_rows"] if E["La"] == E["Lb"] == 1]
    return {"counts": r["counts"], "no_go": r["no_go"], "conflict": r["conflict"], "ledger_preserving_p2_table": table, "internal_vector_higgs": "allowed by the selection rule; needs an action, its physical mode operator, anomaly cancellation and a stable background (OPEN)",
            "status": "EXACT, conditional on the declared class (R2C-01)"}
