"""INT-12 (AXG-04 §5): the C3X ordinary spin-bordism certificate, Omega_7^Spin(B(SU3 x SU2 x U1^2)) = 0.

Data: the Sq^2 matrices over F2 between the monomial bases of H^4, H^6, H^8 of BG (bases below), recorded by AXG-04.
Recomputed here: ranks 2 and 7, and kernel(Sq^2: H^6 -> H^8) = image(Sq^2: H^4 -> H^6) (dimension 2), so the only
total-degree-seven E2 term of the Atiyah–Hirzebruch spectral sequence, E2^{6,1} = H^6(BG; Z/2) (dim 9), dies at E3.
Scope: ordinary Spin-BG category; it removes a torsion-anomaly obstruction, it is NOT the global GS construction (H-23)."""
import sympy as sp

BASES = {"H4": ["h^2", "h*u", "u^2", "C", "L"], "H6": ["h^3", "h^2*u", "h*u^2", "u^3", "h*C", "u*C", "h*L", "u*L", "D"],
         "H8": ["h^4", "h^3*u", "h^2*u^2", "h^2*C", "h^2*L", "h*u^3", "h*u*C", "h*u*L", "h*D", "u^4", "u^2*C", "u^2*L", "u*D", "C^2", "C*L", "L^2"]}
SQ2_H4_H6 = [[0,0,0,0,0],[0,1,0,0,0],[0,1,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,1,0]]
SQ2_H6_H8 = [[1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,1,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0]]

def _rank_f2(M): return sp.Matrix(M).rank(iszerofunc=lambda x: x % 2 == 0) if False else _rank_mod2(M)

def _rank_mod2(M):
    A = [[int(x) % 2 for x in row] for row in M]; r = 0; rows, cols = len(A), len(A[0])
    for c in range(cols):
        p = next((i for i in range(r, rows) if A[i][c]), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        for i in range(rows):
            if i != r and A[i][c]: A[i] = [(a + b) % 2 for a, b in zip(A[i], A[r])]
        r += 1
    return r

def c3x_spin_bordism_certificate():
    d1, d2 = SQ2_H4_H6, SQ2_H6_H8; r1, r2 = _rank_mod2(d1), _rank_mod2(d2); dimH6 = len(BASES["H6"])
    composite_zero = all(sum(d2[i][k] * d1[k][j] for k in range(dimH6)) % 2 == 0 for i in range(len(d2)) for j in range(len(d1[0])))
    kernel_dim = dimH6 - r2; image_dim = r1
    return {"dims": (len(BASES["H4"]), dimH6, len(BASES["H8"])), "rank_H4_H6": r1, "rank_H6_H8": r2, "composite_zero": composite_zero,
            "kernel_dim": kernel_dim, "image_dim": image_dim, "kernel_equals_image": composite_zero and kernel_dim == image_dim,
            "E3_6_1": dimH6 - r1 - r2 if composite_zero else None, "omega7_spin_BG": 0 if (composite_zero and kernel_dim == image_dim) else None,
            "scope": "ordinary Spin-BG bordism; necessary for a GS completion, not the construction"}
