"""INT-04b (HOPF-02 §§2, 7): coprime block projectors of T2 on H_1(X0(143)) and the quaternionic-multiplicity gate.

The characteristic polynomial of T2 on H_1 (26-dim) is x^2 (x+2)^4 g(x)^2 s(x)^2 with g quartic and s sextic.  Since s
and x (x+2) g are coprime, Bezout gives complementary rational projectors P12 (sextic block, rank 12) and P14 (the
elliptic + old + quartic complement, rank 14); T2 being self-adjoint for the Hodge/Petersson form, the projections are
orthogonal and the unit sphere of H^1(R) is the join S^11 * S^13 = S^25.  Matrices act on the package homology basis
(transpose for the dual cohomology basis).  Quaternionic gate: a quaternionic structure (I, J anticommuting complex
structures) commuting with T must act on each real T-eigenspace; a 2-dimensional real eigenspace admits I but no
anticommuting J with J^2 = -1, so the quartic block (four real eigenvalues, each of multiplicity 2) has NO commuting
quaternionic structure; doubling the multiplicity removes the obstruction (a control, not the original modes)."""
import sympy as sp

def _evaluate(poly, x, M):
    out = sp.zeros(M.rows)
    for c in sp.Poly(poly, x).all_coeffs(): out = out * M + c * sp.eye(M.rows)
    return out

def coprime_block_projectors():
    from .. import hecke as H
    x = sp.Symbol("x"); T = sp.Matrix(H.cuspidal_hecke(2))
    g = sum(sp.Integer(c) * x ** j for j, c in enumerate(H.G4)); s = sum(sp.Integer(c) * x ** j for j, c in enumerate(H.H6)); rest = x * (x + 2) * g
    ba, bb, gg = sp.gcdex(s, rest, x)
    if gg != 1: raise ArithmeticError("sextic and complement not coprime")
    P12 = _evaluate(bb * rest, x, T); P14 = sp.eye(26) - P12
    return {"T2": T, "P12": P12, "P14": P14, "sextic": s, "quartic": g, "rank12": P12.rank(), "rank14": P14.rank(),
            "idempotent": P12 * P12 == P12, "complementary": P12 * P14 == sp.zeros(26), "sphere_join": "S^11 * S^13 = S^25"}

def quaternionic_multiplicity_gate(real_eigenspace_dimensions):
    """A T-commuting quaternionic structure needs every real T-eigenspace dimension divisible by 4; returns the obstruction
    witness (the first eigenspace failing) or None.  Quartic block: (2, 2, 2, 2) -> obstructed; doubled: (4, 4, 4, 4) -> not."""
    bad = [d for d in real_eigenspace_dimensions if d % 4]
    return {"obstructed": bool(bad), "witness_dimensions": bad, "rule": "each commuting real eigenspace must have dimension = 0 mod 4"}
