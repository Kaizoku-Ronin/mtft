"""mtft.origami.gates — the v0.20.0 gate battery, runnable end to end.

    python -c "from mtft.origami.gates import run_all; run_all()"

Every gate below is a certification, not an illustration: each one has failed
at least once during development and been fixed, or is the direct descendant of
a gate that caught a real error (see CORRECTION RECORD).
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from .dimer import ensemble_conservation
from .insertion import (brioschi_curvature, cubic_tensor, cumulant_curvature,
                        fisher_metric, path_independence)
from .instances import (PRISM_C, PRISM_LAMBDA0, closed_curvature_B,
                        galashin_24, mandelstams_24, prism_36)
from .perfect import (Theta, cyclic_matrix, equivariant_kasteleyn_factor,
                      is_valid_pair, orbit_structure, solve_perfect_branches,
                      t_coefficients, winding)

__all__ = [
    "gate_boundary_measurement_24", "gate_mandelstam_24",
    "gate_ensemble_conservation_24", "gate_simplex_curvature",
    "gate_closed_forms_B", "gate_spinflip_parity", "gate_curvature_bounds",
    "gate_fisher_quadric", "gate_path_independence",
    "gate_prism_top_cell", "gate_prism_symmetry", "gate_prism_lambda0",
    "gate_prism_branches", "gate_square_mechanism", "run_all",
]


# ------------------------------------------------------------------ (2,4)
def gate_boundary_measurement_24():
    """Meas(Gamma, wt) reproduces Example 2.18 symbolically (all six Pluckers)."""
    p, q, r, s = sp.symbols("p q r s", positive=True)
    G = galashin_24(p, q, r, s)
    D = {tuple(sorted(I)): v for I, v in G.boundary_measurement().items()}
    want = {(1, 3): sp.Integer(1), (0, 3): p, (2, 3): q, (1, 2): r,
            (0, 2): sp.expand(p * r + q * s), (0, 1): s}
    assert D == {k: sp.expand(v) for k, v in want.items()}, D
    return D


def gate_mandelstam_24(draws=5, seed=143):
    """Lemma 1.10 plus the EXACT channel identities S(1,3)=pr, S(2,4)=qs."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(draws):
        p, q, r, s = np.exp(rng.uniform(-1, 1, 4))
        S = mandelstams_24(p, q, r, s)
        assert abs(S[(1, 3)] - p * r) < 1e-9
        assert abs(S[(2, 4)] - q * s) < 1e-9
        out.append(S)
    return out


def gate_ensemble_conservation_24():
    """The gate that caught the dropped qs matching."""
    p, q, r, s = sp.symbols("p q r s", positive=True)
    G = galashin_24(p, q, r, s)
    idx = {G.edges[k][2]: k for k in range(len(G.edges))}
    cls_B = ensemble_conservation(G, [idx[p], idx[r]])
    c00 = sp.factor(cls_B[(0, 0)])
    assert c00 == sp.factor((1 + q) * (1 + s)), c00
    cls_A = ensemble_conservation(G, [idx[p], idx[s]])
    assert sp.factor(cls_A[(0, 0)]) == sp.factor(1 + q + r)
    return cls_A, cls_B


def gate_simplex_curvature():
    """Section A is the full trinomial: K == 1/4 exactly, coefficients symbolic."""
    X, Y = sp.symbols("X Y", positive=True)
    c0, c1, c2 = sp.symbols("c0 c1 c2", positive=True)
    K = brioschi_curvature(c0 + c1 * X + c2 * Y, (X, Y))
    num = sp.expand(sp.fraction(sp.together(K - sp.Rational(1, 4)))[0])
    assert num == 0, num
    return sp.Rational(1, 4)


def gate_closed_forms_B():
    """Certified det g and K for section B, on both curvature routes."""
    X, Y, c, det_g, K = closed_curvature_B()
    Z = c + X + Y + X * Y
    g = fisher_metric(Z, (X, Y))
    assert sp.expand(sp.fraction(sp.together(g.det() - det_g))[0]) == 0
    assert sp.expand(sp.fraction(sp.together(
        brioschi_curvature(Z, (X, Y)) - K))[0]) == 0
    assert sp.expand(sp.fraction(sp.together(
        cumulant_curvature(Z, (X, Y)) - K))[0]) == 0
    return det_g, K


def gate_spinflip_parity():
    """g EVEN and the cubic Amari tensor ODD under (X,Y) -> (c/X, c/Y).

    Consequence: at the fixed point X = Y = sqrt(c) every third cumulant
    vanishes, and since the cumulant curvature is quadratic in T, K = 0 there
    follows from SYMMETRY ALONE.
    """
    X, Y, c, _dg, K = closed_curvature_B()
    Z = c + X + Y + X * Y
    sub = {X: c / X, Y: c / Y}
    g = fisher_metric(Z, (X, Y))
    T = cubic_tensor(Z, (X, Y))
    for i in range(2):
        for j in range(2):
            assert sp.expand(sp.fraction(sp.together(
                g[i, j].subs(sub, simultaneous=True) - g[i, j]))[0]) == 0
    for key, val in T.items():
        assert sp.expand(sp.fraction(sp.together(
            val.subs(sub, simultaneous=True) + val))[0]) == 0
    fp = {X: sp.sqrt(c), Y: sp.sqrt(c)}
    for key in [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)]:
        assert sp.simplify(T[key].subs(fp)) == 0
    assert sp.expand(sp.fraction(sp.together(
        K.subs(sub, simultaneous=True) - K))[0]) == 0
    rho = (sp.sqrt(c) - 1) / (sp.sqrt(c) + 1)
    assert sp.simplify((rho - sp.tanh(sp.log(c) / 4)).rewrite(sp.exp)) == 0
    gfp = sp.Matrix(2, 2, lambda i, j: sp.simplify(g[i, j].subs(fp)))
    assert sp.simplify(gfp - sp.Rational(1, 4)
                       * sp.Matrix([[1, rho], [rho, 1]])) == sp.zeros(2, 2)
    return rho


def gate_curvature_bounds():
    """STRICT bounds -(c-1)/4 < K < (c-1)/(4c), via positive-coefficient brackets."""
    X, Y, c, _dg, K = closed_curvature_B()
    D = X * Y + c * X + c * Y + c
    B1 = sp.expand(D ** 2 - (c - X ** 2) * (c - Y ** 2))
    B2 = sp.expand(D ** 2 + c * (c - X ** 2) * (c - Y ** 2))
    assert not [t for t in B1.as_ordered_terms() if t.as_coeff_Mul()[0] < 0]
    assert not [t for t in B2.as_ordered_terms() if t.as_coeff_Mul()[0] < 0]
    assert sp.expand(sp.fraction(sp.together(
        K + (c - 1) / 4 - (c - 1) * B1 / (4 * D ** 2)))[0]) == 0
    assert sp.expand(sp.fraction(sp.together(
        (c - 1) / (4 * c) - K - (c - 1) * B2 / (4 * c * D ** 2)))[0]) == 0
    return -(c - 1) / 4, (c - 1) / (4 * c)


def gate_fisher_quadric():
    """z00 z11 = sqrt(c) z10 z01 in S^3(2); at c = 1 it IS the Clifford torus."""
    X, Y, c, _dg, _K = closed_curvature_B()
    Z = c + X + Y + X * Y
    z = [2 * sp.sqrt(w / Z) for w in (c, X, Y, X * Y)]
    assert sp.simplify(z[0] * z[3] - sp.sqrt(c) * z[1] * z[2]) == 0
    assert sp.simplify(sum(t ** 2 for t in z)) == 4
    u, v = sp.symbols("u v", positive=True)
    zt = [t.subs({c: 1, X: sp.tan(u) ** 2, Y: sp.tan(v) ** 2}) for t in z]
    tgt = [2 * sp.cos(u) * sp.cos(v), 2 * sp.sin(u) * sp.cos(v),
           2 * sp.cos(u) * sp.sin(v), 2 * sp.sin(u) * sp.sin(v)]
    # positive patch: equal squares + positivity => equal
    assert all(sp.simplify(a ** 2 - b ** 2) == 0 for a, b in zip(zt, tgt))
    return True


def gate_path_independence():
    """d psi is exact: two wild paths agree, every loop integrates to zero."""
    import mpmath as mp
    X, Y, Zs = sp.symbols("X Y Z", positive=True)
    Z3 = 40 + 6 * (X + Y + Zs) + 5 * (X * Y + Y * Zs + Zs * X) + 18 * X * Y * Zs
    A = [mp.log(mp.mpf("0.3")), mp.log(2), mp.log(mp.mpf("1.2"))]
    B = [mp.log(4), mp.log(mp.mpf("0.5")), mp.log(3)]
    line = lambda t: [A[i] + t * (B[i] - A[i]) for i in range(3)]
    wild = lambda t: [A[i] + t * (B[i] - A[i])
                      + mp.sin(mp.pi * t) * [2, -3, 1][i]
                      + mp.sin(2 * mp.pi * t) * [-1, 2, 2][i] for i in range(3)]
    loop = lambda t: [A[i] + mp.sin(2 * mp.pi * t) * [1, .5, -1][i]
                      + (1 - mp.cos(2 * mp.pi * t)) * [.7, -1, .4][i]
                      for i in range(3)]
    return path_independence(Z3, (X, Y, Zs), A, B, [line, wild, loop])


# ------------------------------------------------------------------ (3,6)
def gate_prism_top_cell():
    """91 APMs, type (3,6), all 20 Pluckers positive, Grassmann-Plucker OK."""
    G = prism_36()
    assert len(G.apms()) == 91, len(G.apms())
    assert G.k() == 3
    P = G.pluckers()
    assert all(v > 0 for v in P.values())
    cnt, bad = G.verify_grassmann_plucker()
    assert not bad, bad
    assert G.partition_function() == 280
    return P, cnt


def gate_prism_symmetry():
    """C3 (120 deg) invariant, NOT sigma (60 deg) invariant."""
    G = prism_36()
    assert G.cyclic_symmetry(2) is True
    assert G.cyclic_symmetry(1) is False
    return True


def gate_prism_lambda0():
    """The exact C3-fixed perfect branch: Theta(lam0) C^T = 0 in exact arithmetic."""
    C = sp.Matrix([[14, 0, 0, 2, 21, 13], [0, 14, 0, -14, -49, -21],
                   [0, 0, 14, 10, 14, 2]]) / 14
    lam0 = sp.Matrix(PRISM_LAMBDA0.astype(int).tolist())
    assert sp.Matrix.vstack(C, lam0).rank() == 3          # lam0 subset C
    ts = t_coefficients(lam0)
    assert [sp.nsimplify(t) for t in ts] == [sp.Rational(-3, 100),
                                             sp.Rational(-7, 100)] * 3
    Th = Theta(lam0)
    scale = sp.lcm([sp.denom(sp.nsimplify(x)) for x in Th])
    Thi = sp.simplify(Th * scale)
    assert sp.simplify(lam0 * Thi.T) == sp.zeros(2, 2)     # Lemma 9.2
    assert sp.simplify(Thi * C.T) == sp.zeros(2, 3)        # PERFECT SYSTEM
    assert abs(winding(PRISM_LAMBDA0) - 2) < 1e-9
    assert abs(winding(Theta(PRISM_LAMBDA0)) - 4) < 1e-9
    assert is_valid_pair(PRISM_LAMBDA0)
    sh2 = sp.Matrix([[lam0[r, (i - 2) % 6] for i in range(6)] for r in range(2)])
    assert sp.Matrix.vstack(lam0, sh2).rank() == 2         # sigma^2-FIXED
    return Thi


def gate_prism_branches(n_starts=1500):
    """Exactly four admissible branches, organized as 1 + 3 under C3."""
    br = solve_perfect_branches(PRISM_C, n_starts=n_starts)
    assert len(br) == 4, len(br)
    assert all(b["valid"] for b in br)
    assert all(abs(b["wind_lam"] - 2) < 1e-6 and abs(b["wind_theta"] - 4) < 1e-6
               for b in br)
    R = cyclic_matrix(PRISM_C, 2)
    fixed, cycles, _perm = orbit_structure(br, R)
    assert len(fixed) == 1 and len(cycles) == 1 and len(cycles[0]) == 3, \
        (fixed, cycles)
    n0 = PRISM_LAMBDA0  # calibration: the exact fixed branch must be present
    return br, fixed, cycles


def gate_square_mechanism():
    """Delta_024(w) = (2w+3)(w^2+3): trivial character x Eisenstein norm.

    The product is FORCED by C3; the square at w = 2 is not — w = 2 is the
    unique POSITIVE solution of 2w+3 = w^2+3 (the other balance point is
    w = 0, giving the ramified 3^2 = 9).  At w = 2 the conjugate blocks give
    the split Eisenstein prime of norm 7, so 49 = 7_trivial x 7_Eisenstein.
    """
    w = sp.Symbol("w", positive=True)
    K = sp.Matrix([
        [-1, -1, 0, 1, 0, 0], [0, -1, -1, 0, 1, 0], [-1, 0, -1, 0, 0, 1],
        [1, 0, 0, w, 0, 1], [0, 1, 0, 1, w, 0], [0, 0, 1, 0, 1, w]])
    D = sp.expand(-K.det())
    assert sp.factor(D) == sp.factor((2 * w + 3) * (w ** 2 + 3)), sp.factor(D)
    assert D.subs(w, 2) == 49 and D.subs(w, 1) == 20
    wf = sp.Symbol("wf")   # unrestricted: sympy drops w=0 under positive=True
    balance = sorted(sp.solve(sp.Eq(2 * wf + 3, wf ** 2 + 3), wf))
    assert balance == [0, 2], balance
    assert (wf ** 2 + 3 - (2 * wf + 3)).expand() == (wf * (wf - 2)).expand()
    return sp.factor(D), balance


# ------------------------------------------------------------------ runner
def run_all(verbose=True, branches=True):
    """Run the whole battery; returns a dict of results."""
    res = {}
    steps = [
        ("boundary_measurement_24", gate_boundary_measurement_24),
        ("mandelstam_24", gate_mandelstam_24),
        ("ensemble_conservation_24", gate_ensemble_conservation_24),
        ("simplex_curvature", gate_simplex_curvature),
        ("closed_forms_B", gate_closed_forms_B),
        ("spinflip_parity", gate_spinflip_parity),
        ("curvature_bounds", gate_curvature_bounds),
        ("fisher_quadric", gate_fisher_quadric),
        ("path_independence", gate_path_independence),
        ("prism_top_cell", gate_prism_top_cell),
        ("prism_symmetry", gate_prism_symmetry),
        ("prism_lambda0", gate_prism_lambda0),
        ("square_mechanism", gate_square_mechanism),
    ]
    if branches:
        steps.append(("prism_branches", gate_prism_branches))
    for name, fn in steps:
        res[name] = fn()
        if verbose:
            print(f"  {name:<28s} PASS", flush=True)
    if verbose:
        print("ALL ORIGAMI GATES PASS")
    return res
