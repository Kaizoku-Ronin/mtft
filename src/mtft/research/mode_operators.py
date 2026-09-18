"""INT-09 (AXG-03 §§6–7, AXG-04 §3): index versus cohomology, purity certificates, the elementary-scalar Bochner bound,
and the 6D scalar-Yukawa chirality gate.

* Dirac index of S (x) L on a curve: deg L (h^0 - h^1 of S0 (x) L).  Index is not h^0: O(P) gives (h^0, h^1) = (2, 1) — one
  family plus a mirror pair — while O(P + Q - R) with u(P) != u(Q) gives (1, 0): a pure single mode.
* Purity for S0 (x) O(D), D = sum of W13-fixed points: h^1 = h^0(S0(-D)) = #{a + b u vanishing on D} = 2 - rank of the
  evaluation matrix [1, u(P_i)].
* Elementary charged scalar in a degree-d line bundle with constant curvature: Bochner lambda_min >= 2 pi |d| / A; for
  d = -6 on the compact metric (A = 48 pi R^2) this is 1/(4 R^2) — measured on the magnetic mesh as 0.2504–0.2512
  (KK-TOWER-01/02).  Eighteen Dolbeault H^1 classes are NOT eighteen zero modes of this scalar operator (H-19).
* 6D Yukawa: the ordinary Lorentz-scalar bilinear psi_1 psi_2 H exists only for OPPOSITE 6D chiralities (H-20)."""
import sympy as sp

def spin_dirac_index(degree_L): return int(degree_L)

def s0_twist_cohomology(u_values, degree_shift_points):
    """(h^0, h^1) of S0 (x) O(D) for D = (sum of the points with the given u-values) - (points subtracted), using the pencil
    (1, u) and generic evaluation.  u_values: u at the added points; degree_shift_points: number of subtracted generic points.
    h^1(S0(D_+)) = 2 - rank[1, u(P_i)]; h^0 = deg D_+ + h^1 (index); subtracting a generic point lowers h^0 by one while
    h^1 is unchanged if the evaluation at it is nonzero on the section space."""
    E = sp.Matrix([[1, sp.nsimplify(v)] for v in u_values]); h1 = 2 - E.rank(); h0 = len(u_values) + h1
    h0 -= degree_shift_points; return {"h0": h0, "h1": h1, "index": h0 - h1, "pure": h1 == 0}

def purity_certificate(u_values):
    """h^1(S0(D)) = 0 iff the evaluation matrix has rank 2 (both signs of u present on D)."""
    E = sp.Matrix([[1, sp.nsimplify(v)] for v in u_values]); return {"rank": E.rank(), "pure": E.rank() == 2}

def bochner_bound(degree, area=None, R=1):
    """lambda_min >= 2 pi |d| / A for the elementary scalar; default area the compact metric 48 pi R^2."""
    R = sp.nsimplify(R); A = sp.nsimplify(area) if area is not None else 48 * sp.pi * R ** 2
    return {"lambda_min_bound": sp.simplify(2 * sp.pi * abs(sp.Integer(degree)) / A), "area": A, "zero_modes": 0 if degree != 0 else None}

def six_d_scalar_yukawa_gate(chirality_1, chirality_2):
    """Ordinary 6D Lorentz-scalar Yukawa bilinear between two Weyl fermions exists iff their 6D chiralities are opposite."""
    ok = chirality_1 != chirality_2
    return {"allowed": ok, "reason": "opposite 6D chiralities required for the scalar bilinear" if ok else "same-chirality pair: no ordinary scalar bilinear (AXG-03 §7)"}
