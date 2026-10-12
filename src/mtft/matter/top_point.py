"""TP-01 (2026-10-05): the flux family fails, the constant-jet top survives — a two-point design.

1. Flux families come in twenty-fours (EXACT).  For a block with curve degree a and torus flux b, chi(V*) - chi(V) = 24 b (checked on the surface's intersection
   form: (11, 1) gives -1 and 23), so one unit of torus flux yields 2g - 2 = 24 chiral families — one A-type and 23 Phi-type.  TF-3's flux family is retracted
   (CC-45); chiral families come only from exceptional charges.
2. The two units of a constant-jet top can be bought globally (EXACT structure, ESTIMATE size).  TF-1 priced the condensate's density with jets only.  Its global
   (untwisted) components have Taylor coefficients of every order at the point, and the conj of the second one carries character -2 at the price of their size:
   (r/L)^2 to (2r/L)^2, i.e. 0.11-0.42 at eps = A_E/3, 0.16-0.64 at A_E/2, 0.70 at 0.55 A_E — not the lattice's two-fold anisotropy (0.02-0.08).  So a
   degree-0 top at a point of exceptional charge 1 is marginal, not excluded; TF-2 and LB-4 become caveats (CC-46).
3. The two-point design (STRUCTURAL).  P1 = (cusp, O) with exceptional charge 2: jets(-2) - jets(2) = 3 - 1 = 2 degree-1 families (the spin-1/2 doublet x, w).
   P2 = (cusp, P), a second point on the same cusp fibre (P a torsion point, so the separation is a torus distance, not CC-37's cusp distance), with exceptional
   charge 1: jets(-1) - jets(1) = 1 constant-jet family.  Per block 2 + 1 = 3; the coupling Higgs has exceptional charges (+4, +2), so it stays a torus-constant
   Phi-type pullback.  The jets at P1 need |e_Q| + 1 = 3 units; the top at P2 needs 2, paid globally.
4. A unique Higgs (EXACT).  The Phi-type Higgs space is H^0(C, K_C (x) L_H) with deg = 24 - (a_Q + a_u): empty above 24, thirteen cusp forms at 0 (Matter-04, trivial L_H; a generic flat L_H gives twelve),
   and exactly one — constant along the curve — when a_Q + a_u = 24 and L_Q (x) L_u = K_C.  With b = 0 the curve degrees are free, so the design can carry a
   unique Higgs; D-flatness of the Q block then ties a_Q A_E to 2 eps_1 + eps_2.
OPEN: the jet–top mixing across the torus (V_cb ~ 0.04), the second family's three units, PA-01 stability with two points, the global fraction of the condensate."""
import math


def flux_family_chirality(a, b):
    from mtft.surface.blowup import chi
    return chi((-a, -b, 0, 0)) - chi((a, b, 0, 0))

def global_cost(eps_over_AE, across=True):
    return (4 if across else 1) * eps_over_AE / math.pi

def two_point_families():
    from mtft.matter.background import _jets
    return {"P1": _jets(-2) - _jets(2), "P2": _jets(-1) - _jets(1), "higgs_exceptional": (4, 2)}

def higgs_count(sum_a):
    """Generic h^0(C, K_C (x) L_H), deg K L_H = 24 - sum_a, for the coupling Higgs; '1 iff L_Q L_u = K_C' at sum_a = 24."""
    deg = 24 - sum_a
    if deg < 0: return 0
    if deg == 0: return 1
    if deg == 24: return 13                     # L_H trivial (Matter-04's curve bundles): the thirteen cusp forms; a generic flat L_H gives 12
    return max(deg - 12, 0) if deg < 24 else deg - 12

def d_flat_relation(a_Q, A_E, eps1, eps2):
    """FI of the Q block: a_Q A_E - (2 eps1 + eps2) must vanish without a squark vev."""
    return a_Q * A_E - (2 * eps1 + eps2)

def report():
    return {"flux_chirality": {(a, b): flux_family_chirality(a, b) for (a, b) in ((11, 1), (13, -1), (0, 1))},
            "global_cost": {f: (round(global_cost(f, False), 3), round(global_cost(f), 3)) for f in (1 / 3, 1 / 2, 0.55)},
            "families": two_point_families(), "higgs_count": {s: higgs_count(s) for s in (0, 12, 24, 30)}}
