"""TF-01 (2026-10-05): can any jet configuration host the top?  No — and what to build instead.

1. The minimal-breaking theorem (EXACT, local characters).  Net jet families of a block with exceptional charge e at the point are the degree-(|e| - 1) coboundary
   jets, character |e| - 1.  The Pati–Salam condensate lives in the (4bar, 1, 2) block together with u^c, so its jets have degree at most |e_u| - 1 and its density
   reaches character -(|e_u| - 1) at best.  A quark Yukawa therefore needs at least u = d_Q + d_u + 2 - (|e_u| - 1) = |e_Q| + 1 >= 2 units of rotation breaking,
   whatever the degrees: the standard point (|e| = 3) needs 4 (YA-1), a degree-0 "constant-jet" top at an |e| = 1 point still needs 2.  Global condensate
   components can supply more negative character only at the price of (r/L)^j.  No jet configuration gives an unsuppressed Yukawa.
2. Two units cost (NUMERICAL).  The torus's two-fold anisotropy is the non-holomorphic weight-2 Eisenstein series: E2-hat(tau) = E2 - 3/(pi Im tau) = 0.106 for 143a1
   (zero for square and hexagonal lattices — the lattice is nearly balanced).  At the blow-up radius 0.326 the two-unit breaking is 0.019 (edge) to 0.076 (across):
   even the best jet configuration gives y <~ 0.08 against y_t = 0.70, nine times short.  With LB-01 this closes every jet realisation of the top.
3. The two-plus-one design (STRUCTURAL).  Put two families as degree-1 jets (|e| = 2 at the blown-up cusp: the spin-1/2 doublet x, w) and the third as a flux
   family: a block with curve degree a and torus flux b has global chiral index a b, so (a, b) = (1, 1) on Q and (-1, -1) on u^c gives one global family each.
   Then b_H = -(b_Q + b_u) = 0 keeps the Higgs torus-constant and Phi-type, and the flux family's Yukawa is a global cup product — a magnetised-torus overlap,
   not a point-localised amplitude — so no character obstruction applies.  The jet families need u = |e_Q| + 1 = 3 units, the right ballpark for the second
   family.  To check: anomalies, D-flatness with torus flux, PA-01 stability, the jet–flux mixing (V_cb ~ 0.04).
OPEN: building and auditing the two-plus-one background."""
import math
import mpmath as mp


def min_breaking_units(e_Q, e_u):
    d_Q, d_u = abs(e_Q) - 1, abs(e_u) - 1
    c_min = -(abs(e_u) - 1)
    return d_Q + d_u + 2 + c_min

def two_fold_anisotropy(dps=30):
    from mtft.research.theta_torus import tau_143a1
    from mtft.matter.lattice_breaking import blowup_radius
    mp.mp.dps = dps
    tau = mp.mpc(tau_143a1(dps)); q = mp.exp(2j * mp.pi * tau)
    E2 = 1 - 24 * mp.fsum(mp.mpf(sum(d for d in range(1, n + 1) if n % d == 0)) * q ** n for n in range(1, 60))
    E2hat = mp.re(E2) - 3 / (mp.pi * tau.imag)
    G2hat = (mp.pi ** 2 / 3) * E2hat * tau.imag                   # unit-area lattice
    r = blowup_radius()
    return {"E2hat": float(E2hat), "G2hat": float(G2hat), "delta2_edge": float(abs(G2hat) * r ** 2 / 2), "delta2_across": float(abs(G2hat) * (2 * r) ** 2 / 2)}

def flux_family_design(a_Q=1, b_Q=1):
    a_u, b_u = -a_Q, -b_Q
    return {"Q": (a_Q, b_Q, a_Q * b_Q), "u_c": (a_u, b_u, a_u * b_u), "higgs_torus_flux": -(b_Q + b_u),
            "jet_families_at_point": 2, "jet_exceptional_charge": 2, "jet_breaking_units": min_breaking_units(2, 2)}

def report():
    from mtft.matter.lattice_breaking import top_yukawa_at_crossing
    tf = two_fold_anisotropy(); yt = top_yukawa_at_crossing()["y_t"]
    return {"units": {e: min_breaking_units(e, e) for e in (1, 2, 3)}, "two_fold": tf, "best_jet_bound": tf["delta2_across"],
            "shortfall": yt / tf["delta2_across"], "design": flux_family_design()}
