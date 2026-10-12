"""Pati–Salam first look (2026-10-02): the minimal fermion flux of Matter-04 regrouped as SU(4) x SU(2)_L x SU(2)_R.  Index counts only.

At the Scherk–Schwarz minimum (exceptional fluxes phi = (c 0, L -3, B 3, A1 0, A2 3), classes (a, 0, e) with a = -e/3 at eps = A_E/3) the stack classes
coincide pairwise: A1 with the colour stack c, A2 with B.  Without gauge Wilson lines separating them, the blocks between coinciding stacks are trivial
bundles, their gauge bosons are massless, and the stack group U(3) x U(2) x U(1)^3 enhances to U(4) x U(2) x U(2):
    4   = c (colour, 3) + A1 (lepton number as the fourth colour)
    2_L = L
    2_R = B (d-type) + A2 (u-type)
Blocks between the groups (each with b = 0, so under the twist only exceptional jets survive, e^2 per rank pair, |e| net):
    (4, 2, 1)       <- Hom(L, 4):   Delta phi = 3   : 4*2*9 = 72 Weyl, net 24 = three families of (q, l)
    (4bar, 1, 2)    <- Hom(4, R):   Delta phi = -3  : 72 Weyl, net 24 = three families of (u^c, d^c, e^c, nu^c)
    (1, 2, 2)       <- Hom(R, L):   Delta phi = -6  : 2*2*36 = 144 Weyl in a real representation (no protected chirality); its bosons hold the
                                                     Matter-03 H_d/H_u pair — the Pati–Salam electroweak bidoublet
Total 288 = the Matter-04 minimum.  Inside each group Delta phi = 0: no fermion zero modes (the twist removes them), massless gauge bosons.
The four right-handed components land exactly where Pati–Salam puts them: (c, B) = d^c, (c, A2) = u^c, (A1, B) = e^c, (A1, A2) = nu^c — a right-handed
neutrino per family, without being asked for.  Gauge Wilson lines on 143a1 that separate A1 from c and A2 from B are flat, so they do not change these
indices: they are the Wilson-line breaking SU(4) x SU(2)_R -> SU(3) x U(1)_{B-L} x U(1)_{3R}.
OPEN (the next investigation): the breaking to SU(3) x SU(2) x U(1)_Y (Wilson lines versus (4,1,2) scalars from the bosonic zero modes), the extra U(1)s and
their Stueckelberg masses, anomalies, the Yukawa structure (4,2,1)(4bar,1,2)(1,2,2) and its mass relations, and the vector-like jet pairs."""
from fractions import Fraction

PHI = {"c": 0, "L": -3, "B": 3, "A1": 0, "A2": 3}
RANKS = {"c": 3, "L": 2, "B": 1, "A1": 1, "A2": 1}
GROUPS = {"4": ("c", "A1"), "2L": ("L",), "2R": ("B", "A2")}

def group_phi():
    out = {}
    for g, stacks in GROUPS.items():
        vals = {PHI[s] for s in stacks}
        assert len(vals) == 1, "stacks in one Pati–Salam group must share their flux"
        out[g] = vals.pop()
    return out

def block_counts():
    """(rep, Delta phi, total Weyl = dim * Delta^2, net = dim * |Delta|) for the three inter-group blocks."""
    gp = group_phi(); dim = {"4": 4, "2L": 2, "2R": 2}
    rows = []
    for (g1, g2, rep) in (("4", "2L", "(4,2,1)"), ("4", "2R", "(4bar,1,2)"), ("2L", "2R", "(1,2,2)")):
        d = gp[g1] - gp[g2]; D = dim[g1] * dim[g2]
        rows.append({"rep": rep, "delta_phi": d, "weyl_total": D * d * d, "net": D * abs(d), "families": Fraction(abs(d) * D, dim[g1] * dim[g2]) if rep != "(1,2,2)" else None})
    return rows

def right_handed_assignment():
    return {("c", "B"): "d^c", ("c", "A2"): "u^c", ("A1", "B"): "e^c", ("A1", "A2"): "nu^c"}

def first_look():
    rows = block_counts()
    return {"groups": GROUPS, "group_phi": group_phi(), "blocks": rows, "total": sum(r["weyl_total"] for r in rows), "right_handed": right_handed_assignment()}
