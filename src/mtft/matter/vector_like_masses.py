"""VM-01 (2026-10-05): what gives the vector-like pairs their masses — not the Pati–Salam condensate, at any order, in our vacuum.

1. The condensate cannot do it holomorphically (EXACT).  An SU(4) invariant built from fundamentals needs (#4 - #4bar) = 0 mod 4.  Sigma psi psi' for a (4bar,1,2)
   pair has (1, 2): forbidden.  Sigma Sigma'-bar psi psi' (2, 2) is allowed but needs the conjugate condensate, and PB-01's D-term vacuum has only one orientation
   (Sigma-bar = 0).  Sigma^4 is SU(4)-allowed but vanishes for a single rank-one condensate (identical bosons, PB-01).  So the condensate gives no holomorphic
   mass at any order.
2. Nothing else in the vacuum does either (EXACT, given Matter-05 and PB-01).  The Dirac mass psi psi' is gauge-invariant by itself; in W = Tr Phi [A, A] it comes
   from adjoint Phi-type vevs, and every adjoint boson has a positive one-loop mass (Matter-05): no vev.  Torus Wilson lines only move the gamma-localised pairs
   (c phi(x) = w0), they do not lift them.  At tree level every vector-like pair is massless, left and right alike.
3. Consequences (STRUCTURAL).  TH-01's reconciliation needs a mass source the model does not yet have.  And massless coloured vector-like states are excluded
   outright (light vector-like quarks): the gamma-charged-matter design (17 pairs per block) needs a mass source regardless, and that source decides the
   left–right asymmetry.
4. The natural candidate and what it implies (ESTIMATE).  A conjugate condensate Sigma'-bar would give m ~ v_R^2 / M_* to right- and left-handed pairs alike, with
   different O(1) Clebsches (two SU(4) channels, 1 and 15).  Reaching m_R ~ M_PS / 5.6 (17 pairs) then needs M_* ~ 5.6 M_PS: the cutoff, or the Kaluza–Klein scale,
   within a factor ~6 of the Pati–Salam scale — where the towers of TH-4 also enter.
OPEN: a second D-flat orientation of the condensate; one-loop (SUSY-breaking) masses with a chirality source; the Clebsch asymmetry."""


def su4_invariant_possible(n4, n4bar):
    """Necessary condition for an SU(4) singlet in a product of fundamentals and antifundamentals."""
    return (n4 - n4bar) % 4 == 0

def operators():
    return {"Sigma psi psi'": su4_invariant_possible(1, 2), "Sigma Sigma'bar psi psi'": su4_invariant_possible(2, 2), "Sigma^4": su4_invariant_possible(0, 4),
            "psi psi' (Dirac)": su4_invariant_possible(1, 1)}

def vacuum_sources():
    """What the established vacuum supplies (PB-01: one orientation; Matter-05: adjoint bosons positive)."""
    return {"conjugate_condensate": False, "rank_one_condensate": True, "adjoint_vev": False}

def tree_level_mass_possible():
    v = vacuum_sources(); ops = operators()
    via_conj = ops["Sigma Sigma'bar psi psi'"] and v["conjugate_condensate"]
    via_sigma4 = ops["Sigma^4"] and not v["rank_one_condensate"]
    return via_conj or via_sigma4 or v["adjoint_vev"]

def cutoff_ratio(pairs=17):
    """M_* / M_PS if m_R = v_R^2 / M_* must equal the TH-01 requirement M_PS / exp(n l / pairs)."""
    import math
    from mtft.matter.thresholds import required_nl
    return math.exp(required_nl(7 / 3) / pairs)

def report():
    return {"operators": operators(), "tree_level_mass": tree_level_mass_possible(), "cutoff_ratio_17": cutoff_ratio(17), "cutoff_ratio_19": cutoff_ratio(19)}
