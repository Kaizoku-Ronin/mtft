"""GS-01 (2026-10-05): Green–Schwarz relaxation, the slepton direction, and what protects the vector-like pairs.

1. Relaxation undoes the condensate (EXACT).  Both stack U(1)s are anomalous (PS1-2) and their FIs are Kaehler data, xi = J.F.  If the Kaehler moduli are dynamical,
   V = (1/2)[(xi - v^2)^2 + (5/8) v^4] must be minimised over xi as well as v: the unique minimum is v = 0, xi = 0 (the U(1) term vanishes along v^2 = xi and the
   nonabelian quartic pulls v to zero along it).  The FI-driven condensate of PS-01 exists only if the moduli are stabilised at xi != 0 by something outside the
   D-term sector; with fixed moduli CO-1 (v^2 = (8/13) xi) stands.  PS1-3 becomes a caveat.
2. The slepton direction (EXACT).  The condensate's SU(4) D-term, -v^2 (T^A)_44, is non-anomalous and cannot relax: it gives the (4,2,1) scalars +(1/8) g^2 v^2 in the
   colour directions and -(3/8) g^2 v^2 in the lepton direction, whatever the U(1)s do.  Only a U(1) FI on the Q block, with fixed moduli, can lift the slepton
   doublet (VS-2's xi_4L >= (11/13) xi_R4); if the U(1) D-terms relax, the slepton doublet is tachyonic and the electroweak group breaks at the Pati–Salam scale.
3. What protects the vector-like pairs (STRUCTURAL).  Give the adjoint Phi charge -2 and every A-type mode charge +1 under U(1)_V: every zero-mode cubic Tr Phi[A, A]
   and every D-term is then invariant, and the Dirac mass A A' carries charge 2 — the pairs are massless to all orders in the zero-mode theory (instantons aside).
   The Kaluza–Klein tower breaks it: its supersymmetric masses pair Phi_n with A_n (charge -2 + 1 != 0).  Radiative masses are therefore Kaluza–Klein-mediated,
   m ~ (g^2 C / 16 pi^2) M_SS F(M_SS / M_KK), and with VS-4's split partners the right-handed ones are further suppressed.
OPEN: what stabilises the Kaehler moduli at nonzero FI; the tower's loop function F."""


def d_potential(v2, xi):
    return 0.5 * ((xi - v2) ** 2 + (5 / 8) * v2 * v2)

def minimum_fixed_xi(xi=1.0):
    return (8 / 13) * xi

def minimum_dynamical_xi(n=300, span=2.0):
    best = None
    for i in range(n + 1):
        for j in range(n + 1):
            v2, xi = span * i / n, span * j / n
            val = d_potential(v2, xi)
            if best is None or val < best[0] - 1e-15: best = (val, v2, xi)
    return {"v2": best[1], "xi": best[2], "V": best[0]}

def su4_shift(direction):
    """m^2 / (g^2 v^2) from the condensate's SU(4) D-term for a (4,2,1) scalar."""
    return -(0.75 if direction == "lepton" else -0.25) / 2

def u1v_charges():
    return {"Phi": -2, "A": 1}

def invariant(charges):
    return sum(charges) == 0

def report():
    q = u1v_charges()
    return {"fixed_xi": minimum_fixed_xi(), "dynamical_xi": minimum_dynamical_xi(), "slepton": su4_shift("lepton"), "squark": su4_shift("colour"),
            "cubic_PhiAA": invariant([q["Phi"], q["A"], q["A"]]), "dirac_AA": invariant([q["A"], q["A"]]), "kk_Phi_A": invariant([q["Phi"], q["A"]])}
