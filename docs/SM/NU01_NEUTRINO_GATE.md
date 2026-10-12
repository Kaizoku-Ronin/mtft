# NU-01 — the neutrino gate of the Pati–Salam branch (2026-10-02)

Module `mtft.matter.neutrino_gate`; tests `tests/test_nu01_neutrino_gate.py`.

1. **A cheaper Pati–Salam-breaking flux (EXACT).** gamma-charges (1, 1, -1, -1, -1) on (c, A1, L, B, A2) with delta = -8 kappa: traceless, still D-flat with the
   (4,1,2) block alone, no gamma-flux on the bidoublet; 256 gamma-localised extras per orientation against 384 for PS-01's (1, 1, 0, -2, -2).
2. **Neutral inventory (EXACT, index level, per rank-one entry).** Lepton-doublet entry: 3 + 16 against 6 + 16 (net 3; 19 vector-like pairs). nu^c entry:
   3 + 16 against 6 + 16 (net 3, B-L = +1; 19 vector-like (nu^c, N) pairs). Bidoublet entries: 21 against 15, B-L = 0 — the only B-L-neutral neutral
   fermions. No B-L-neutral singlet fermions: the twist removes adjoint fermions.
3. **Localisation (EXACT on the product).** On each fibre {x} x 143a1 a gamma-charge-c block is the flat bundle of c phi(x): the 4c^2 bosonic extras sit on
   phi^-1(E[c]), the 4c^2 fermionic ones on phi^-1({Q : cQ = P_eta}) — disjoint sets. For c = 1 the bosonic set is the four cusps. Classes with disjoint
   supports have zero cup product, so the leading holomorphic Majorana coupling of the gamma-localised N through the gamma-localised condensate vanishes: a
   selection rule, the first structural source of lepton-number suppression in the construction. Jet-sector couplings at the blown-up point are not covered.
4. **Atkin–Lehner at the cusps (EXACT).** 143a is the (+,+) newform (rank 1, root number -1), so phi o W_Q = phi and W_Q x id_E preserves the graph, gamma, K and
   the spin structure. (Z/2)^2 permutes the cusps simply transitively: the four cusp-localised modes of a c = +-1 block carry its regular representation, one
   mode in each sector including (-,-), which holds no weight-2 newform (Paper 37's residue sector). The single-point blow-up and O(infinity) break the symmetry;
   a symmetric background would blow up the whole Atkin–Lehner orbit of (infinity, O).
5. **Selection rules (EXACT).** B-L is gauged in SU(4); the condensate has |Delta(B-L)| = 1; Majorana masses need two insertions (dimension five, ~ v_R^2/M_*).
6. **Requirements (NUMERICAL, one-loop SM running to M_PS = 1.1 x 10^6 GeV; Y_nu = Y_u there).** m_D = 121 GeV, 0.41 GeV, 0.85 MeV. Normal ordering with
   m_3 = 50 meV (Paper 37's scale), m_2 = 8.7 meV: type I needs M_R(3) = 2.9 x 10^14 GeV (excluded in this branch); an inverse seesaw with M = M_PS needs
   mu(3) = 4.3 MeV and mu(2) = 65 GeV — suppressions mu/M of 3.8 x 10^-9 and 5.7 x 10^-5 relative to the dimension-five estimate.
OPEN: the B-L-conserving Dirac pairing M, jet-sector Majorana couplings, PMNS (Y_nu = Y_u mixes little), Majorana versus Dirac, the Atkin–Lehner-symmetric
background.
