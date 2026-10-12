# AL-01 — Atkin–Lehner symmetry and one family per cusp (2026-10-02)

Module `mtft.matter.cusp_families`; tests `tests/test_al01_cusp_families.py`. Index-level statements under the Scherk–Schwarz twist; every cusp maps to O, so
all blown-up points (cusp, O) lie on the graph of phi.

1. **Full Atkin–Lehner symmetry forbids three families (EXACT).** On the blow-up of all four cusps, a line-bundle background symmetric under the whole group —
   acting on the gauge group through any homomorphism into {id, D, C, DC} — has no Pati–Salam block with three net families: an element acting trivially on
   the gauge group forces even net chirality; one acting through charge conjugation forces a vector-like spectrum. All sixteen homomorphisms checked over
   exceptional charges in [-4, 4]. So three families plus a residue cannot come from the four Atkin–Lehner sectors.
2. **One involution as left–right parity (EXACT).** Blowing up {infinity, W_Q(infinity)} with W_Q acting as D gives the minimal Pati–Salam indices with 216 fermion
   zero modes, but left- and right-handed families sit at different cusps: the tree-level Yukawa vanishes identically.
3. **One family per cusp (EXACT).** Blow up three of the four cusps; L-charge -1 and R-charge +1 at each. (4,2,1) and (4bar,1,2): exactly three net families
   each with no vector-like partners; bidoublet net 6 with 12 pairs; 96 fermion zero modes (288 for one point). D-flat at eps_i = A_E/3; ample iff A_X > A_E.
   The empty cusp can be any of the four (all related by Atkin–Lehner): the background is unique up to isomorphism. Each family's left- and right-handed modes
   share a cusp. [CC-36: the tree-level Yukawas of jet fermions vanish identically; the diagonal-Yukawa / CKM = 1 statement is retracted — see YK-01.]
OPEN: the Higgs sections of this background and their cusp values (the family hierarchy); the CKM-generating corrections; Pati–Salam-breaking gamma-flux brings
back localised vector-like pairs; no (nu^c, N) pairs in the jet sector, so inverse-seesaw partners must be gamma-localised and their Dirac pairing involves the
same separation as mu.
