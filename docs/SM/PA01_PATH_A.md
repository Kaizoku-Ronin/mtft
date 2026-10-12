# PA-01 — Path A: the matter-untwisted completion is a one-loop maximum (2026-10-02)

Module `mtft.matter.path_a`; tests `tests/test_pa01_path_a.py`.

1. **Content (EXACT).** Per matter entry on the three-cusp background: R = 13 A-type (H^1(O(-P)) x 1) + 12 Phi-type (H^0(K(-P)) x dz); Rbar = 12 A-type (H^1(O(P)) x 1)
   + 1 A-type (H^0(O(P)) x dzbar) + 2 jets + 13 Phi-type (H^0(K(P)) x dz): 25 vector-like pairs and the net 3.
2. **Tree-level masses (EXACT).** All are zero modes. A gauge-invariant vector-like mass needs a vev of a stack-U(1) modulus (Phi-type or Wilson-line); the Higgs links
   (4,2,1) to (4bar,1,2) and the Pati–Salam condensate links (4,2,1) only to the 48 bidoublet jets.
3. **Stability (EXACT at one loop, product approximation).** The untwisting Wilson line is the point theta_eta of the relative U(1) T = (1 on the quark–lepton group,
   -1 on L and R). Only the two matter blocks carry T-charge (2), so V(theta) = -2c sum_{odd w} K(w) cos(q theta.w): minimal at theta = 0 (the jet vacuum),
   maximal at theta_eta. The Hessian there is exactly minus the one at the jet vacuum; V(1 - s) = -V(s) along the straight path. The matter-untwisted background
   rolls back to the jet vacuum.
4. **Convergence with Path B (EXACT).** In the jet vacuum the 25 pairs per entry sit at the Scherk–Schwarz scale M = 2 pi dist(P_eta, O)/A_E — M sqrt(A_E) = 3.106 (real
   spin structure), 3.537 (complex pair) — and the 3 net jets stay massless. They are the messengers through which the light families get Yukawas:
   y_eff ~ y_1 y_2 v_R / M.
OPEN: the couplings jet–messenger–Higgs and messenger–jet–condensate (the effective Yukawa matrices); or a rigid realisation of the twist that removes the T Wilson
line (almost-commuting holonomies, orbifold projection).
