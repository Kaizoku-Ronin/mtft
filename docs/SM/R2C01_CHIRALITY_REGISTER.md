# R2C-01 — 6D chirality assignments for M1: the elementary-scalar Higgs is excluded, the internal-vector Higgs survives
(Astra, 2026-09-18; independently reproduced 2026-09-19; EXACT, conditional on the declared class)

Class: M1's ten oriented bifundamentals with one 6D Weyl field each, original flux m = (0,-3,3,3,0), convention n_L - n_R = eps d.
Selection rule (Spin(6,C) = SL_4(C)): scalar bilinears need opposite chiralities (4 (x) 4^* contains 1); an odd Clifford insertion
needs equal chiralities.  Ranks (own gamma matrices): scalar/transpose 0 for equal, 4 for opposite; vector 4 for equal, 0 for opposite.

Exhaustive scan (own enumeration): 1024 assignments; 16 preserve the six family indices (forcing the six family chiralities to +1);
4 preserve the full ledger; 64 permit all four scalar Yukawa contractions (eps_cL = -eps_ca = -eps_cb, eps_Ld = -eps_ad = -eps_bd);
0 do both — the up-quark pair (cL, ca) alone is the contradiction; 8 of the 16 keep eps_La = eps_Lb (the 4D doublet balance).
Irreducible gravitational term for the four ledger-preserving assignments: [p2] I8 = -n_grav/1440 with n_grav = 24, 22, 18, 16.
Pushforward identity: sum_k m_k dI8/df_k = sum eps (m_i - m_j) I6 — three internal modes multiply the REDUCED anomaly, not the bulk count.

Consequences for route 2:
- An elementary 6D scalar Higgs cannot couple the three arithmetic families of M1.  The Higgs must be an internal vector (gauge)
  component — the flux-recombination picture (KK-05 tachyon, SM-08 Higgs as H^1 of L_x L_y^{-1}) — for which the selection rule
  is satisfied with all family chiralities equal.  OPEN: the action, its physical mode operator, anomaly cancellation (the
  irreducible p2 term of the declared spectrum), and a stable background.
- The gate `six_d_yukawa_chirality` in `research.pipeline` now records this (FAIL for scalar Higgs, OPEN for vector Higgs).
- Not withdrawn: the three-mode arithmetic construction (h^0 = 3 pure), the section-product tensors, the M1 record.
Tools: `research.chirality`; tests `tests/test_route2_r2c01.py` (counts and ranks recomputed, not asserted from the study).
