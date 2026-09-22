# R2C-08 — 8D chirality selection rule on X0(143) x 143a1; CC-31 (2026-09-19; EXACT)

## CC-31 (my error in R2C-07)
The register claimed a factorised Yukawa Y = I x T for two families of type (3, -1).  Both such families have the same Kuenneth parity
(H^0(X) (x) H^1(E)), and the gauge-vertex Yukawa is a (0,2)-form int psi_1 ^ psi_2 ^ A_H valued in K_S: it needs q_1 + q_2 + q_H = 2 with
q_H = 1, i.e. OPPOSITE parities.  The (1,0) Higgs polarisation cannot help (its form would be valued in K^2).  The naive M4 has no Yukawa.

## Selection rule (EXACT)
Two families have a gauge-vertex Yukawa iff their internal (Kuenneth) chiralities differ; the two blocks then carry the same 8D chirality,
exactly as R2C-01 found in 6D.  `product_surface.yukawa_selection_rule`.

## Exhaustive scan of family types with |index| = 3
- curve x curve (both multiplicities from the three CM points): 8 opposite-parity pairs, every Higgs block has a zero degree on one factor
  ((+-6, 0) or (0, +-2)): NEVER slope-free.  The arithmetic families cannot both be Yukawa partners of a light Higgs on the product surface.
- mixed origin: 8 solutions (4 up to conjugation): Q = (3, +-1) from the curve, u^c = (+-1, -3) from the torus (three theta functions of degree 3
  on 143a1), Higgs (-4, 2) at A_X = 2 A_E or (-2, 4) at A_X = A_E/2.  Then everything of R2C-07 holds: polystable locus, massless Higgs modes
  (56 for (-2, 4)), electroweak scale = Kaehler deviation.

## Flavour consequence (EXACT structure)
With Q's family index on the curve and u^c's on the torus, M_ij(v) = (A v B^T)_ij with A the curve overlaps and B the torus theta triple
products: rank M <= rank v.  One Higgs VEV -> rank one -> only the top massive at leading order; the light masses are set by the
subleading singular values of the VEV matrix.  This is the first hierarchy MECHANISM in the programme (the curve alone gave anarchy);
its size is not predicted until the VEV structure is.  The torus factor is exactly computable (theta functions with characteristics on 143a1).

## Status of M4 after R2C-08
Record updated: naive M4 (both (3,-1)) — no Yukawa; M4' (mixed origin) — Yukawa, slope-free Higgs, rank-1 leading order; OPEN: 8D anomalies,
moduli stabilisation near the locus, the torus theta overlaps (next exact computation), the origin of u^c's multiplicity as arithmetic (level-3
theta functions of the conductor-143 curve rather than the CM points).
