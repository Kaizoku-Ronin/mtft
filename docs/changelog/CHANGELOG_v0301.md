# mtft v0.30.1 — step 5 machinery: Riemann–Roch spaces with CM poles (2026-09-13)

- Frozen weight-2 basis extended to 400 coefficients (first 131 identical to the previous freeze); the exact
  Yukawa tensor keeps precision 130 (Sturm-sufficient) via `al_eigenbasis(prec=...)`.
- New `surface.rrspace`: CM classes (20 W143-fixed, 4 W13-fixed, W11-partners explicit), eigenform values at
  points of H via height-maximising reduction, and the three-family space H^0(S0 (x) O(P1+P2+P3)) as G/(f1 f_K)
  with G in H^0(K^2) under 12 exact cusp conditions and 21 point conditions: rank 33 certified (smallest
  singular value 2e-4 vs 1e-10 noise), dimension 3, third section with simple poles at exactly the same-sign
  pair of u — the purity theorem's bookkeeping, reproduced constructively.
- Exact fact recorded and verified (1e-10): div(f1 dz) = the 24 Atkin–Lehner fixed points, f1 the 143a1
  differential; every W13-even differential vanishes at the W13 points, every W143-even one at the W143 points.
- The two tests left stale by 0.29.1 (condensation_energy dict; petersson_gate keys) fixed in-tree to match
  the repaired repo; SM-01 report corrected (nu^c range).  Process: full fast suite run before this candidate;
  built from a clean tree (no __pycache__).
Next: H^0(K(2 SigmaP)) (18-dim) from cubic products with double vanishing, the 3 x 3 x 18 up-type Yukawa tensor
of M1, compact-metric normalisation.  Base: PyPI 0.30.0.  Pin four-way.
