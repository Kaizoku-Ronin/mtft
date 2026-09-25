# mtft v0.29.0 — the KK physics toolkit: flux sectors, Wilson lines, Yukawa tensors, condensation (2026-09-14)

Consolidates the unpushed 0.28.1–0.28.5 candidates (CC-24, CC-25, sparse spectral assembly,
SPEC-02 arithmetic census, WL-01 Wilson-line spectrum) and adds three GP-free instruments:

- `surface.yukawa` — exact Yukawa tensor H^0(K)⊗H^0(K^3)→H^0(K^4) of the canonical ring from the
  frozen weight-2 basis (`_data/x0143_weight2_basis.json`: mfbasis coefficients to q^130, W11, W13):
  AL eigenbasis by exact projectors, modular ranks, exact rational coordinates; gates 60/84,
  sectors (12,18,17,13)/(24,18,19,23), 0 selection-rule violations in 7,543 entries, generic rank 13.
- `surface.condensation` — tachyon mass² −|d|/(24R²), Atiyah–Bott condensation energy, and the
  extension long exact sequence whose connecting map is the Yukawa matrix: h^0(V)=48, h^1(V)=0,
  sector ranks (1,6,5,1) for AL-invariant classes.
- `surface.petersson` — hyperbolic quadrature on the glued Manin mesh with height-maximising
  reduction (one-step AL lookahead) and q-expansion evaluation; Petersson Gram matrices of the
  differential, cubic and quartic bases; gate: Riemann bilinear relations from frozen periods
  reproduced to 0.5% diagonal / 1.4% off-diagonal (DIAGNOSTIC).

Physics chain now in the package (all tagged): KK-01 anomaly-quantised chiral content t·12·(5,−4,1);
WL-01 Wilson-line mass of the vector-like pair = Hodge norm / area (verified 0.5%); YUK-01 exact
Yukawa tensor; PET-01 normalisations; COND-01 Higgs mechanism = bundle recombination with
h^0(V) = 48.  Not supplied anywhere: a physical R, the gauge coupling, the position of the VEV in
Ext^1.  Sampler note unchanged: `metropolis_sweep` ensembles before 0.28.0 must be rerun (CC-21, CC-23).
Tests: `test_surface_yukawa.py` (exact, 15 s), `test_surface_petersson.py` (slow), plus the 0.28.x tests.
Base: live 0.28.3 (SHA-256 746d7164…).  Pin four-way.
