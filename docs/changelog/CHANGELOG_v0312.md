# mtft v0.31.2 — CC-26: HYM normalisation contraction fix; retraction of SM-08…14 mass-ratio claims; SM-15 (2026-09-14)

CC-26 (Astra, Atlas v0.31.1 audit V0311-C01): `hym.up_sector_mass_ratios` and the frozen M2 tensors used (L^{-1})^T
instead of L^{-1} in the family normalisation.  Fixed; `hym.normalise_yukawa` is the single implementation and
`hym.normalisation_is_basis_invariant` the gate (tolerance 1e-6 over a measured 1e-8).  `x0143_m2_tensors_h02.npz`
regenerated with the correct contraction and now ships the raw tensors and all Gram matrices (M1 and M2) for
re-audit.  Corrected: no leading-order hierarchy — M1 up m1/m3 ~ 0.26, m2/m3 ~ 0.62; M2 up 0.19, 0.56; M2 down
0.25, 0.61; CKM anarchic (all moduli ~0.54 median).  Retracted: hierarchy classes, the 7.8e-5 "spin-structure
invariant" as physics, CKM typicality, the SM-13 headline.  Unaffected: all algebraic results.
SM-15: `hym.w13_graded_normalisation`, `hym.w13_epsilon_scan`; degeneracy theorem verified in the corrected
normalisation (even VEV -> (1, 1, 0)); a W13-symmetric vacuum implies m_c = m_t — the W13 grading cannot be the
origin of the quark hierarchy.  `rrspace.w13_grading_and_texture` now also returns the eigenbases.
Registers: outputs/SM/CC26_NORMALISATION_RETRACTION.md, SM15_W13_DEGENERACY_REPORT.md.  Base: PyPI 0.31.1.
KK-TOWER-01: `surface.magnetic` — magnetic Bochner Laplacian of O(D) on X0(143) in the compact metric (P1 Peierls FEM,
coexact connection from the dual-graph solve with invisible divisor lumps); validated against Riemann–Roch (degree 30:
18-fold Landau level +0.6%; degree 15: 3-fold).  M1 towers recorded (DIAGNOSTIC): near-degenerate 4th family-tower
state (+0.023) and 19th Higgs mode (+0.05).  Slow test `test_surface_magnetic.py`.  Register: outputs/SM/KKT01_MAGNETIC_TOWERS_REPORT.md.
