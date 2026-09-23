# R2C-02 — the internal-vector Higgs of M1: mode operator, tachyon mass, Yukawa as gauge interaction (2026-09-19)

## EXACT
- Mode operator on the (0,1)-forms valued in L = L_x L_y^{-1} (degree d_L = -6): m^2 = nabla^*nabla_{L (x) Kbar} - 2(2 pi |d_L|/A) + K_gauss.
  nabla^*nabla is the Bochner Laplacian of the degree-30 bundle (the SM-03 Higgs target space); with K = -1 and A = 2 pi (2g-2) the curvature
  terms cancel exactly (2 pi deg K / A = 1 = -K), giving m^2_LLL = -2 pi |d_L|/A = -1/4 in curvature units, multiplicity h^1(L) = 18, all other
  modes massive.  Flat limit: the Nielsen–Olesen/Bachas tachyon -|F_L|.  Same magnitude, opposite sign to the elementary-scalar bound (H-19).
- The Yukawa is the 6D gauge interaction psi-bar Gamma^zbar A_zbar psi: y_ijk = g_6 I_ijk = g_4 sqrt(A) I_ijk (unit-normalised modes,
  g_4^2 = g_6^2/A).  The SM-03 section-product tensors are these overlaps; the equal-chirality vector selection rule (R2C-01) is satisfied by
  the family-preserving assignment.  This supplies the "specified 6D interaction" the v0.31.1 review found missing — for the vector route.

## DIAGNOSTIC (h = 0.25 magnetic mesh; A_c = 151.7)
- Vector-Higgs spectrum: 18 modes at m^2 = -0.243 (exact -0.2485 in mesh units; +0.9% Landau offset), spread 2.5e-2; first massive mode +0.35.
- Dimensionless overlap I~ = sqrt(A) Y: rms entry 0.27; y_t/g_4 = kappa s_max(v): 0.42 / 0.63 / 0.87 (5/50/95% over Higgs directions), max 1.18.
  Observed y_t/g ~ 0.7 near 1e16 GeV lies inside the typical band.  kappa is an O(1) Clifford-convention factor of the reduction (to be pinned).
- The same order holds for the third-family down and lepton Yukawas (anarchic overlaps): tree-level third-family Yukawa unification, which
  phenomenologically requires large tan beta with the two Higgs doublets M1 has.

## Not claimed
The Higgs direction remains unselected (the ratio is a distribution); no vacuum, no absolute scale, no RG running; the irreducible p2 term of
the declared spectrum (R2C-01) still needs a tensor/gravity sector before the parent stands alone; kappa unfixed.
Tools: `research.vector_higgs` (exact identity test; slow FEM test with margins 6% over ~2% measured).


## Addendum (2026-09-23, v0.33.0): pointwise form and sign restriction
The mode operator equals 2 d*d - |B| identically on the (0,1)-forms of a NEGATIVE-degree block, for every Kaehler metric (the curvature cancels
pointwise); the formula as written applies to d < 0 only (CC-30).  `vector_higgs.kodaira_forms`; `V0330_COMPENDIUM_NOTES.md` item 2.
