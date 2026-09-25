# mtft v0.30.6 — SM-08/09: compact metric, HYM normalisations, certified up-sector hierarchy (2026-09-14)

`surface.hym` (new): Liouville solver for the compact constant-curvature metric of X0(143) on the SPEC mesh
(Gauss–Bonnet gate: 48 pi to 0.2–0.4%), Green's functions of the compact Laplacian at CM points and cusps
(component-pinned SPD Neumann solve; sub-cell log correction), HYM metrics and kinetic Gram matrices by
quadrature with the compact area element, and `up_sector_mass_ratios`: the bilinear-normalised up-type Yukawa
tensor of M1 and its singular-value ratios over Higgs directions.  Two-mesh certification: family Gram spectrum
(7.9e-5, 6.9e-4, 1) to 1%; m2/m3 in 1e-4–1e-3, m1/m3 ~ 1e-8 for every Higgs direction.  Corrections: SM-08's
rank-4 flag (printout artefact), conjugated normalisation (fixed), bordered Neumann solve (replaced).
Tests: `test_surface_hym.py` (2).  Full fast suite before build; clean tree.  Base: PyPI 0.30.0 + 0.30.1–0.30.5.
