# SM-09 — certification of the up-sector hierarchy; corrections to the pipeline (2026-09-14)

## Corrections made in this wave (all recorded)
1. The rank-4 flag of SM-08 was a printout artefact: the algebraic Sym²(3) → 18 map has full rank 6 with singular
   values 1.45e4 … 0.027 (spread 2e-6, far above the 1e-12 residual).  The geometry is intact.
2. The normalised tensor of SM-08 was built with a conjugated factor on the second family index, breaking the
   symmetry of the bilinear coupling (0.38).  Both family indices are now orthonormalised linearly; symmetry 1e-14.
3. The bordered Neumann solve for Green's functions was numerically singular on the finer mesh (entries down to
   1e-9 in the constraint row): replaced by a pinned-node SPD solve, one pin per connected component of the mesh
   graph — needed because at h = 0.15 the centroid filter leaves two isolated sliver triangles (mass 6e-7).
4. Sub-cell logarithmic correction of the point-source Green's functions at the CM points implemented; it moves
   m2/m3 by 3% and the Gram spectrum not at all — the singularity is not what drives the hierarchy.

## Two-mesh certification (h = 0.2, nx = 8  vs  h = 0.15, nx = 14; metric, Green's functions and Grams re-solved)
   compact area:            151.37 (0.38%)          151.12 (0.21%)      -> 48 pi = 150.80
   family Gram spectrum:    (7.80e-5, 6.82e-4, 1)   (7.90e-5, 6.89e-4, 1)   agreement 1%
   m2/m3 (5, 50, 95%):      (2.9e-5, 1.6e-4, 1.6e-3)  (5.3e-5, 2.6e-4, 1.7e-3)
   m1/m3 (5, 50, 95%):      (3e-9, 1.5e-8, 1.2e-7)    (3e-9, 1.3e-8, 1.3e-7)
The HYM norms of the three up-sector family wavefunctions span a factor 1e4, and this spread is converged to 1%;
it is the origin of the hierarchy.  Certified statements: for every Higgs direction the up sector has one heavy
family; m2/m3 lies in 1e-4–1e-3 (median 2–3e-4); m1/m3 ~ 1e-8.  The median of m2/m3 moves by 60% between
meshes, reflecting the less-converged 18 x 18 Higgs Gram (condition 3e6): the order of magnitude is certified,
the second digit is not.
Comparison (orientation only): observed m_c/m_t ~ 2–4e-3, m_u/m_t ~ 6e-6.  The pattern matches; the second
family sits a factor ~10 low at the median (the 95th percentile reaches the observed value), the first family
comes out two to three orders of magnitude lighter than observed.  Both statements are at the mercy of the
Higgs direction (a modulus) and of the down-sector-style twists not yet applied to the up sector.

## Tools
`surface.hym`: `assemble`, `solve_liouville`, `solve_neumann` (component-pinned), `greens_functions`,
`quad_with_nodes`, `hyp_dist`, `gauss_bonnet_gate`, and the one-call pipeline `up_sector_mass_ratios(h, nx)`
returning percentiles, the Gram spectrum and the normalised tensor.  Tests: Gauss–Bonnet within 1%; the
component-pinned Neumann solve on a disconnected toy graph.  Next: the same pipeline on the twisted (down,
lepton) sectors, which needs a vectorised Gamma_0-only reduction at the quadrature points.
