# KK-TOWER-01 — magnetic Bochner spectra of the M1 bundles on X0(143) in the compact metric (2026-09-14; DIAGNOSTIC)

## Construction (`surface.magnetic`)
P1 finite elements with Peierls phases.  The Chern connection of O(D) with its HYM metric is coexact in the unitary
frame ((1/2)*d log h_D), so on the mesh it is the dual-graph solution for edge phases with prescribed triangle
fluxes: uniform smooth curvature -2 pi deg(D)/A_c per unit area plus invisible 2 pi n lumps in the triangles carrying
the divisor points (cusp points: a boundary-adjacent triangle).  No Green's functions are needed — the lump solve IS
the Green's function, without the P1-gradient error that broke a first attempt.  Stiffness conformally invariant;
mass matrix with e^{2u}; complex Hermitian shift-invert eigensolve.

## Validation: holomorphic sections appear as a Landau level at 2 pi deg/A_c with multiplicity h^0
   degree 30 (h^0 = 18 for every bundle):  18 modes at 1.2525 vs 1.2452 (+0.6%), spread 1.4%, gap 0.50 to the 19th
   degree 15 (h^0 >= 3 for every bundle):  3 modes at 0.627 vs 0.623 (+0.7%), spread 0.3%
   O(p1+p2+p3), degree 3 (h^0 = 1):        1 mode at 0.1263 vs 0.1245
   O(-p1-...-p6), degree -6:               1 mode at 0.2504 (the section of the conjugate bundle)
The mesh reproduces Riemann–Roch.  The +0.7% offset is the h = 0.2 discretisation (h = 0.3: +1.2%).

## Errors made and caught on the way (all in the register)
1. Face-constant gradients have zero curl: the flux lives on the edges (averaged normal jump) — caught by the flux check.
2. Point flux spread over a star of triangles is visible; it must sit in one triangle (invisible mod 2 pi).
3. Cusp flux entering through the truncation circles was absorbed as a uniform shift by a mean subtraction, giving a
   non-integer effective degree — caught by the missing Landau cluster.
4. The Higgs modes of a degree -6 bundle are H^1, i.e. sections of K (x) L^{-1} of degree 30 — the KK operator lives there.
5. The family sections live in S0 (x) O(sum P), degree 15, not "degree 3" (that is the Dirac index).
6. A generic degree-15 bundle has h^0 = 3 (Riemann–Roch), not 1 — the control was right and my expectation wrong.

## The M1 towers (curvature units: the compact metric has K = -1, area 48 pi; h = 0.2)
   family bundle S0(sum P), degree 15:  0.6255, 0.6261, 0.6298 | 0.6498, 0.6862, 0.7420, 0.8206, 0.9011, ...
   spin structure S0, degree 12:        0.5018, 0.5138 | 0.5619, 0.6111, 0.6529, 0.6943, ...
   Higgs bundle K(2 sum P), degree 30:  18 modes 1.238–1.260 | 1.3037, 1.8044, 1.8971, 2.0052, ...
Spinor and gauge-field KK masses follow by constant Weitzenböck shifts (constant curvature, constant flux), so the
GAPS are the physical content: the first family-tower excitation lies only 0.023 above the zero-mode level
(control bundle: 0.066), and the Higgs tower has a 19th mode 0.05 above the Landau level before the true gap of 0.5.
Both near-degenerate states are DIAGNOSTIC (single mesh; cusp points modelled at the truncation circle, which
shifts the S0 Landau pair by 1–3%) and are the first things to certify on a second mesh: a light vector-like
fourth family-like state and a 19th light Higgs direction would be physical predictions of M1.

## What this enables (the wave's dependency order)
KK spectra and eigenmodes -> overlaps of zero modes with KK modes -> (i) one-loop Coleman–Weinberg potential on the
Higgs moduli space (vacuum selection), (ii) wavefunction-renormalisation matrix of the families (radiative
hierarchy), (iii) KK decay widths into Standard Model fermions (collider-type predictions).  Eigenvectors are
available from the same solve; the overlap machinery exists.
