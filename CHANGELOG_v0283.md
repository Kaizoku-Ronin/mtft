# mtft v0.28.3 — SPEC-02: arithmetic census of the hyperbolic spectrum (2026-09-12)

`surface.spectral.arithmetic_census`: point evaluation of FEM eigenmodes across faces (reduction to
the coset fundamental domain, chart rotation by the dart index, P1 interpolation), Atkin–Lehner
parities via the explicit W_Q matrices, and the T2 Hecke eigenvalue (u(2z)+u(z/2)+u((z+1)/2))/√2.
Two independent classifiers agree on every mode at N = 11 and N = 143: Y-unstable pseudo-modes
satisfy the Eisenstein formula a2 = 2cos(t log 2), t = √(λ−¼) (1.881/1.875, 1.828/1.814, 1.752/1.740
at 143), Y-stable modes do not.  X0(143) stable modes: λ = 0.392 (−,−) a2 = 0.87; 0.451 (−,−) −0.11;
0.563 (+,+) −1.32; 0.972 (−,+) 0.58; 1.007 (+,+) −1.80; 1.011 (+,+) 1.16; 1.069 (+,+) −1.26;
1.202/1.204 (w13 = −1, W11 mixing: oldform pair); 1.245 (+,−) −0.76; 1.345 (+,−) 1.65.  All |a2| ≤ 2.
Physics reading (overlay): the fiber's KK tower is organized by conserved discrete AL charges and
Hecke flavour eigenvalues; the Eisenstein continuum begins at λ = 1/4 (Selberg threshold) and every
bound state sits above it, distinguished from the continuum by its Hecke eigenvalues.
Status: DIAGNOSTIC at h = 0.2 (parities ±1 to 0.5%, eigenvalues to ~1%); certification needs the
finer run.  Tests: `tests/test_surface_spectral_census.py`.  Base: unpushed 0.28.2.
