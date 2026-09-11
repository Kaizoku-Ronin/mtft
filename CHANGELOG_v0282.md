# mtft v0.28.2 — SPEC-01: hyperbolic Laplace spectrum on X0(N) (2026-09-12)

`surface.spectral`: cusp-truncated P1 finite elements on the glued Manin ideal triangles in the
upper-half-plane chart (stiffness = Euclidean by conformal invariance; mass carries 1/y²), with
truncation heights scaled by cusp width and gluing keyed by log-height.  `maass_candidates`
separates Y-stable Maass cusp-form candidates from migrating Eisenstein pseudo-modes.
Results: Gamma_0(11) stable {4.40, 6.47, 9.15, 9.54, 11.13, …} (mesh-converged to 0.3%);
X0(143) stable {0.392, 0.451, 0.563, 0.972, 1.007, 1.011, 1.069, 1.202, 1.204, 1.245, …}, so
lambda_1 ≈ 0.39 ≥ 975/4096 (Kim–Sarnak, proven) and ≥ 1/4 (Selberg).  Near-degenerate stable pairs
are the oldform signature.  Recorded design error: the first prototype truncated all cusps at one
height and produced no Y-stable modes.  N = 143 numbers are DIAGNOSTIC at h = 0.2 pending a finer run
(memory-bound here); N = 11 is CERTIFIED_NUMERICAL.  These eigenvalues are the KK masses of the
fiber in curvature units; no physical scale.  Also carries KK-01 (see KK01/ report): anomaly-free
chiral flux spectra on X0(143) are t·(g−1)·(5,−4,1) on charges (1,2,3); the (5,−4,1) is arithmetic
of the integers, the unit g−1 = 12 and the Hecke flavour multiplets are the curve's.
Base: live 0.28.0 + unpushed 0.28.1.  Tests: `tests/test_surface_spectral.py` (1 fast, 1 slow).
