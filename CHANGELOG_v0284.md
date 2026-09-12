# mtft v0.28.4 — SPEC corrections from Astra's KK-01/KK-03 audit (2026-09-13)

`surface.spectral`: sparse per-face scatter (the dense ravel of all faces cost ~GB; payload now MB),
enabling finer meshes (Astra: N=143 h=0.15 lambda_1 = 0.3910 vs 0.3923 at h=0.2; 0.02% between
truncation heights).  Tag corrected to DIAGNOSTIC everywhere (KK-A07: no continuum enclosure);
scope stated (KK-A06: untwisted scalar Laplacian on the cusped Y0(N); scalar m² in radius units;
not the compact twisted Dirac operator).  Interpolation fallbacks are counted.  KK-01 report carries
an append-only correction appendix (A01 false |j|<=2 claim; A02 122 zero modes; A04 AL sectors
(12,18,17,13) replace the wrong Hecke-multiplet claim).  KK-04 (Astra) recorded: compact twisted Dirac
D_B with exact levels K³: ±2/R (36), ±√6/R (13); S⁵: ±√3/R (24), ±2/R (s = h⁰(S)); O: ±√λ_j; first
nonzero mode from O with λ₁R² <= 4/3 (HYY via the genus-1 AL quotient, degree 8) — verified
arithmetic; AL-fixed theta characteristics 4096 (3072 even, 1024 odd).  No physical R.
Base: live 0.28.3 (SHA-256 746d7164…).  Pin four-way.
