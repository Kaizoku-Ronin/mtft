# mtft v0.31.1 — CI fix: rank-tolerance calibration in rrspace (2026-09-12)

`w13_grading_and_texture`: rank tolerances 1e-8 -> 1e-6 in the even/odd Higgs rank gates.
The v0.31.0 gate summed the 8 even Higgs directions with tol = 1e-8 x sigma_1 = 1.6e-8; the
summed third singular value is construction noise (2.2e-9 on one BLAS, >= 1.6e-8 on the CI
runner), so the rank-2 assertion was environment-dependent and failed in CI (publish blocked;
PyPI never received 0.31.0).  Rank 2 for an even Higgs VEV is exact by the W13 selection rule;
the new tolerance sits between the certified forbidden-entry floor (< 1e-6, measured 4e-8) and
the O(1) allowed scale, with >= 700x margins both ways.  Auditor-diagnosed, auditor-applied.
Base: v0.31.0 (tag superseded before publication).  Pin four-way.
