# Replay the T2 memory approximation experiment

Start with `MTFT_T2_Memory_Approximation_Report.md` and the exact-identity note
`TRUNCATION_IDENTITIES.md`. The research outcome is **no validated compression
for the prescribed joint response/time budgets**. A passing replay verifies
that this result and its controls reproduce; it does not reverse that outcome.

With Python 3.12 and the listed scientific dependencies:

```bash
python -m pip install -r requirements.txt
python run_study.py
```

The runner verifies the frozen input and protocol hashes, computes training
selections, hashes them before a separate held-out invocation, runs the
independent implementation and symbolic controls, renders the figure, and
writes `VALIDATION.json` and `run_log.txt`. The two primary phases can also
be run separately:

```bash
python memory_approximation.py --phase train
python memory_approximation.py --phase test
```

Replay overwrites result files beside the scripts. Work in a copy to preserve
the original measurements and checksum manifest. Last-place floating-point
values, archive metadata, paths, and timings may vary. Numerical cross-route
and exact-control comparisons use tolerance 1e-9; approximation decisions
use the separate frozen 1%, 5%, and 10% budgets.

`inputs/T2_blocks.npz` and `inputs/T2_independent_geometry.npz` are compact
extracts from the completed T2 operator study. `EXTRACTION_PROVENANCE.json`
records the full parent matrix-file hashes and selected keys. The parent
provenance and validation report are included. This replay recomputes the
approximations from frozen geometry; it does not rebuild upstream Hodge
matrices. The unchanged MTFT 0.26.2 source archive is bundled for provenance
and its hash is checked, but no installed MTFT version is used.

`selection_frozen.json` and `selected_models.npz` preserve every training
subset score, basis, chosen rank, and training error. `approximation_results.json`
summarizes all 18 prescribed family/rank models. `approximation_pointwise.npz`
retains spectral and time errors, raw norm-transfer values, and conditional
absolute bounds. All ranks below five fail the joint budgets; all r=5 entries
are full-model controls.

The independent files reconstruct fresh sector bases from the full real26
operator and localization projector, score response by Schur solves, and
independently freeze their own selections before holdout. Their NPZ includes
the corresponding error arrays and ambient retained-subspace projectors.
The independent routes share upstream numerical geometry.

The 44 exact symbolic controls use synthetic matrices and do not certify
T2's floating-point coefficients. `RUN_NOTES.md` preserves the pre-holdout
serialization fix and the protocol's appended notation clarification.

No physical clock, Hamiltonian, mass assignment, global approximation bound,
or measured runtime speedup is established. The short-prefix diagnostic is
not a substitute for the failed full-domain approximation target.
