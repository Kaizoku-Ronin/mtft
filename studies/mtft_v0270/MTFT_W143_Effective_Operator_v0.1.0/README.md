# Replay the W143 effective-operator study

Read `MTFT_W143_Effective_Operator_Report.md` for the assessment and
`SCHUR_FEEDBACK_IDENTITIES.md` for the conditional algebraic proof.
`PROTOCOL.md` was frozen before measuring the new canonical modes and responses.

With Python 3.12 and the listed scientific dependencies:

```bash
python -m pip install -r requirements.txt
python run_study.py
```

The runner extracts the bundled, unchanged MTFT 0.26.2 archive into a temporary
directory, recomputes the primary and independent routes, runs the symbolic
control, renders the figure, and writes `VALIDATION.json` and `run_log.txt`.
No installed MTFT version is used. Alternatively, provide a matching extracted
source tree with `python run_study.py --source /absolute/path/to/mtft-0.26.2`.
The runner verifies all archived source files before using that tree.

Outputs are written beside the scripts and overwrite the supplied result files.
Run a copy of this directory to preserve the original results and manifest.
Floating-point values, timing, PNG metadata, and absolute paths may vary across
runs; numerical comparisons use the documented tolerance rather than file hashes.

The two numerical routes share frozen Hodge-frame data from the preceding
`MTFT_Active_Module_v0.1.0` study. Its two input arrays and provenance are bundled
under `inputs/previous`. The earlier full-study manifest is retained there as
`ORIGINAL_STUDY_SHA256SUMS.txt`; it describes that earlier complete study, not
just the input subset included here. This runner reuses those frozen inputs and
does not rebuild their upstream construction.

`effective_operator_results.json` contains all 907 regular points, all 22
conditioning stress points, real compression-pole evaluations, and six direct
matrix-exponential checks. `effective_operator_matrices.npz` retains the full
and block matrices. `canonical_check_results.json` and its NPZ preserve the
independent positive-eigenspace route. `feedback_identity_results.json` records
27 exact symbolic checks on a synthetic generic channel.

The registered response gate applies to the regular grid. The raw Schur
calculation loses accuracy in the conditioning stress test; this is preserved
and reported, not turned into a passing regular-grid result. No physical
Hamiltonian, energy scale, time unit, or QFT identification is inferred.
