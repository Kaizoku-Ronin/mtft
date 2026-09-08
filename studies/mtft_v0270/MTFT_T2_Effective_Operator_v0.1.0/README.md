# Replay the T2 effective-operator experiment

Read `MTFT_T2_Effective_Operator_Report.md` for the findings and
`GENERAL_FEEDBACK_IDENTITIES.md` for the conditional proof.
`PROTOCOL.md` is the locally frozen measurement specification.

With Python 3.12 and the listed scientific dependencies:

```bash
python -m pip install -r requirements.txt
python run_study.py
```

The runner extracts the bundled supplied MTFT 0.26.2 archive into a temporary
directory and verifies all archived source files. It then runs the primary
complex block calculation, independent real-frame calculation, exact symbolic
controls, figure rendering, and cross-route comparisons. It does not use an
installed MTFT version or modify the release source.

To use an existing matching extracted source tree:

```bash
python run_study.py --source /absolute/path/to/mtft-0.26.2
```

Outputs overwrite the recorded result files beside the scripts. Run a copy of
this directory to preserve the original results and manifest. Timing, local
paths, image metadata, and last-place floating-point values can vary on replay;
the numerical comparison tolerance is 1e-9.

`inputs/previous` contains two frozen arrays and provenance from the preceding
MTFT_Active_Module_v0.1.0 study. The direct and localization projectors share
upstream Hodge-frame geometry; these are independent computation routes, not
independent physical observations. This replay consumes those frozen arrays
and does not rebuild their earlier geometric construction. The retained
`ORIGINAL_STUDY_SHA256SUMS.txt` describes the entire earlier study, not only
this subset. `inputs/W143_REFERENCE.json` records the prior baseline used in
the comparison table; the W143 study is not rerun here.

`t2_effective_operator_results.json` retains all 1,807 regular points, 55 pole
stress points, five real-pole null checks, 257 dynamics points, six short-time
measurements, and six matrix-exponential controls. Its companion NPZ stores
the complex matrices, residues, dark basis, and projected response matrices.
`independent_T2_results.json` contains the exact arithmetic characteristic
polynomial and independent checks; its NPZ stores the real26 construction and
all ambient complex projected responses. `independent_comparison.npz` stores
the cross-route response errors. `general_feedback_identity_results.json`
records 40 exact symbolic checks on synthetic Hermitian controls.

`VALIDATION.json` reports replay consistency separately from research outcomes.
The independent-pair decomposition fails for T2, initially dark directions gain
delayed coupling, and direct Schur evaluation loses accuracy in the pole stress
test. These findings remain in the bundle rather than being treated as passing
structural or stress-accuracy hypotheses. No physical Hamiltonian, clock,
energy scale, irreversible limit, or Standard Model correspondence is inferred.
