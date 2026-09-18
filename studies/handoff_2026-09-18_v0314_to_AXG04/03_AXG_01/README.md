# AXG-01 reproduction bundle

Read `MTFT_AXG01_Report.md` for the derivations, physical interpretation and limits.

This is an additive audit of the supplied MTFT 0.31.4 candidate, not a modified
MTFT release. It uses the declared representative M1 charge/flux data to test
a proposed two-axion effective extension.

## Run

Use Python 3.12 with the packages listed in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
python run_all.py
```

The run itself requires no network access or installed MTFT package. The exact
runtime used for the delivered outputs is recorded in `results_summary.json`.
Expected total: **159 exact assertions passed, zero failed**. The full MTFT test
suite is not part of this targeted audit.

## Files

| File | Role |
|---|---|
| `anomaly_audit.py` | Full polynomial factorization, descent, lattice, normalized vector masses; compares with copied source |
| `operator_audit.py` | Higgs/Yukawa/matter charges and integer axion dressings |
| `vacuum_audit.py` | Minimal kinetic extension, radius monotonicity, conditional D-term flavor degeneracy, sign-changing control |
| `run_all.py` | Runs all three and records the combined status |
| `*_results.json` | Exact subaudit outputs and individual checks |
| `results_summary.json` | Check counts, runtime and copied-source hash |
| `input/smflux_v0314.py` | Unmodified copy of `src/mtft/surface/smflux.py` from the supplied candidate |
| `input/MTFT_LICENSE.txt` | Original package's MIT license and copyright notice |
| `SCOPE.md` | Scope and conventions recorded for this investigation |
| `manifest.json` | SHA-256 hashes of packaged files, excluding the manifest itself |

The source snapshot retains its upstream license. Mathematical results are
conditional on the stated field content, axion charge map and background
assumptions. Passing assertions establishes the displayed algebra, not a
physical vacuum, experimental agreement, global anomaly cancellation or a
six-dimensional UV completion. No parameters are fitted.

## Reading order for hand replication

1. Report §2: multiply the charge matrix by the surviving generators.
2. §3: expand and factor the anomaly polynomial.
3. §§4–5: distinguish anomaly cancellation from kinetic normalization.
4. §6: test allowed operators by adding charge vectors.
5. §7: differentiate the restricted radius potential and examine the Hessian.
