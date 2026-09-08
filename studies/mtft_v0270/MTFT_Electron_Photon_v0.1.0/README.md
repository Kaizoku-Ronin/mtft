# MTFT electron–photon investigation

Read `MTFT_Electron_Photon_Report.md` first. `Electron_Photon_Results.png` summarizes the finite computations. This study distinguishes the MTFT internal geometry from supplied QED structure and physical reference inputs.

The bundle includes a corrected local source archive, `mtft-0.26.0-sm-audit-20260907.tar.gz`, and the reviewable `SM_Audit_Corrections.patch`. This is an unpublished audit variant of 0.26.0, with two corrections and their retired historical accessors. It does not replace or rename the upstream release.

## Reproduce

Use Python 3.12 and the packages listed in `requirements.txt`, plus MTFT's declared dependencies. Extract the source archive and make its source directory importable:

```bash
tar -xzf mtft-0.26.0-sm-audit-20260907.tar.gz
export PYTHONPATH="$PWD/mtft-0.26.0-sm-audit-20260907/src${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
python run_all.py
python -m pytest -q mtft-0.26.0-sm-audit-20260907/tests/test_sm_audit_corrections.py mtft-0.26.0-sm-audit-20260907/tests/test_predictions.py::TestHiggsSector
```

Alternatively apply the unified patch with `patch -p1 < SM_Audit_Corrections.patch` from a clean 0.26.0 source directory. Do not apply it again to the already corrected archive. `inputs/original/` preserves the two original source files for comparison.

The exact projectors and frozen X0(143) input are included, so no GP, Sage, external database, or network is needed to rerun the experiments. Execution uses the provided frozen data, not a fresh modular-symbol reconstruction. Floating-point results may differ slightly by BLAS/LAPACK build; interpret them against their stated gates.

## Files

- `EXPERIMENT_PLAN.md` and `LOCAL_COVARIANCE_PLAN.md`: plans recorded before the respective experiments.
- `check_corrections.py`: public parameter consistency, SU(N) draws, and local SU(3) action invariance.
- `charge_study.py`: exact sector identities, 256 charge/mixing samples, and stripped tree Ward checks.
- `local_qed_covariance.py`: 12 finite Euclidean Wilson–Dirac covariance configurations.
- `photon_atomic.py`: sampled Fourier-symbol Maxwell tests and radial Coulomb convergence.
- `magnetic_and_coupling.py`: existing alpha alternatives and the Pauli freedom test.
- `render_results.py`: regenerate the results figure from the JSON files.
- `*_results.json`, `charge_exact_certificate.json`, and `charge_operators.npz`: data and reusable operators.
- `PROVENANCE.json`, `INDEPENDENT_REVIEW.md`, `targeted_tests.txt`, and `SHA256SUMS.txt`: source/runtime identification and validation record.

`SHA256SUMS.txt` describes the delivered files. Rerunning simulations or rendering can change file bytes; it does not update that delivery manifest. The report is the written interpretation of this recorded run, not an automatically regenerated claim about later edits.
