# MTFT 0.26.2 claim audit

Read `MTFT_v0262_Mistral_Assessment.md`. The report assesses the supplied Mistral review against the attached release and includes exact controls, finite computations, and explicit evidence limits. No release source was changed.

`CORRECTION_TO_PRIOR_HIGGS_AUDIT.md` corrects our previous wording about the two Higgs coupling conventions. `CLIFFORD_GATE_NEXT.md` proposes the next bounded operator-construction experiment; it is not a completed CG-01 result.

## Reproduce

The original uploaded source archive is included under the safe filename `inputs/mtft-0.26.2.tar.gz`. In a Python 3.12 environment with the package dependencies and the scientific libraries listed in `requirements.txt`:

```bash
tar -xzf inputs/mtft-0.26.2.tar.gz
python run_audit.py --source mtft-0.26.2
```

The runner configures the source import path and one OpenBLAS thread, executes the five bounded scripts, and runs three targeted release tests. No Sage, GP binary, network lookup, full genus-13 Pfaffian sum, or theta-function evaluation is required. The two-frame theta input census uses data already packaged with MTFT.

Individual scripts:

- `release_mellin_accounting.py`: convention conversion, SU(N) checks, finite Mellin normalization, prediction-row dependency accounting, and the trace-based ruler.
- `d4_check.py`: numerical D4 fingerprint and an exact standard triality control.
- `bimodule_control.py`: synthetic exact non-involutive example and TH2 input census.
- `finite_graph_check.py`: exact Ising DOS/graph controls and finite diffusion diagnostics.
- `render_finite_controls.py`: plot from the recorded exact DOS and numerical heat traces.

The JSON files contain source hashes, exact versus numerical classifications, and the full results. `SHA256SUMS.txt` is the manifest of the delivered snapshot. Reruns may change timing, runtime-specific paths, plot metadata, and floating-point bytes without changing the conclusions. They do not regenerate the written assessment or its manifest.

The original 0.26.2 source remains intact. The previous local 0.26.0 correction patch is not needed for this release.
