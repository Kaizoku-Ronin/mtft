# Matter-03 reproducibility guide

This focused study establishes a classical paired Higgs branch in the fixed
B3 bulk parent, recomputes its charged deformation complex, and implements
normalization with explicit physical inputs. It leaves 237 degree-one
charged classes per sign, has no selected Higgs scale, and has no computed
B3 kinetic metrics or observed mass predictions.

From a checkout containing these unreleased additions:

```bash
python -m pip install -e ".[dev,viz]"
python studies/matter03/run_all.py
```

The runner prints verdicts and writes each log and JSON result under
`results/`. It does not execute the full release suite or hours-long studies.
Node.js is optional; without it the lab logic checks are explicitly SKIPPED.
`--skip-tests` also records the requested skip. The HTML works offline.

| File or directory | Purpose |
|---|---|
| `vacuum/` | 57 exact arithmetic checks; frozen canonical vectors and stability/existence proof |
| `modes/` | 34 exact finite checks; six map ranks and the Leray proof controlling extension terms |
| `kinetics/` | 117 checks: 11 exact symbolic, 106 numerical or invalid-input examples |
| `viz/` | Offline lab and Node/mock-DOM checks; no actual browser layout certification |
| `plot_matter03.py` | Rebuild the standalone scientific overview |
| `results/run_summary.json` | Actual stage statuses, timings, interpreter and dependency versions |

The public source regressions include the new kinetic, vacuum and charged
APIs and the relevant Matter-01/02, atlas and B3 tests. Mathematical existence
also invokes the cited surface theorem; a finite arithmetic script does not
prove that theorem or solve the D-term PDE.

Read `docs/SM/MATTER03_REPORT.md` first. Follow the linked derivations and
`MATTER03_PHYSICAL_INPUTS.md` for the exactness boundary, missing metric data,
and the distinct Standard Model metastability question. Historical archives
are kept intact; this source still carries the unreleased base version 0.34.0.
