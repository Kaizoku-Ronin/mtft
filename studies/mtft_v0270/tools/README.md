# Small integration tools

These standalone tools add no public package API and do not apply the historical v0.26.0 correction patch.

- `check_release_contracts.py` checks the release's explicit Higgs conventions, SU(N) group membership, reproduction of the historical determinant defect, and one finite local gauge-action covariance control.
- `replay_electron_photon.py` takes an explicit source tree and frozen study directory, copies the required scripts and inputs to a new output directory, and replays the finite experiments with the compatible convention checker. It refuses an existing output directory and any output inside the source or frozen study. Original reports and historical manifests are not copied into the new replay as if they described fresh results.

Python 3.12 with NumPy, SciPy, SymPy, and MTFT's dependencies was tested. Matplotlib is required only for `--render`.

```bash
OPENBLAS_NUM_THREADS=1 python replay_electron_photon.py \
  --study /path/to/MTFT_Electron_Photon_v0.1.0 \
  --source /path/to/mtft-0.26.2 \
  --output /path/to/new-electron-photon-replay --render
```

The finite experiments do not derive QED or identify the Standard Model from MTFT. They preserve the original study's supplied QED ingredients, fixed reference inputs, and numerical scope. A new source version is identified in the fresh provenance record and must pass the executable gates on its own.
