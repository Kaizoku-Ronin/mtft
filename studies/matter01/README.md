# Matter-01 reproduction

Run from the repository root:

```bash
python -m pip install -r studies/matter01/requirements.txt
python studies/matter01/run_all.py
```

Read [the report](docs/MATTER01_REPORT.md), [QED derivation](docs/QED_DERIVATION.md),
and [vacuum audit](docs/VACUUM_AUDIT.md). The runner retains all JSON and logs
in `results/` and rebuilds the standalone lab and scientific overview in `viz/`.
Node.js is optional for the JavaScript runtime gate. This focused gate does
not start MTFT's full release suite or long archived studies.

Electron mass and alpha are supplied inputs. The physical MTFT electron,
photon kernel, full mode spectrum and normalized Yukawas remain OPEN.
