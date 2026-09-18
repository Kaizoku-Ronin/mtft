# Reproducing the MTFT v0.31.4 audit

This bundle adds an audit; it does not modify the supplied MTFT source.

Contents:

- `MTFT_v0314_Spin_and_Vacuum_Audit.md`: findings, derivations, limitations and source links.
- `audit_v0314.py`: exact eta multipliers, D8/Q8 group checks, canonical-square compatibility, AL fixed-point counts, symbolic anomalies and frozen tensor reconstruction.
- `results.json`: recorded results and source-data hashes.
- `run.log`: audit program output.
- `TEST_RESULTS.md`: targeted shipped-test results.
- `input/mtft-0.31.4.tar.gz`: the supplied input archive, unchanged.
- `manifest.json`: SHA-256 hashes of bundled files.

Create a Python environment using `requirements.txt`, extract the input archive to a directory, and run:

```bash
python audit_v0314.py --source /path/to/mtft-0.31.4 --output results.json
```

The script imports MTFT from that source directory. It needs no network access and does not invoke PARI/GP. It checks the eta transformations independently of the frozen AL matrices, then uses those matrices as a separate consistency check. Numerical normalization reconstruction is explicitly diagnostic. The seven targeted tests are not a rerun of the reported full release suite.

Use `source/src` on `PYTHONPATH` to repeat the shipped tests named in `TEST_RESULTS.md`.
