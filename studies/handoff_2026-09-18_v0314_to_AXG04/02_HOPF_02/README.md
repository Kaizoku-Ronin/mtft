# HOPF-02 reproduction bundle

Read **MTFT_HOPF02_Report.md** first. It includes hand derivations, group
matrices, operator conventions, the interpretation of S^143, and limitations.

To reproduce the exact computations:

```bash
python -m pip install sympy
python investigate.py
```

Python 3.10 or newer is required. No MTFT installation or external data is
needed. The script writes `results.json` alongside itself and stops on a
failed check. Expect the rational Manin-symbol reconstruction and centralizer
calculation to take some time before the next progress message.

- `PREREGISTRATION.md`: declared calculation scope, written before the run.
- `investigate.py`: exact symbolic verification code.
- `results.json`: full matrices, transformations, permutations, and check results.
- `input/hecke_v0314.py`: unchanged computational source from the supplied archive.
- `input/MTFT_LICENSE.txt`: license for that source.
- `manifest.json`: SHA-256 hashes of deliverables and source provenance.

The study is separate from the MTFT release. The positive-control trace metric
is constructed here and is not asserted to equal the Petersson metric. The
new topological arguments and the inherited spin-lift inputs are distinguished
in the report.
