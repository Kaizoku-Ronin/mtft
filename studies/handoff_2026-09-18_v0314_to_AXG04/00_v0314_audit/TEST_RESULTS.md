# Targeted shipped tests

Run date: 2026-09-16. Source: the supplied mtft-0.31.4 archive.

Python 3.12.14; NumPy 2.3.5; SciPy 1.17.0; SymPy 1.14.0; mpmath 1.3.0; pytest 9.1.1.

Two test commands were run with the extracted `src` directory on `PYTHONPATH`:

```text
python -m pytest -q tests/test_surface_review_v0311.py tests/test_surface_arithspin.py
4 passed in 6.21s

python -m pytest -q tests/test_surface_m2.py
3 passed in 0.59s
```

All seven tests passed. They cover the shipped section-space D8 test, the existing theta/cuspidal/unit gates, the anomaly/Majorana additions, and the corrected M2 normalization checks.

The independent audit adds the canonical-square compatibility condition that the shipped D8 test does not cover. Passing that original test is consistent with finding that its W11 lift does not preserve the specified theta-square map.

The reported full 810-test release suite was not repeated. No FEM mesh, quantum vacuum, supersymmetric parent or six-dimensional gravitational anomaly calculation was certified by this targeted run.
