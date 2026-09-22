# MTFT graded geometry investigation

Read `INVESTIGATION.md` for conclusions, derivations, source links, and limits.

Run from this directory using Python 3.10+ with NumPy:

```bash
python verify_graded.py
python audit_operators.py
python audit_bracket.py
```

The operator audit alone uses only the Python standard library. No MTFT installation or network connection is required. Results are written to `graded_results.json`, `operator_results.json`, and `bracket_results.json`. Input fixtures are hash-pinned and retain the upstream MIT license.

The bundle investigates mathematical structures in MTFT 0.32.0. It does not modify or publish the package, and does not establish a physical Monster symmetry, Feigenbaum renormalization law, coupling constant, or gravitational model.
