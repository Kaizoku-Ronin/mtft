# MTFT Standard Model compatibility assessment

Read `MTFT_Standard_Model_Assessment.md` first. The assessment is dated 7 September 2026 and applies to the attached MTFT 0.26.0 release. It records the state of the implementation and the mathematical/physical bridge; it does not modify the package or claim a complete audit of all papers.

To reproduce the numerical audit, make the attached package importable, then run:

```bash
python audit_sm_readiness.py
```

The recorded environment is Python 3.12, NumPy 2.3.5, and MTFT 0.26.0. The audit also uses MTFT's installed dependencies. It checks the mass pipeline, temporary changes to the VEV inputs, absence of a PDG lookup on that pipeline, a catalog mass/calibration round trip, the Higgs quartic relation, two MW/MZ formulas, group membership for two SU(3) samplers, and the 19 frozen surface gates. Temporary overrides exist only inside the running process and are restored.

`audit_results.json` retains source hashes and all numerical results. `package_honest_report.txt` is the package's existing prediction report, with its own stored reference values and tolerances; it is not a new fit against the 2026 PDG release. The evidence status of each conclusion is identified in the assessment. This is an appended audit record, preserving earlier results and source files.
