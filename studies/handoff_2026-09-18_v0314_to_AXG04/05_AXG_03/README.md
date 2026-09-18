# MTFT AXG-03 reproduction

This additive investigation tests the explicit AXG-02 six-dimensional matter
assignment. It does not modify MTFT or construct a complete anomaly-free parent.
Read `MTFT_AXG03_Report.md` for the assumptions, hand calculations and conclusions.

## Run

Python 3.12 was used. From this directory:

```bash
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_all.py
```

The four scripts overwrite only their generated JSON results and the combined
summary. No network is needed after installing the dependencies. The scripts
use the bundled inputs; no sibling investigation directory is required.

| Script | Question |
| --- | --- |
| `tensor_factor_audit.py` | Does the residual 6D anomaly fit the proposed tensor pairing? |
| `tensor_flux_audit.py` | What integer axion charges follow from tensor flux reduction? |
| `higgs_parent_audit.py` | Can elementary scalar Higgs fields be massless and have the proposed 6D Yukawas? |
| `discrete_operator_audit.py` | What are the candidate discrete selection rules and fermion anomaly? |

The tensor audit reconstructs the complete chiral anomaly polynomial from
Chern characters and compares it to the archived AXG-02 polynomial. Its main
inertia result uses a restriction without the ambiguity of a general quartic
Gram presentation. The scalar bound is an analytic theorem; the script checks
its inputs and exact consequences, not a numerical approximation to a spectrum.

`input/smflux_v0314.py` is copied from the supplied MTFT 0.31.4 source distribution;
its original license is included. `input/axg02_parent_anomaly_results.json` is
the prior investigation's archived result. `manifest.json` records SHA-256
hashes of bundle files except itself. The manifest is a record of the delivered
version; running the scripts with a different Python patch version can change
the environment metadata in `results_summary.json`.

No physical mass in GeV, quantum-consistent discrete remnant, stabilized vacuum,
or dynamical gravitational completion is asserted by these checks.
