# MTFT R2C-01 reproducible study

Exact six-dimensional chirality investigation following the v0.32.0 release.

Read MTFT_R2C01_Report.md for the derivation and HANDOFF.md for integration.
The baseline is the supplied archive, identified in PROVENANCE.json.

## Run

From this directory, using Python 3.10 or newer:

~~~bash
python -m pip install -r requirements.txt
python chirality_audit.py --output results
~~~

Verified using Python 3.12.14 and SymPy 1.14.0. No MTFT installation is required.
The frozen Python input files are provenance snapshots; the script only executes
the standard-library-only parent record and parses the stack constants with AST.
It reconstructs anomalies independently rather than calling MTFT's anomaly code.

## Expected outcome

| Quantity | Exact result |
|---|---:|
| Sign assignments | 1,024 |
| Preserve six named family indices | 16 |
| Preserve complete M1 sector-index ledger | 4 |
| Permit all four elementary-scalar bilinears | 64 |
| Preserve families and permit those bilinears | 0 |
| Exact checks | 293 passed, 0 failed |

The changed-flux scan contains four scalar-compatible family branches and all
64 completions back to the ten-sector field list. The anomaly restrictions
assume the operator class and tensor-charge quantization stated in the report.
They are not no-go theorems for arbitrary extra matter, new fields or inflow.

The surviving vector/form route is an allowed interaction structure, not a
completed action, vacuum, gravity model or mass prediction.

## Integrity

SHA256SUMS.txt covers the bundled deliverable files except itself. To check it:

~~~bash
sha256sum -c SHA256SUMS.txt
~~~

Results can be regenerated in a separate directory to preserve the distributed
JSON/CSV files:

~~~bash
python chirality_audit.py --output reproduced_results
~~~

CSV booleans are written as True/False. Sector order is fixed and recorded.
Every anomaly coefficient is a rational SymPy expression; matrices use strings
such as I for the imaginary unit. No stochastic or floating-point calculations
are involved.

The upstream source snapshots retain their original authorship and serve only
to identify and validate the supplied research inputs.
