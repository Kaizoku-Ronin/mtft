# AXG-02 reproduction bundle

Read **MTFT_AXG02_Report.md** for the findings, conventions, proofs and scope.
This is an additive investigation using the supplied MTFT 0.31.4 M1 inputs.

## Reproduce

With Python 3.12:

    python -m pip install -r requirements.txt
    python run_all.py

The run itself is offline and does not require installing MTFT. It writes two
exact result files and results_summary.json, which records versions and counts.
The delivered run has 84 passing exact assertions and no failures.

## Contents

| File | Purpose |
|---|---|
| parent_anomaly.py | 6D anomaly via characters and independent Cartan weights; exact flux pushforward; conditional chirality repair |
| flux_and_scalar.py | Smooth-scalar flux obstruction, integer axion basis, gauge representation check, conditional potential classification |
| run_all.py | Runs both audits and records combined status |
| *_results.json | Exact expressions and individual checks |
| results_summary.json | Counts, dependency versions and input hash |
| input/smflux_v0314.py | Unmodified input module from the supplied archive |
| input/MTFT_LICENSE.txt | Original module copyright and MIT license |
| input/SM01_PREREGISTRATION.md | Original declared parent and counting assumptions |
| input/REVIEW_V0311_RESPONSE.md | Original product-group clarification |
| SCOPE.md | Explicit assumptions and limits of this investigation |
| manifest.json | SHA-256 hashes of packaged files except the manifest itself |

The target is necessary consistency and vacuum conditions. No new parent,
supersymmetric completion, physical mass scale, or experimental agreement is
claimed. The arbitrary-unit negative-potential example is a mathematical control.
