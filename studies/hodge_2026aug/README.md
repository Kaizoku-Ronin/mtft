# studies/hodge_2026aug — the du03 Hodge-session records (2026-08-06/07)

Measurement records of the du03 Hodge/metric session: the transported-well
no-go, the Grok-proposal triage, the commutant computation, and the CC-02
metric correction. Files are byte-identical to the session artifacts:

du03_dispersion (EARLY version — transported wells, ILL-POSED verdicts;
superseded by studies/du03_dispersion.py, kept here as the session record
and because du03_grok_triage / du03_metric import its `geometric_well`)
-> du03_grok_triage (three candidate escapes tested and closed)
-> du03_commutant (60/5/1 dims, Jordan structure, 1|12|12|1 blocks)
-> du03_metric + CC-07_du03_metric_correction (the two-metrics result:
Euclidean 5 vs polarization 43, joint admissible space = scalars only;
the du03 obstruction is metric-independent; anchor count stays 2).

Run note: these scripts import each other and du02_cycle_space_map by bare
name; run with the repo studies/ dir on PYTHONPATH. du03_metric's
eigenvector route to the invariant metric is LAPACK-fragile with the
doubled-line degeneracies — the K3 setup report (2026-08-07) certifies all
decisive claims via an SVD/Dykstra route that needs no eigendecomposition.
