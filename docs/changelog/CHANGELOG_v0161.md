# v0.16.1 — sdist completeness hotfix (2026-08-18)

One-line packaging bug, found by Claude's post-release artifact review,
auditor-verified against both published sdists.

## The bug

`MANIFEST.in` read `recursive-include studies *.py *.md` — every `.txt`,
`.json`, and `.log` under `studies/` was silently dropped from the
sdist, and `recursive-include tests *.py` dropped
`tests/zeros_gamma_T100.npy`.

**Impact on v0.16.0:** the sdist shipped `studies/ci_2026aug/` with 13 of
23 files — all prose and `ci_verify_kimi.py` survived, but the 7 wave
data files, the certificate JSON, `ci_verify_kimi.json`, and
`sha1_manifest.txt` did not, so the byte-preservation claim was
unverifiable from the sdist and the auditor script raised
`FileNotFoundError` there. **The git repository was always complete**
(post-push byte-integrity 303/303 over the git tree); this was an
sdist-only omission.

**Older than this arc:** the v0.15.0 sdist has the same class of
omission — the v0150_2026aug bundle's auditor JSONs
(`m8_verify_kimi.json`, `c9b_*.json`), the V6 period data
(`period_matrix_manin_v6.txt`, `x0143_period_data_v6.json`, GP
transcript), `run_v6.log`, and the ledger JSONs were all dropped at
v0.15.0. And the missing `tests/zeros_gamma_T100.npy` explains the
v0.15.0 gate's 611-passed/2-skipped sdist run versus 612/1 in-repo: the
weil E2 test silently skipped in the sdist. Now diagnosed and fixed.

## The fix

    recursive-include studies *.py *.md *.json *.txt *.log
    recursive-include tests *.py *.npy

Verified by building the sdist locally and diffing its `studies/` and
`tests/` listings against the git tree: complete.

## Note on sector-ordering conventions

Two orderings of the AL joint sectors are in circulation, both correct:
the CI reports use (+,+), (+,−), (−,+), (−,−) → (1, 6, 5, 1); the
v0.16.0 changelog and the bundle README use (+,+), (−,+), (+,−), (−,−)
→ (1, 5, 6, 1). Same data. The convention will be pinned explicitly in
`mtft.canonical` (v0.17.0 wave); a note is also added to the ci_2026aug
bundle README.

## Next

`mtft.canonical` (v0.17.0): B, the adapted q-expansions, the four class
matrices and the descent table as frozen certified data under
`src/mtft/canonical/_data/` (the `crypto/_data` pattern — sidesteps this
bug class), with P1–P9, the ten bundle-rank tests and the nine-row
projection table as test gates. Registered and still unrun: the disc²
cross-level falsification (AG-D5), and the CC-09 downstream sweep for
claims keyed to h(−143) = 7 or degree-7 CM minimal polynomials.
