# Handoff studies: MTFT v0.31.4 candidate -> AXG-04 (Astra, 16–18 September 2026)

Frozen import of the research handoff `MTFT_v0314_to_AXG04_Claude_Handoff` (baseline: the v0.31.4 candidate,
SHA-256 2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5).  Original reports, scripts and result
ledgers are unchanged; `SHA256SUMS` and `provenance.json` are the handoff's own integrity records.

Excluded from this import: the two embedded copies of the baseline sdist (`studies/00_v0314_audit/input/` and
`studies/01_SC7_01/input/`, 3.7 MB each — the baseline is this package's own history) and the teaching edition
(held by the author for discretionary use; not part of the software).

Study scripts execute at module scope and write JSON beside themselves; they are regression references, not library
code.  Rerun them in a separate working directory (`reproduce.py --study <ID> --workdir <dir>`), never in place.
Historical check counts (HOPF-02 85, AXG-01 159, AXG-02 84, AXG-03 149, AXG-04 322) are the ledgers' own records;
AXG-04 was rerun during import (322 passed, 0 failed).

Reconciliation of the handoff's H-01…H-30 with the project correction register: `docs/SM/HANDOFF_RECONCILIATION.md`.
