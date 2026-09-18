# Integrity and reproduction

Python 3.12 was used for the delivered study records. The integrity verifier and runner use only the standard library. The mathematical studies need their original dependencies, recorded in their `requirements.txt` files and JSON metadata.

## Verify the handoff without running mathematics

From the extracted handoff root:

```bash
python tools/verify_bundle.py
python tools/reproduce.py --list
```

`SHA256SUMS` covers every delivered payload file except itself. `provenance.json` also records the original source hashes and the mapping from each original ZIP member to its extracted copy. The verifier checks the payload hashes, archive CRCs and exact archive-member correspondence. This is integrity verification, not a proof of the mathematics.

The delivered studies, original ZIPs and supplied baseline are unchanged. New handoff documents, tools and their separate packaging-validation record are identified by their locations; no baseline source was patched.

## Rerun one study in a separate working directory

Use a Python environment with the required numerical libraries. For AXG-04, for example:

```bash
python -m venv ../mtft-handoff-venv
../mtft-handoff-venv/bin/python -m pip install -r studies/06_AXG_04/requirements.txt
../mtft-handoff-venv/bin/python tools/reproduce.py --study AXG-04 --workdir ../mtft-rerun-axg04
```

On Windows, use the environment's `Scripts/python.exe` path instead. The runner uses its own Python interpreter for child scripts. It performs no dependency installation or network request. It copies the selected original study to a new working folder so its result-writing scripts do not alter the archive's frozen records. An existing study output folder is not overwritten.

For the initial audit and SC7-01, the runner also extracts the supplied v0.31.4 archive into the separate working folder and passes the source path explicitly. It does not use an unrelated installed MTFT version. It validates archive paths and rejects links/special entries before extraction.

```bash
python tools/reproduce.py --study V0314-AUDIT --workdir ../mtft-rerun-audit
python tools/reproduce.py --study SC7-01 --workdir ../mtft-rerun-sc7
python tools/reproduce.py --study STUDY-EDITION --workdir ../mtft-rerun-hand-math
python tools/reproduce.py --study HOPF-02 --workdir ../mtft-rerun-hopf
python tools/reproduce.py --study AXG-01 --workdir ../mtft-rerun-axg01
python tools/reproduce.py --study AXG-02 --workdir ../mtft-rerun-axg02
python tools/reproduce.py --study AXG-03 --workdir ../mtft-rerun-axg03
```

`--study all` runs all eight in catalog order once their dependencies are installed. HOPF-02's exact Manin/centralizer reconstruction and SC7's exact local-system scan can be slower than the small AXG polynomials. Each study has a log; `rerun_status.json` records process status and elapsed time. A zero process return code does not change an archived physical gate from FAIL to PASS.

The `STUDY-EDITION` choice runs the 187 hand calculations only. It does not rebuild or rerender the Word document. The original [paper reproduction instructions](teaching/paper_work/README_Reproduction.md) explain that separate workflow and its Pandoc/font/layout dependencies.

## Expected historical records

| Study | Expected recorded outcome |
|---|---|
| V0314-AUDIT | Symbolic/diagnostic program PASS; seven selected shipped tests are documented separately |
| SC7-01 | Exact computation PASS; required-spectrum physical gate FAIL; restricted radius gate FAIL |
| STUDY-EDITION | 187 exact hand checks |
| HOPF-02 | 85 passed, zero failed |
| AXG-01 | 159 passed, zero failed |
| AXG-02 | 84 passed, zero failed |
| AXG-03 | 149 passed, zero failed |
| AXG-04 | 322 passed, zero failed |

Runtime versions, source paths and environment metadata can differ in a rerun. Compare the mathematical outputs under their declared conventions rather than demanding byte equality of regenerated environment records. Archived outputs retain their original hashes.

No top-level command here reruns the claimed 810-test release suite. For the selected original package tests, follow `studies/00_v0314_audit/TEST_RESULTS.md`. A future integrated package release needs its own fresh-install and full-suite verification.
