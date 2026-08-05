# mtft v0.11.2 — packaging fix: sdist self-containment (2026-08-05)

Scope: **packaging only.** No code changes; the wheel's contents and the
package's import surface are unchanged from v0.11.1.

## Fixed — source distribution could not reproduce its own release gate

Found by external review of the PyPI tarball (ChatGPT, on the v0.11.1
sdist) and independently reproduced by the auditor (Kimi K3) before any
change was made:

- `mtft-0.11.1.tar.gz` shipped `tests/test_x0143_particle_box.py` but
  **not** the three study engines it imports
  (`studies/x0143_particle_box{,_v02,_v03}.py` — the tessellation,
  Hecke/period, and dissipative-capture engines of the particle box).
- Running the suite from the published sdist therefore gave
  FileNotFoundError on the studies path. Targeted reproduction
  (auditor, on the actual PyPI artifact):
  `1 failed, 1 passed, 6 errors` in `tests/test_x0143_particle_box.py` —
  every non-pass the same missing-file cause, none a numerical failure.
- Root cause: the repo had no `MANIFEST.in`; setuptools' default file
  list picked up `tests/` but not `studies/`.

This mattered because the project presents the particle-box tier as part
of the release gate — the published source artifact must be able to
reproduce that gate on its own.

## Changes

- **`MANIFEST.in` added**: `studies/` (engines + README), `tests/`,
  `CHANGELOG_*.md`, `CITATION.cff`, and the publishing docs are now
  explicit members of the sdist. (The wheel still ships only the
  `mtft` package — unchanged and correct.)
- **Publish gate hardened** (`.github/workflows/publish.yml`): after
  building, CI now unpacks the just-built sdist and runs the particle-box
  test tier *from the unpacked tarball* — the exact reproduction path
  that failed for v0.11.1. A future packaging regression of this class
  fails the gate before anything reaches PyPI.

## Verification (auditor, pre-push)

- Rebuilt sdist inspected: `studies/` present (all engines), `tests/`
  complete, changelogs and docs included.
- Full suite run from the unpacked, freshly built sdist:
  454 passed, 1 skipped, 0 failed — the tarball now reproduces the gate.

Audit lineage: external review (ChatGPT) → reproduction + fix (Kimi K3)
→ Roger.
