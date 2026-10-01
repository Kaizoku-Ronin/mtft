# Reconciliation — v0.34.0 release-review corrections into the Matter-03 candidate (Claude, 2026-10-01)

Base: `current/repo_candidate/` of `MTFT_v0340_to_Matter03_Claude_Handoff.zip` (Matter-01/02/03 on v0.34.0; bundle integrity: 1130 payload
files match their SHA-256 manifest). Version unchanged (0.34.0); nothing pushed or published.

## 1. The three review patches
All three apply cleanly to the candidate. Every file they touch is byte-identical between pristine v0.34.0 and the candidate except
`docs/changelog/UNRELEASED.md`, whose hunk does not overlap the Matter sections. Applied in the order parent/geometry → S1 → CC-35.

Each finding was checked before acceptance; four of them correct code or wording written in the PV-01/PV-02/CC-35 handoff (Claude) and are
accepted as such:
- **CC-35 stiffness wording** (my error): `mu_N` is not the Hessian of the implemented single-trace Polyakov action. Independent recomputation:
  at the SU(3) centre-symmetric point (kappa = 1, y = 0.18174) the tangent gradient norm is 0.01928985 (not stationary) and the actual Hessian
  eigenvalues are 3.60067e-8 and 0.0984541, against 0.1476812 for both arithmetic combinations; at y = 0.3 the actual Hessian has the eigenvalue
  -0.00242941 while mu_3 = 0.05782 > 0. Identical to the review's numbers. The SU(2), n = 2 exception to the CC-35 lemma is also correct.
- **Anomaly gate** (my error): `u1(q)` content was silently dropped and returned a passing verdict; now OPEN with the content named.
- **B1 scope** (my error): `a_c*b_c = 0` was imposed, not forced; the review's U(8) assignment with `a_c*b_c = 1` gives exactly `SM_NET(3)` + 63
  neutral singlet indices (now a regression test). Box-3/4 counts and the finite trilemma stand for the imposed subclass only.
- **Blow-up API** (my errors): equal-stack and identically-zero slope loci; `collinear_family(a_L)` D-flat at eps/A_E = 1/3 and ample for r > 1/3
  for every nonzero a_L (my docstring scaled both with a_L); nine doublet pairs at a_L = 1; weak-triplet labels kept in `sm_net_content`.
- **S1 gluing ranks** (pre-existing v0.33.1 code): the count is 32 - 15*pairing_rank and independence does not fix the rank; the historical values
  are the generic ones.

## 2. Release navigation and the missing artifacts
- `CHANGELOG_v0340.md` created from the notes that shipped as `UNRELEASED.md` (verbatim under a release title); `UNRELEASED.md` now holds only
  post-v0.34.0 work (review corrections, Matter-01/02/03); changelog index and README "current line" updated.
- Root cause of the three artifacts missing from the v0.34.0 sdist: `MANIFEST.in` shipped only `docs/**/*.md`, two named `viz/` files, and no
  `scripts/cosmology/`. Rules added: `docs *.png`, `viz/*.html viz/*.png`, `scripts/cosmology *.py`. Restored from the PV/CC-35 handoff:
  `docs/SM/figs/cc35_polyakov_kernel.png` (sha256 35b2c044…) and `viz/arithmetic_susy_picture.html` (sha256 629ce143…).
  `scripts/cosmology/cc34_cosmology_audit.py` exists in no sdist since it was introduced; it must be taken from the repository.
- **Latent publish-gate failure in the Matter-03 candidate:** `tests/test_boundary_and_gprun.py::test_manifest_covers_every_studies_extension`
  failed because the Matter studies add `.html`, `.png` and `.cjs` files that the studies rule did not ship. Rule extended; the test passes;
  a rebuilt sdist (7.2 MB) contains the restored figures, all labs and the 21 Matter study visual files.

## 3. Gates on the reconciled tree
- Affected regressions (review tests, PV-01/02, CC-35, parent ledger, Matter-01/02/03, docs, Legend, release pin): 246 passed.
- Complete routine suite (`-m "not slow"`), before the studies-manifest fix: 1120 passed, 1 failed (the manifest test above), 4 skipped,
  40 deselected. After the fix, full rerun: **1121 passed, 4 skipped, 40 deselected** (6 min 50 s, Python 3.12, no PARI/GP).
- Not run: the 40 slow cases (the PV-01 box-3/4 scans passed in the release review), PARI/GP-dependent tests, browser layout of the labs.

## 4. Not reconciled here
- The panel-6 coherent-modes extension of `viz/arithmetic_susy_picture.html` (`MTFT_Visualization_Extension_v0340.zip`, from the same R&D
  sessions) is not in the handoff bundle; the original page is restored, the extension awaits its zip.
- TOT/ADS studies remain standalone, as the handoff specifies.
