# Unreleased — after v0.34.0

## v0.34.0 release-review corrections (reconciled 2026-10-01)

- Search scope: `anomaly8d.b1_scan` imposes `a_c*b_c = 0` and two individual mixed antitriplet blocks of index -3. Its 60/156
  assignments and the finite trilemma are statements about that imposed subclass, not an exhaustive B1 classification; an explicit
  U(8) assignment with `a_c*b_c = 1` has exactly `SM_NET(3)` plus 63 neutral singlet indices (`tests/test_parent_edge_cases_review.py`).
- `parent.gates`: unsupported content (e.g. `u1(q)`) makes the 8D anomaly verdict OPEN with the content named, instead of being dropped;
  identically zero slope conditions are reported as the whole positive-r locus.
- `surface.blowup`: equal-stack slopes give the whole ample cone in `polystable_locus`; `collinear_family(a_L)` has its D-flat locus at
  eps/A_E = 1/3 and ampleness r > 1/3 for every nonzero a_L (the released docstring's a_L/3 scaling was wrong), nine doublet pairs at a_L = 1.
- `anomaly8d.sm_net_content` keeps weak-triplet labels.
- `product_surface.s1_colour_recombination_dims`: gluing dimension 32 - 15*pairing_rank; independence alone does not fix the rank;
  the historical 17, 2, 2 values are the generic counts (`tests/test_s1_recombination_rank.py`).
- CC-35 wording: mu_N(y) is a defined arithmetic stiffness diagnostic, NOT the Hessian of the implemented single-trace Polyakov action
  (at the SU(3) centre-symmetric point, kappa = 1, y = 0.18174, the tangent gradient is 0.01928985 and the actual Hessian eigenvalues are
  3.60e-8 and 0.0984541; at y = 0.3 one is -0.00242941). The CC-35 lemma has one exception: SU(2), n = 2 (Tr U^2 = chi_1 - 1).
- Release navigation: `CHANGELOG_v0340.md` added (the shipped notes, verbatim); README points at the v0.34.x line; MANIFEST.in now ships
  `viz/*.html`, `viz/*.png`, `docs/**/*.png` and `scripts/cosmology/*.py` — the cause of the three artifacts missing from the v0.34.0
  sdist (they were excluded by the manifest). `docs/SM/figs/cc35_polyakov_kernel.png` and `viz/arithmetic_susy_picture.html` restored
  from the PV/CC-35 handoff; `scripts/cosmology/cc34_cosmology_audit.py` must come from the repository (not in any sdist).
- Record: `docs/SM/RECONCILIATION_V0340_REVIEW_MATTER03.md`.

## Matter-01 (2026-10-01): U(1) reference engine and conditional bridge

- New `mtft.matter` CLI and massive explicit-spinor tree-level QED amplitude
  for electron-positron annihilation into two photons; full-sphere event
  factor, Ward/Bose/gauge-shift/boost checks and independent scalar oracle.
- Frozen PDG reference mass/coupling inputs with explicit override provenance.
  These are not MTFT predictions; no quantum-state time evolution is added.
- Conditional same-charge Gram normalization, electroweak matching and uniform
  Einstein-frame breathing gauge prefactor; C3X and B3 remain distinct.
- Exact B3 flux-class audit and RR cohomology lower bounds. The archived
  minimum-index-extras witness can have much larger cohomology bounds than
  the simpler collinear control. Actual spectra, normalized Yukawas, exotic
  masses and the photon axion/Stueckelberg kernel remain OPEN.
- Offline interactive lab, scientific overview, atlas evidence ledger without
  mapping promotions, focused runner and [report](../SM/MATTER01_REPORT.md).
- Existing results and release version are preserved. Full release tests and
  long archived studies are separate from this milestone gate.

### Matter-02 collinear representative and holomorphic coupling

- Specify the cusp-at-infinity blow-up and trivial elliptic flat factor.
- Add exact cohomology dimensions, rational cusp-jet certificates and a
  globally valid rank-three lepton–Higgs holomorphic residue witness.
- Append the opposite archived/BHV chirality conventions without changing
  released index functions or source data.
- Audit hypercharge traces and explicitly INPUT axion/Higgs mass kernels;
  physical photon and normalized particle masses remain open.
- Add the local interactive mode lab, scientific overview, independent
  certificates and a focused reproducibility runner.

### Matter-03 Higgs branch, deformed light modes and kinetic inputs

- Derive the arbitrary-Hermitian-metric obstruction to a nonzero single
  directed Higgs on the unchanged full B3 wall background.
- Certify an explicit reverse field and odd exceptional spectral valuation
  seven; canonical-twisted surface stability supplies classical paired
  F/D-solution existence. No physical metric or Higgs scale is computed.
- Add a charged deformation-complex census: six aa ranks total 41 and leave
  237 degree-one classes per sign for the fixed pair; all aligned pairs in
  these spaces leave at least 160. Observed family selection remains open.
- Add lower-index field-coordinate kinetic normalization with explicit Grams,
  canonical Higgs direction, vev and parent factor, including complex-basis
  covariance tests. Historical section normalization remains unchanged.
- Add derivations, source APIs, focused runner, offline lab and scientific
  overview. Standard Model metastability is kept distinct from this branch.
- Preserve v0.34.0 release pins, original arithmetic data and prior archives.
