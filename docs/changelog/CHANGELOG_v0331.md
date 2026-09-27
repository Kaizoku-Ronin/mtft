# mtft v0.33.1 — SM interaction atlas, repository navigation, and the v0.33.0 review corrections (2026-09-26)

## From the release review of v0.33.0 (24 September 2026)

- Theta precision: `theta_torus.J_143A1` was an `mpf` rounded at module-import precision and both tau routes consumed it, so raising `dps`
  later could not recover the lost digits; agreement of the two routes was agreement on a rounded input. Now the rational j is rebuilt at each
  routine's working precision (`J_143A1_NUMERATOR/DENOMINATOR`, `_j_143a1_at_working_precision`); the historical constant is kept for callers.
  Regression: `tests/test_v0330.py::test_tau_precision_does_not_depend_on_import_history` (a 15-digit import followed by a 50-digit request).
  This repairs the input only; it is not an interval certificate and the series truncations are unchanged.
- S2 mass ratio: the sentence "m_c/m_t depends only on the four torus Higgs directions" (R2C-09, `theta_torus` docstring, compendium Part E
  App. R 55) contradicted the factorisation M = A (x) b(v)^T: the torus direction only scales M, and the ratio of the two nonzero singular
  values is that of the curve pairing A (undefined when b(v) = 0). Corrected in the docstring, `up_mass_rank_bound`, an R2C-09 addendum and
  `s2_ratio_independence`; test `test_s2_ratio_is_curve_only`. Identification with a physical charm/top ratio still requires the vectorlike
  states and the light spectrum.
- Notation (compendium Part E §0.12): the family type's curve factor is F_X = S0 (x) M with M = O_X(D) the flux bundle of degree a; the spin
  twist is not inserted twice. Part E's App. U row for V.3 is superseded by Lemma E.2 (the degree-six count is closed). Manuscript-only.
- Paper records: CC-33, tensor survivors, gaugino divisibility, point counts, gonality, the quotient certificate, the rank statements and both
  closed-form theta matrices refresh in about a second with the reviewer's `export_v0330_updates.py` (handoff, `01_release_review/`); the
  older GP/terminal kit stays pinned to 0.32.1 and must not label a rerun as 0.33.x unchanged.

## Exact block census of split surface backgrounds (27 September 2026)

- `product_surface.block_census(r, model)`, `census_scan`, `curve_h0_flux`, `curve_h1_flux`: V.1–V.3 applied leg by leg to every gauge block
  of the split S1/S2 bundles on X0(143) x 143a1 — exact ground levels (units 2 pi/A_E), Kuenneth multiplicities, Dolbeault-closedness, the
  flat-torus ladders, critical ratios, and a lower bound on the Morse index. Reproduces, from arithmetic alone, the handoff's S2 wall numbers
  (279 negative directions, 90 Dolbeault-closed, 4 + 56 massless Higgs coefficients) and gives S1's wall 99 / 9. No split background is
  polystable at any r. Legend `surface_block_census`; `tests/test_v0331.py`; notes item 8.

## SM interaction atlas and repository organisation (25 September 2026, commit ad457373)

- Add `mtft.interactions`: a pinned, offline reference containing all 153
  vertices of MadGraph's stock tree-level SM model, with source hashes,
  particle/color/Lorentz/coupling data, parameter definitions and its license.
  Explicitly record the light-quark Yukawa omissions and other source assumptions.
- Add searchable HTML, per-vertex SVG and JSON exports, static declarative UFO
  ingestion, and a hash-bound MTFT evidence ledger. Initial mappings are all
  `unmapped`. Promotion requires references; amplitude-test records require a
  complete context and passing numerical comparison. No MTFT amplitudes are
  claimed or computed by this feature.
- Add an optional MadGraph process card and a reference reproduction script.
- Register the atlas in Legend with an external `GIVEN` source root. External
  SM inputs are not described as quantities derived from MTFT arithmetic.
- Move all 67 historical `CHANGELOG_*.md` files into `docs/changelog/`, preserving
  their bytes and release history; update current navigation and sdist inclusion.
- Replace the obsolete module/test-count README with installation instructions,
  a capability map covering the current tree, research-status links, the GP
  runner, visualization entry points and separate routine/long-run guidance.
- Add a reproducible full module index and focused atlas regression coverage.

The release version remains 0.33.0 until the next release is chosen. No
historical numerical results, frozen study artifacts or physics engines change.
