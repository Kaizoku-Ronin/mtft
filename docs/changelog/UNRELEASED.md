# Unreleased — SM interaction atlas and repository navigation

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
