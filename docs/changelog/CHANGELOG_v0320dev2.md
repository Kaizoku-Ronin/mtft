# mtft 0.32.0.dev2 — Wave B of the v0.32.0 plan: exact engines (2026-09-18)

New pure-function modules extracted from the frozen handoff studies (INT-02…INT-06, INT-09): `surface.spin_circle`
(rank-one local-system cohomology on the spin circle bundle, fiber-acyclicity, connection metric, radius potential),
`surface.hopf_geometry` (half-form pencil, Hopf pullback), `surface.crt_dessin` (CRT projective line, dessin, cusp widths,
genus), `surface.hodge_blocks` (P12/P14 projectors, quaternionic gate), `research.anomalies` (anomaly polynomials,
polarised traces, primitive shift matrix), `research.charge_lattices` (kernel, Smith remnants, integer dressings, vector
masses), `research.mode_operators` (index/cohomology certificates, Bochner bound, 6D Yukawa chirality gate).
Tests `test_wave_b_engines.py` use independent constructions (theorems, Riemann–Hurwitz, direct charge sums, the
package's own Manin permutations).  Register: `docs/SM/WAVE_B_ENGINES.md`.  Base: 0.32.0.dev1.  Pin four-way.
