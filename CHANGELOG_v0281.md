# mtft v0.28.1 — CC-24, CC-25 and exact-input hardening from Astra's v0.28.0 audit (2026-09-11)

CC-24 `surface.marked`: readout denominator was the (13,17)-specific 320; now (q+1)² − a_q² via
`readout_coefficient` (reversed sector had read 9α/16).  CC-25 `lattice.metropolis_sweep`: unit
periodic extents rejected before mutation (staple would contain the updated link).  Hardening:
`oldsector` requires exact integer primes/coefficients; `marked.decomposition` rejects Float/complex.
Policy: none of the QCD-07…11 constructions (H36, scalar trace assignment, metric replacements,
potentials, unit scales) is adopted as a default.  See MTFT_v0280_INTEGRATION_ASSESSMENT.md.
Base: live 0.28.0 (SHA-256 3d3b25fe…).  Tests: `tests/test_cc24_cc25.py` (3).  Pin four-way.
