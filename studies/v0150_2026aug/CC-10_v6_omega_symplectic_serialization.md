# CC-10
## v6 period JSON: derived field `Omega_symplectic_13x26` serialization defect

STATUS: append-only correction. Primitives unaffected; derived field retained
for provenance, permanently rejected by the API.

WHAT IS WRONG. In `x0143_period_data_v6.json`, 15 of the 338 entries of the
derived field `Omega_symplectic_13x26` are corrupt: entries of true magnitude
~1e-58..1e-60 are stored as their mantissa with the exponent dropped (O(1) to
O(10) values). Max discrepancy vs the reconstruction Omega_cusp @ S is
9.956824444577826731 at (row 13, col 13) — literally the orphaned mantissa of
a numerically-zero entry printed by GP as `9.956824444577826731 E-60`.

MECHANISM (Cert, token-level). PARI/GP prints small reals with a space before
the exponent (`mantissa E-60`). The historical text-to-JSON exporter split on
whitespace and took a fixed field as the value, dropping the exponent token.
Verified: the set of corrupt JSON entries is exactly the set of raw-text lines
containing a split exponent token (15 = 15, biconditional), and the raw GP
text itself, parsed exponent-aware, matches Omega_cusp @ S on all 338 entries
to max 6.4e-50.

WHAT IS NOT WRONG. The primitive fields `Omega_cusp_13x26`, `E_intersection`,
`Q_intersection_inverse`, `S_symplectic`, and the frozen `tau_13x13` are
clean. Reconstructing Omega_sym := Omega_cusp @ S and tau := A^{-1}B replays
the frozen tau to ~6e-50 (dps=50); first Riemann bilinear residual ~2.9e-50;
tau symmetry ~2.9e-50. The 19 PASS / 1 FAIL status of the v6 run (C9b) is
untouched.

BLAST RADIUS. Zero known consumers: repo-wide grep finds no reader of
`Omega_symplectic_13x26` outside `mtft.periods.core.legacy_omega_symplectic`
(forensic accessor, v0.21). No downstream number changes.

REMEDIATION (shipped in v0.21.0). `mtft.periods` never reads the derived
field; `omega_symplectic()` is reconstructed at call time as Omega_cusp @ S;
`gate_period_reconstruction` asserts the legacy discrepancy REMAINS > 1 so a
silent "fix" of the frozen artifact cannot pass unnoticed. The defective
field stays shipped verbatim for provenance.

LESSON (recurring class). Same family as the v0.15-era frame lesson and
MANIFEST omissions: derived artifacts must be re-derived at call time from
primitives; text-format parsers must be tested against the printer's actual
token grammar (GP's space-before-exponent) — a parser test with a synthetic
`1.23 E-45` line would have caught this at export time.
