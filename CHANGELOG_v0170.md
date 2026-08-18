# v0.17.0 — `mtft.canonical` (draft 2026-08-18)

The canonical-ideal arc moves from prose to gated code. One new
subpackage, one new test module, no behaviour change elsewhere.

## `mtft.canonical`

Frozen certified data plus an API, with every pre-registered number
**re-derived from that data at call time** rather than stored as an
assertion. If the shipped q-expansions were altered, the gates fail.

* `_data/` — six wave artifacts, byte-identical to `studies/ci_2026aug/`,
  with `PROVENANCE.txt` carrying SHA-1s that match that bundle's
  `sha1_manifest.txt`. Shipped as package data under `src/`, not under
  `studies/` — the `jacobian_order` pattern, which is structurally immune
  to the glob class of bug fixed in v0.16.1.
* API: `s2_qexpansions`, `adapted_basis`, `adapted_qexpansions`,
  `ideal_basis`, `ideal_by_sector`, `ci_a_quadric`, `monomial_sector`,
  `MONOMIALS`, `COORDINATE_LABELS`, `DESCENT`, `PREDICTIONS`.
* **Sector ordering pinned.** `SECTOR_ORDER = ((+,+), (+,-), (-,+), (-,-))`
  giving S_2 dims (1, 6, 5, 1). The v0.16.0 CHANGELOG quotes the same data
  in the order ((+,+), (-,+), (+,-), (-,-)) as (1, 5, 6, 1). Both correct,
  not independent results. `reorder_sectors()` converts; a test asserts the
  round trip so the two can never silently diverge.

## `mtft.canonical.gates` — eight gates, all recomputed

| gate | re-derives | cost |
|---|---|---|
| `gate_petri` | P1–P4: 13, 91, 36, 55; residual **exactly 0** for all 55 | 0.4 s |
| `gate_sector_grading` | 26/5/4/20, support confinement, residual 0 | 0.1 s |
| `gate_bundles` | ten product-rank tests vs deg L = (0, 6, 5, 1) | 0.2 s |
| `gate_projection` | the nine-row table incl. all three honest excesses | 0.1 s |
| `gate_descent` | eta identity for the ghost lines; bidouble; Riemann–Hurwitz | 0.0 s |
| `gate_ci_a` | Q\* residual 0, confined to (+,+), a = −7²·13·1957² | 0.0 s |
| `gate_generation` | P5–P8: 455, 60, 395 (Enriques–Babbage) | 23 s, `slow` |
| `gate_route2` | Jac(X0(143)/W143) ~ 143a1 × 11a1, prime by prime | 24 s, `slow` |

Ranks use integer elimination mod two distinct primes, per the arc's
policy. Residuals are exact integer comparisons against zero.

`gate_descent` and `gate_route2` need **no external curve data**:
a_p(11a1) comes from the exact eta product η(τ)²η(11τ)², and a_p(143a1)
is read off the shipped e1 = 72·f_143a1.

## Scope note recorded, not papered over

The wave's Route 2 ran **32** primes in [3, 149]. The shipped q-expansions
stop at q^140, so a_p(143a1) is unavailable at p = 149 and the in-package
gate reproduces **31 of 32** (`ROUTE2_PRIMES_EXPECTED = 31`). The 32nd
used external curve data and stays in `X0_143_CI_D_REPORT.md`; it is not
re-derived here. Raising the cap needs deeper q-expansions, not a
different gate.

## Tests

`tests/test_canonical.py` — 16 tests (14 fast, 2 `slow`), including a
provenance test that SHA-1s the six data files against `PROVENANCE.txt`.

## Not done

`mtft.periods` (tau in H_13 still only in PARI/GP output), the disc²
cross-level falsification, and the CC-09 downstream sweep for material
keyed to h(-143) = 7 or a degree-7 CM minimal polynomial.
