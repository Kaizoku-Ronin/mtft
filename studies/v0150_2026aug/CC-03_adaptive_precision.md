# CC-03 (append-only) — precision could be lowered silently

**Filed:** 2026-08-12, v0.15.0 integration. **Class:** [Pr].
**Source:** `INTEGRATION(1).md` CORRECTION CC-03 section (byte-preserved in this
directory); fix landed in `src/mtft/curvature.py`, tests in `tests/test_area_geometry.py`.

## Failure mode

`finite_atom_curvature` took a fixed floor of 0.7*beta + 45 digits and honoured a
caller's lower `dps`. Both halves were unsafe. The required precision is
**support-dependent, not a function of beta**: when the support contains a
collinear triple the leading metric triangle vanishes, det g drops to the next
triangle product, and the cancellation deepens far past any beta-based estimate.
For `{1,2,4,8}` at beta = 80 the old floor of 101 digits returned about
**-2e+19 instead of 1/4** — silently, no warning.

## Fix

The precision loop is now adaptive: the evaluation doubles its working precision
until the result stabilizes (rtol 1e-20), and a caller-supplied `dps` may only
*raise* the working precision, never lower it. Regression tests pin the
`{1,2,4,8}` beta = 80 case to the exact value 1/4.

## Auditor note

The label CC-03 is unchanged from the wave's internal numbering — no collision
existed for this slot.
