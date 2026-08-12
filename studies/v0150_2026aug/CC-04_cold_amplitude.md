# CC-04 (append-only) — the cold amplitude was wrong

**Filed:** 2026-08-12, v0.15.0 integration. **Class:** [EXACT] closed form; E2-certified
with the area-geometry wave. **Source:** `INTEGRATION(1).md` CORRECTION CC-02
section (internal numbering; renumbered CC-04 here to resolve the collision with
the corpus CC-02 w-series correction).

## Correct value (closed form)

    c = (9 L5^2)/(25 L2 L3) [1 - (9/5) L5/L3 + (4/5) L5/L2]
      = 0.27012646530542495706433719670365...

with L2 = log 2, L3 = log 3, L5 = log 5.

## What it replaces

`CURVATURE["cold_amplitude"]` was `0.270126465305424759517602` — wrong in the
16th digit. The old value came from a six-atom extraction at beta = 200, where
atom 6 contaminates at scale (5/6)^200 = 1.46e-16; the observed deviation
1.98e-16 matches that error model. The **cold core is {1,2,3,5}**, not the first
six integers; the cold rate is 6/5 = 36/30.

## Status

Landed in `src/mtft/curvature.py` (`cold_amplitude()` closed form, retracted
value kept under `cold_amplitude_RETRACTED_CC04`) and pinned by
`tests/test_area_geometry.py`.
