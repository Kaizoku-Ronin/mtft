# mtft v0.25.1 — correction release (2026-09-03)

## CC-19 — alpha^-1 comparison target was CODATA 2018
`falsify` predictions #1 and #2 compared against 137.035999084(21) (CODATA 2018).
The current recommended value is CODATA 2022: 137.035999177(21) (NIST wallet
card 2022).  Now `constants.ALPHA_INV_CODATA2022`; the 2018 value is retained as
`ALPHA_INV_CODATA2018_RETIRED`.  Effect on registered verdicts: prediction #2
(3-term) 12.1 sigma -> 7.7 sigma, status unchanged (FAIL on sigma, PASS on its
pre-registered 1 ppm band).  Bookkeeping note, not a registration: the
`_MTFTGauge.alpha_inv_4term` property (137.035999165) sits at 3.9 sigma from the
2018 value and 0.55 sigma from the 2022 value.  It is NOT added to the falsify
engine here — a prediction may only be registered with its derivation date
documented relative to the CODATA 2022 release (May 2024); that is Kimi's call.
Found by Sol (Wave M1 handoff).

## CC-20 — LN_MONSTER truncated
`constants.LN_MONSTER` was the 9-decimal literal 124.126423366.  It is now
computed from `MONSTER_ORDER_FACTORIZATION` (124.12642336632464); the truncated
literal is retained as `LN_MONSTER_TRUNCATED_RETIRED`.  `alpha_inv_monster` is
unchanged at the ppm level and remains rejected (3.54 ppm, ~2.3e4 sigma).
Found by Sol (Wave M1 handoff).

## Tests
`tests/test_cc19_cc20.py` (3 tests).
