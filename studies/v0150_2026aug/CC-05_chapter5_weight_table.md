# CC-05 (append-only) — Chapter 5 weight table

**Filed:** 2026-08-12, v0.15.0 integration. **Class:** [Pr], machine-verified by two
independent routes agreeing to 1e-12. **Source:** `Tano_Weights_Fact_Ledger.md`
§12.1 (byte-preserved in this directory).

## The misstatement

MTFT_Chapter5's worked example and Table 1 give:

    w_6 ~= 1.060,  w_8 ~= 1.040,  w_9 ~= 0.732,  w_10 ~= 1.036

## Correct values

    w_6  = 1.011404  = (4 log 2 + 3 log 3)/6
    w_8  = 0.953077  = 11 log 2 / 8
    w_9  = 0.610340  = 5 log 3 / 9
    w_10 = 0.898720  = (6 log 2 + 3 log 5)/10

## Verification (two routes, agreement 1e-12)

(i)  direct divisor sieve of w_n = sum_{d|n} (log d)/d to N = 2e5;
(ii) closed forms via f = sigma * Lambda, w_n = f(n)/n.

Note the Chapter 5 prose lists the correct summands for w_6 but sums them wrong.
Prime and prime-power entries n = 2, 3, 4, 5, 7 in the table are correct.
Append-only annotation per protocol; nothing rewritten.
