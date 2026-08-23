# RT-1B.6 — executing, incomplete, NO VERDICT

Filed 2026-08-18. **Do not read a result from this file.**

## The XMAX mismatch, fixed

Sol is right and it is mine again: I lowered `XMAX` from 160 to 95 while
debugging the row extraction and staged that copy, while the freeze note
claimed 22/22 at N <= 160. The staged artifact could not reach 111, 123,
129, 141, 143 or 159.

`rt1b_vanishing_support.gp` now runs to 160, carries a vacuity guard, and
additionally **hard-errors if the valid-level count is not 22** — so a
silent change in eligibility logic cannot pass as agreement. Rerun clean:

| claim | reproduced |
|---|---|
| 22 valid levels | yes |
| vanishing support at every level | 22/22 |
| U_p surviving pair dependent only at 55, 77 | yes |
| U_q surviving pair never independent | rank 0 at ten p=3 levels, rank 1 at twelve p>=5, never 2 |

Every sentence in the freeze note now reproduces from the attached script.

## Pre-run clarification appended, then hashed

The sampling frame was appended to the frozen pre-registration
(N = pq, odd 3 <= p < q, 160 < N <= 400, genus >= 3, hypotheses apply only
after Gate 2), timestamped as pre-run with no holdout data observed. Then:

    4f551782f5eb7ef4d7a0335879f9ea793fd00084eba55991bbe76dcd6fc325ce  rt1b6_disjoint_range.gp
    063c345d148b17eaf66d985ef31237ace6b7800dc028163bbb841b64baeff22c  PREREG_census_disjoint_range.md
    17ad1938a2cc75c872606a8a6f1660c731f45ab678d4f358ea7ab7c942eeeb0b  rt1b_vanishing_support.gp
    3f9655a16abc3ab7d0f3400f5e41beae0cc454a9e0c21a7c806619d9440b52a9  rt1b5_cusp_leakage_matrices.gp

Full list in `FREEZE_HASHES.txt`. The runner hash was re-verified
immediately before execution and is unchanged.

## Status: Attempt 1 INTERRUPTED after 8 of roughly 45 levels

`rt1b6_PARTIAL_run.log` is the raw partial output, attached for
provenance. It contains **8** completed level rows — 177, 183, 201, 213,
219, 237, 249, 267 — not 7; the log advanced after the first draft of this
note was written. The process is confirmed dead; Attempt 1 is definitively
interrupted.

**Provenance language, corrected.** It is no longer accurate to say the
holdout is untouched. The accurate statement is:

> The hypotheses, eligibility rules, gates, decision rules and executable
> runner were frozen and hashed **before any holdout observation**.
> Attempt 1 was interrupted after eight ordered cases. No code,
> hypothesis, gate or decision rule was changed afterward.

A completed run is therefore **Attempt 2 — an exact byte-identical rerun
of frozen hash 4f5517...25ce**, not "one execution". The pre-registration
remains meaningful because the plan is fixed and unchanged; only the
execution history needs stating.

**No hypothesis is being called.** Two reasons, and the second is the
important one:

1. The runner emits verdicts only after the loop, so a partial run cannot
   produce a spurious HOLDS. Good — but it also means there is nothing to
   report yet.
2. **The loop ordering front-loads H1.** The outer loop is
   `forprime(p = 3, ...)`, so every N = 3q level is evaluated first. The
   visible prefix is therefore *entirely* H1 test cases, while Q1b has
   been sampled thinly and Q3's elliptic levels at p >= 5 have not been
   reached at all. Reporting from this prefix would be reporting a
   structurally biased sub-sample — the exact error that killed Q2, where
   the discovery set happened to contain no p = 3 level.

What is visible is consistent with H1 so far. That sentence is the most I
am willing to write, and it should not be quoted as evidence.

## To complete

    gp -q rt1b6_disjoint_range.gp

Verify the hash first. Report whatever it returns without modification:
one counterexample kills the corresponding hypothesis, H1 may not be
weakened using this range, and a gate error means halt and re-derive
rather than reinterpret.

## Worth fixing before any future holdout

The front-loading is a design weakness even though it did no damage here.
A one-shot holdout runner should either interleave the sampling frame or
checkpoint per level so that an interrupted run yields an unbiased partial
sample rather than a prefix sorted by the very variable one hypothesis
turns on.
