# RT-1B.7 harness design — for the NEXT holdout, not this one

**B.6 is frozen and is not to be altered.** This design applies to
whatever holdout comes after it. Recording it now so the lesson is not
lost, and so it cannot be mistaken for a change to the running experiment.

The weakness B.6 exposed: its outer loop starts at p = 3, so an
interrupted run yields a prefix sorted by the very variable H1 turns on.
No damage was done — the runner emits verdicts only after the loop, and
the prefix was refused as evidence — but an interrupted holdout should
leave an unbiased partial sample, not a biased one.

Sol's design, adopted:

1. **Precompute and freeze the eligible level list** before execution, as
   a separate hashed artifact.
2. **Fix a deterministic, hypothesis-blind order.** Not "increasing N",
   which still correlates with genus and cost; use a frozen permutation

       pi = sort_by( SHA256( study_id || N ) )

   Reproducible, known before any data is observed, and uncorrelated with
   p, q, genus, or any hypothesis variable.
3. **Checkpoint each completed level atomically** — one file per level,
   containing the level's inputs, its signature, and a hash of the runner.
4. **Emit no global verdict until every eligible level completes.** B.6
   already does this and it is the property that saved the partial run
   from being readable as a result.
5. **Restart skips only levels whose checkpoint hash matches** the current
   runner hash, so expensive completed levels are never recomputed and a
   silently edited runner cannot reuse stale checkpoints.

Property 2 is the new one. Properties 4 and 5 together mean an interrupted
holdout can be resumed rather than restarted, which is what made Attempt 1
expensive to lose.
