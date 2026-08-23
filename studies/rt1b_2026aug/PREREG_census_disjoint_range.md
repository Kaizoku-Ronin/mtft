# Pre-registration — RT-1B.6, holdout range 160 < N <= 400

Filed 2026-08-18. **FROZEN. Nothing in this range has been computed.**
Runner: `rt1b6_disjoint_range.gp`. The 22 levels with N <= 160 are the
discovery sample and are excluded by the runner's own range check.

## H1 — with its eligibility domain stated

    N = 3q,  q > 3,  and mu_2 : Sym^2 H^0(K) -> H^0(2K) surjective
      ==>  U_q ( H^0(2K) )  contained in  H^0(2K).

The surjectivity clause matters: the framework *skips* levels where mu_2
fails (Gate 2 treats hyperellipticity as an exclusion, not a fault), so
without it a future skipped level would leave H1's status ambiguous.

Adopted in Sol's stronger invariance form rather than as r_cusp = 0,
because the discovery data show the whole signature vanishing, (0,0), not
merely its cusp part.

## Q1b — an intrinsic test set, replacing "every good prime"

The earlier wording said "every good prime l" while the harness tested
only the first three. Those are different hypotheses and the mismatch is
corrected here by fixing the set in advance:

    for every good prime l <= L with L = 20, Lambda_l = (4, nu_2 + nu_3).

On the smoke range this yields ~7 primes per level rather than 3. Finite
sampling will not be reported as if it tested a universal statement.

## Q3 — unchanged

When nu_2 + nu_3 > 0, the smaller level prime reaches the elliptic
quotient and the larger does not.

## Not pre-registered, deliberately

Q2 and Q4 are dead and are **not** computed by the runner — no variables,
no output. The r_cusp(U_p) = 1 observation at N = 55, 77 (both q = 11) is
two data points and is left alone.

## Gates — all hard errors

Gate 1 (b_N = 4 + nu_2 + nu_3), Gate 3 (dim S_4^(2) = d4 - 4) and Gate 4
(r_cusp in [0,4] and r_ell in [0, nu_2+nu_3] for **every** signature,
good and bad) call `error()`. Gate 2 remains a logged skip.

**Vacuous-truth guard.** An earlier draft printed HOLDS for all three
hypotheses on a range containing zero levels. The runner now reports
"NO TEST CASES -- vacuous, not a pass" with the case count n, so an empty
or over-filtered range cannot be read as confirmation. Verified on an
empty range.

## Decision rules

- One counterexample kills a hypothesis. Report; do not repair.
- If H1 fails, no weaker version is to be proposed from this range.
- A gate error means my construction is wrong: halt, re-derive, do not
  reinterpret output.
- Run once.

---

## Pre-run clarification, appended 2026-08-18

**No holdout data has been observed.** This appends the sampling frame so
that the prose and the runner are logically identical; it changes no
hypothesis.

> **Sampling frame.** N = pq with odd primes 3 <= p < q and
> 160 < N <= 400. Levels with genus < 3 are excluded. The canonical-leakage
> hypotheses H1, Q1b and Q3 apply only to levels that pass Gate 2
> (mu_2 surjective); levels failing it are logged as
> INVALID_FOR_CANONICAL_LEAKAGE and are neither support nor
> counterexamples.

This makes Q1b's and Q3's domain as explicit as H1's already was.
