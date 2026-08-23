# RT-1B freeze note — corrections, a theorem candidate tested, two bugs

Filed 2026-08-18. Everything below is discovery-range or structural work;
**the holdout range 160 < N <= 400 has not been touched.**

## The four corrections, done

1. **Stale comment.** `rt1b5` printed "the identity is not expected to
   hold" for U_p immediately before demonstrating that it does. Replaced
   with the actual statement: a_1(U_p F) = a_p(F) holds for level primes
   too; what fails is transporting it through W_d, since U_p does not
   commute with W_d.
2. **Q1 -> Q1b.** "Every good prime" is now "every good prime l <= 20",
   an intrinsic set, fixed before the run. ~7 primes per level, not 3.
3. **A separate frozen runner.** `rt1b6_disjoint_range.gp` computes only
   H1, Q1b, Q3 over 160 < N <= 400. No Q2/Q4 variables, no discovery
   output. The discovery census stays untouched as a historical artifact.
4. **H1's domain** now states mu_2 surjectivity explicitly.

## Sol's theorem candidate — tested, holds 22/22

    im( cusp leakage of U_p )  contained in  span{ e_d : d | N, p does not divide d }

Checked directly by testing whether the rows d = p and d = pq of the
4 x dim matrix are identically zero, at every valid level N <= 160.
**True at all 22.** For N = pq the right side is 2-dimensional, so this
would give **r_cusp(U_p) <= 2 for every level prime of every squarefree
two-prime level, with no census at all** — exactly as Sol argued.

That relocates the question, as intended. The zero coordinates become
structural; the content moves to the surviving pair:

| | surviving pair independent? |
|---|---|
| **U_p** (smaller prime) | independent (rank 2) at every level **except N = 55 and 77**, both with q = 11 |
| **U_q** (larger prime) | **never** independent — rank 1 for p >= 5, rank 0 for p = 3 |

So the larger level prime's two surviving cusp functionals are always
dependent, 22/22. That is a sharper statement than the old rank tables and
it is the natural thing to try to prove alongside the vanishing support.

## Two bugs of mine, both caught by contradiction

**`matconcat([A; B]~)` does not stack rows.** It built a **1 x 182**
matrix instead of 2 x 91, so the "surviving rank" column was rank-1 by
construction. Caught because it reported rank 1 for U_11 at N = 143 where
RT-1B.5 had already established rank 2. The vanishing-support checks were
unaffected (direct row tests), and the U_q column *appeared* correct by
luck: a 1 x n matrix has rank 0 or 1, which happened to be the right
answers. Fixed with explicit `matrix(2, n, ...)` construction; the
corrected ranks now agree with the census signatures at every level —
two independent routes to the same numbers.

**Vacuous truth reported as a pass.** A range with zero levels printed
HOLDS for all three hypotheses. The runner now prints "NO TEST CASES --
vacuous, not a pass" with the count n. Verified on an empty range. This
would have been a genuinely dangerous failure in a holdout run.

## RT-1B.5 relations are now gates

`rank = 2, rows d=11 and d=143 zero` for U_11; `rank = 1, rows d=13 and
d=143 zero, c_1 - c_11 = 0`, plus a nonzero-image witness, for U_13; and
`rank = 4` for T_2, T_3, T_5, T_7. All PASS.

**Normalization recorded**, per Sol: c_d(F) = a_1(F | W_d) with W_d from
`mfatkininit`. The vanishing support and the dimensions are
normalization-invariant; the coefficient 1 in c_1 = c_11 is pinned to
that choice and must not be quoted as invariant.

## State

| | |
|---|---|
| B.1–B.5 | established, as Sol's summary has them |
| vanishing support | tested 22/22, theorem candidate |
| **B.6** | **pre-registered and frozen, NOT RUN** |

Next action is a single execution of `rt1b6_disjoint_range.gp`. Whatever
it returns is out of sample.
