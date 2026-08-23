# RT-1B.4 — Boundary Census. Two hypotheses died.

Filed 2026-08-18. Pre-registered in `PREREG_boundary_census.md` **before**
the run. Script: `rt1b4_boundary_census.gp`.

## Two corrections first

**The good-pair count was wrong.** I claimed 31; the staged artifacts
reproduce **22**. I verified Sol's count by parsing my own scripts:
11 valid good pairs from the cross-level script, 13 from the 143 sweep,
overlap {(143,2), (143,3)}, union 22. "Computed somewhere" is not
"reproducible from the staged certificate", and I was on the wrong side
of that line. The census supersedes the count: **66 good (N, l) pairs
across 22 valid levels, all from one script.**

**"The two level primes are never symmetric" is false**, contradicted by
the table printed directly above it: N = 55 gives U_5 = U_11 = (1,0) and
N = 77 gives U_7 = U_11 = (1,0). Withdrawn.

**Notation.** `B_N = B_cusp (+) B_ell` implied a canonical splitting that
does not exist. Correct statement: the short exact sequence

    0 -> S_4^(2)/H^0(2K) -> S_4/H^0(2K) -> S_4/S_4^(2) -> 0

The two ranks are well defined as successive quotients without a
splitting. Adopted.

## The census

All odd squarefree N = pq <= 160 with g >= 3. Gates enforced in code:
boundary formula, canonical-multiplication surjectivity, dim S_4^(2),
and 0 <= rank <= b_N. **22 valid levels, 66 good (N, l) pairs, 0 gate
violations.** Three levels (33, 35, 39) correctly halted as
INVALID_FOR_CANONICAL_LEAKAGE.

| result | verdict |
|---|---|
| **Q1** every good prime gives (4, nu_2 + nu_3) | **HOLDS** — 66/66 |
| **Q2** r_cusp(U_q) = 1 for the larger level prime | **FALSIFIED** |
| **Q3** smaller level prime reaches the elliptic quotient | **HOLDS** — 9/9 |
| **Q4** r_ell <= 1 for both level primes | **FALSIFIED** |

### How Q2 died

Not marginally. Every level with p = 3 — that is 51, 57, 69, 87, 93, 111,
123, 129, 141, 159, ten of the twenty-two — gives

    U_q = (0, 0)

The larger level prime does **nothing at all**: no cusp leakage, no
elliptic leakage. The seven-level sample simply contained no p = 3 level,
so the hypothesis was formed on a biased sub-family.

### How Q4 died

r_ell = **2** at N = 93, 111, 129 (U_3) and N = 145 (U_5). The bound of 1
held on all seven original levels and fails on four of the fifteen new
ones.

### Q1 is now the strong result

66 good pairs, 22 levels, boundary dimensions 4, 6 and 8, elliptic
channels arising from nu_2 and from nu_3 separately — and every single
good prime saturates **both** quotients exactly. That is the statement
worth trying to prove.

## The pattern I am NOT claiming

The p = 3 split is clean: r_cusp(U_q) = 0 at all ten p = 3 levels and
= 1 at all twelve p >= 5 levels. 22 for 22.

**My own pre-registered decision rule forbids proposing it from this
data** — "if Q2 fails, no replacement rule is to be proposed from the
same data". So it is recorded as an observation, explicitly post hoc, and
goes into a fresh pre-registration to be tested on the **disjoint** range
160 < N <= 400. If it survives there it is worth a mechanism; if it dies
it cost nothing.

Likewise r_cusp(U_p) = 1 occurs only at N = 55 and 77, both with q = 11.
Two data points. Not a pattern.

## Sol's next object, endorsed

Ranks say that something happens; the row space says what. For the cusp
quotient the four functionals c_d(F) = a_1(F | W_d), d | N, are explicit,
so the image of C . U_p restricted to H^0(2K) is computable as an actual
matrix — exact linear relations among c_1, c_11, c_13, c_143 at level 143,
where the partial AL result already says the two lines sit in opposite
W_13 parities. That converts (2, 0) from a statistic into a description.
The elliptic quotient needs a local-functional basis first; that is the
harder half and I have not built it.

## Status

| item | status |
|---|---|
| Q1, good primes saturate both quotients | computation, 66 pairs / 22 levels |
| Q3, smaller level prime reaches elliptic | computation, 9 levels |
| Q2, Q4 | **FALSIFIED, recorded, not repaired** |
| p = 3 split | observation, post hoc, pre-registered for disjoint range |
| cusp row space at 143 | next, computable |
| elliptic local functionals | not built |
