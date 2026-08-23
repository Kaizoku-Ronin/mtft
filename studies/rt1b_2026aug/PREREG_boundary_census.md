# Pre-registration — RT-1B.4 Boundary Census

Filed 2026-08-18 **before** running the enlarged census. Four hypotheses,
stated by Sol from the seven-level table, now to be tested against every
odd squarefree N = pq in range that survives the gates.

## Hypotheses

| # | statement |
|---|---|
| **Q1** | every good prime l satisfies Lambda_l = (nu_infinity, nu_2 + nu_3) — full saturation of both channels |
| **Q2** | for N = pq with p < q, r_cusp(U_q) = 1 always |
| **Q3** | when nu_2 + nu_3 > 0, it is always the **smaller** level prime that reaches the elliptic quotient |
| **Q4** | r_ell(U_p) <= 1 for every level prime |

Current evidence: Q1 on 7 levels, Q2 on 7, Q3 on 3, Q4 on 7. All four are
7-or-fewer-point patterns. None is claimed as a law.

## Gates enforced in code, not by reading output

1. `bN = dim S_4 - (3g-3)` must equal `4 + nu_2 + nu_3`, else the level is
   reported MISMATCH and skipped.
2. `rank(mu_2)` must equal `3g-3`, else `INVALID_FOR_CANONICAL_LEAKAGE`
   (canonical multiplication not surjective — hyperelliptic) and skipped
   before any leakage is computed.
3. `dim S_4^(2)` must equal `d4 - 4` (four cusps, AL free transitive).
4. Every computed rank must satisfy `0 <= r <= bN`, else GATE VIOLATION.

## Decision rules, fixed in advance

- A single counterexample kills the hypothesis. Report it; **do not repair
  the statement post hoc to exclude the counterexample.**
- If Q2 or Q3 fails, the smaller/larger framing is wrong and no
  replacement ordering rule is to be proposed from the same data.
- Levels failing gate 1 or 3 indicate an error in my own formula, not a
  discovery; halt and re-derive rather than reinterpret.

## Also being corrected in this round

- The staged good-pair count was reported as **31**; the artifacts
  reproduce **22**. Verified by parsing the staged scripts. The census
  supersedes the count entirely.
- "The two level primes are never symmetric" is **false** and was
  contradicted by the table printed directly above it: N = 55 gives
  U_5 = U_11 = (1,0) and N = 77 gives U_7 = U_11 = (1,0). Withdrawn.
- `B_N = B_cusp (+) B_ell` implies a canonical splitting that has not been
  established. The correct statement is the short exact sequence
  `0 -> S_4^(2)/H^0(2K) -> S_4/H^0(2K) -> S_4/S_4^(2) -> 0`; the two ranks
  are well defined as successive quotients without choosing a splitting.
