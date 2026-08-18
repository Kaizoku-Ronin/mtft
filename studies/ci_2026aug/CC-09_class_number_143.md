# CC-09 (append-only) — h(-143) = 10, not 7 (AG Pr 7.8.1 Minkowski slip)

**Filed:** 2026-08-18, v0.16.0 integration, first raised as a **proposal** in
`X0_143_CI_C_REPORT.md` sec.5 (wave files byte-preserved).
**Class:** [EXACT] — independent reduced-forms enumeration by the auditor,
agreeing with the wave's own correction (`ci_verify_kimi.py`,
`ci_verify_kimi.json`).

## The claim being corrected

Arithmetica Generale **Pr 7.8.1** states that the class number of
Q(sqrt(-143)) is **7**, that j(tau_0) is an algebraic integer of degree 7,
and supports this with "the Minkowski bound is sqrt(143)/pi = 3.81, so we
check primes 2 and 3" and "the discriminant -572 has 7 reduced forms".

## Correct values

    h(-143) = 10     (class group cyclic of order 10)
    h(-572) = 10
    (auditor's enumeration also: h(-11) = 1, h(-44) = 3, h(-52) = 2)

## Diagnosed mechanism

The Minkowski bound for an imaginary quadratic field is
(2/pi)·sqrt(|d|), not sqrt(|d|)/pi — a factor of 2. The text used 3.81
where the correct bound is 7.61, tested only p = 2, 3, and undercounted
10 as 7.

## Downstream

- The degree of j(tau_0) over Q in Pr 7.8.1 is **10**, not 7. Any claim
  keyed to a degree-7 CM minimal polynomial at the discriminant -143 CM
  point needs re-checking.
- The fixed-point count #Fix(W_143) = h(-572) + h(-143) = 20 is
  **unaffected**: it already uses the corrected values and lands on the
  Riemann-Hurwitz answer (auditor: both routes recomputed, 0/4/20).
- The CI-A reading "7 = h(-143)" of the factor 7^2 in
  a = -7^2·13·1957^2 is dead. The disc^2 = 1957^2 factor stands as
  flagged-not-claimed (AG-D5); 7^2 and 13 remain unexplained and are not
  claimed.
