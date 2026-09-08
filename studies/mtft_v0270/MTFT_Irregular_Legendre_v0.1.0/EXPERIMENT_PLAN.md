# Irregular primes in square intervals: analysis plan

Written on 2026-09-06 before computing the new labels above 1,024.
The prior conversation already inspected every odd prime below 1,024.
This is an exploratory experiment, not an external preregistration.

## Scope

- Census: every prime in `1 < p < 101^2 = 10201`, covering square intervals
  `n = 1,...,100`. The exceptional prime 2 has no regularity label.
- An odd prime is irregular exactly when some even `k`, `2 <= k <= p-3`,
  satisfies `B_k = 0 (mod p)`. Record every such index, not just a Boolean.
- The full exploratory comparison uses `n = 6,...,100`: 37 is the first
  irregular prime, so the earlier intervals are kept in the census but are
  excluded from label-exchangeability comparisons.
- A separate comparison uses the previously unexamined range `n = 32,...,100`.
  Each comparison estimates and conditions on its own observed label totals.
- Do not expand the range or alter the tests in response to a p-value.

## Exact arithmetic

Compute Bernoulli residues in F_p by inverting the truncated series
`A(t) = sum_{j=0}^{p-3} t^j/(j+1)!`. Its inverse has coefficients `B_k/k!`
through degree `p-3`. All denominators are invertible modulo p in this range.
Use integer polynomial convolutions with an explicit int64 overflow bound.
Save all even Bernoulli residues, all witnesses, and their indexing.

Validation: multiply each reconstructed inverse by A; compare all residues
for primes below 1,024 with exact rational Bernoulli numbers; independently
check every positive witness and a fixed sample of negative entries using
`B_k = (sum_{a=1}^{p-1} a^k / p) (mod p)`, evaluated modulo p^2. Check prime
enumeration by trial division, and reconcile overlapping prime records with
the previous census when it is supplied.

## Null models

Every randomization leaves prime values, interval membership, and the total
number of primes in every interval unchanged. It never creates or moves a prime.

1. Global: uniformly permute regular/irregular labels within the selected
   analysis range, preserving its total number of irregular primes.
2. Stratified: uniformly permute labels within each combination of residue
   modulo 30 and the fixed square-index band `floor((n-6)/20)`; preserve
   irregular totals in every stratum. This conditions on residue composition
   and broad variation with size. The model is a robustness comparison, not a
   theorem about how irregular primes are generated.

Use 49,999 randomizations per model/range and NumPy PCG64 with base seed
20260906. Store the actual randomization statistics and interval envelopes.

## Three comparison statistics

1. Number of intervals without an irregular prime: upper-tail test for excess.
2. Longest consecutive run of such intervals: upper-tail test for excess.
3. Difference between the irregular fraction among primes in thin intervals
   and the fraction among primes in the other intervals: two-sided test.

Define thin intervals without consulting regularity labels: take the bottom
quarter (floor, ties broken by increasing n) of the selected intervals ranked
by `m_n / ((2n+1)/log(n*n+n+1/2))`, where `m_n` is their actual prime count.
The logarithmic expression is only a density reference for this ranking;
it is not assumed to be an accurate theorem for each short interval.

Compute the exact conditional expected zero probability by the hypergeometric
formula, or a product of those formulas across strata. Center the two-sided
statistic at its analytic null expectation, not the Monte Carlo sample mean.
Include ties and use `(1 + tail_count)/(1 + repetitions)` for Monte Carlo
p-values. Apply Holm correction across all 12 model/range/statistic comparisons.
Report effect sizes, null ranges, and raw and adjusted p-values. No claim of
randomness or independence follows from a failure to detect a departure.

## Interpretation boundaries

An irregular-free interval can contain regular primes. Prime-count conditioning
makes every known empty/nonempty status identical in every shuffle: these tests
assess label organization, not Legendre's conjecture itself. Infinitude of the
irregular family does not establish local coverage or infinitude of arbitrary
intersections. This artifact provides exact finite arithmetic and exploratory
statistical evidence; it is not a Lean-checked proof.
