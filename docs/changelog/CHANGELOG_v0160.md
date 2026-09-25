# v0.16.0 — the canonical-ideal arc (draft 2026-08-18)

One new study bundle (`studies/ci_2026aug`, 17 wave files byte-preserved
plus auditor artifacts), two new corpus corrections (CC-08, CC-09), two
API-hygiene/ledger fixes, and a full independent exact audit of the
2026-08-16 canonical-ideal wave.

## The arc (all claims auditor-verified exact; see bundle README)

* **Canonical embedding.** X0(143) ↪ P^12, degree 24, dim I_2 = 55 =
  (g−2)(g−3)/2 (Petri); quadrics generate the canonical ideal; not
  hyperelliptic/trigonal/plane-quintic. Pre-registered P1–P9 all hit.
* **Atkin–Lehner adapted basis.** W_11, W_13 simultaneously diagonal;
  det B = −1078272; I_2 graded (26, 5, 4, 20) with support confinement;
  the AL grading is the whole of the available sparsity (honest negative:
  LLL vs HNF tradeoff is intrinsic).
* **Projection table.** Nine projections computed by exact ranks, six
  matching Riemann–Roch on the quotients; three honest excesses
  (newspace 3, f2+f3+old 2, f1+f2 1) — the newspace deficiency: a
  3-dimensional piece of H^0(2K) is reachable only through the level-11
  ghost (the newforms have no (−,−) sector).
* **Descent.** π_*K_X = four eigen-line-bundles on E = X0(143)* = 143a1
  of degrees (0, 6, 5, 1), two derivations (sectors and ramification).
  Ten bundle-rank tests 10/10. L_-- = O_E(P), P = (4,−7) a generator of
  E(Q) (rank 1, trivial torsion). Jac(X0(143)/W_143) ~ 143a1 × 11a1,
  32/32 primes (Route 2, replayed with both sides independent).
  Quotient genera 7/6/2, X* genus 1; #Fix(W_Q) = (0, 4, 20) by both
  Riemann–Hurwitz and the class-number formula; W_11 acts freely.
* **Modular degree of 143a1 = 4** — exact deduction from certified inputs
  (Hom(J0(143), E) = Z from the certified charpoly factorization;
  deg π = 4; deg φ₀ = 1 excluded by genus).

## Corpus corrections (append-only ledger)

* **CC-08** — AG Pr 3.7.5: the Atkin–Lehner split of the f2/f3 Galois
  orbits is false. W_p eigenvalue = −a_p ∈ {±1} is rational, hence
  constant on orbits. Exact auditor replay: q4 uniform (−,+), q6 uniform
  (+,−); S_2 joint sectors (1, 5, 6, 1).
* **CC-09** — AG Pr 7.8.1: h(−143) = 7 is false; h(−143) = h(−572) = 10
  (Minkowski bound is (2/π)√|d|, not √|d|/π). j(τ₀) has degree 10.
  #Fix(W_143) = 20 unaffected. The "7 = h(−143)" reading of the CI-A
  coefficient factor 7² is dead; disc² = 1957² stays flagged-not-claimed.

## Fixes

* `modular_curve.HomologyData.intersection_matrix` — docstring now states
  the field ships the template block-J, not the computed pairing of
  H_1(X_0(N), Z) (zero internal consumers; API hygiene only).
* `legend.py` — dropped the stale "(pending drop)" pointer for
  `eisenstein_congruences.py` (shipped in `studies/promotion_2026aug/`).

## Auditor artifacts

* `studies/ci_2026aug/ci_verify_kimi.py` (+ `.json`) — self-contained
  stdlib-only exact replay of the full arc: Cremona-endpoint W_Q
  construction calibrated against `hecke.star_involution()`, certified
  integer ranks (two-prime mod-p + Bareiss minor), exact curve arithmetic
  over Q(i), and the 32-prime Route 2 with independent point counts on
  both sides.
* C3 (h(−572) = 13) is the wave's own honest miss, recorded per policy;
  tracing it produced CC-09.
