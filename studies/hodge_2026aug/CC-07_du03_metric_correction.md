> **RENUMBERING NOTE (2026-08-12):** this correction was originally filed as CC-02 before the corpus correction ledger existed. The corpus ledger's CC-02 slot is the w-series shift-chain correction (W1_weil_compression_study §1). To keep every correction addressable by a unique id, this du03 metric correction is renumbered **CC-07**; its content below is unchanged. The renumbering is recorded in the v0.15.0 integration notes.

---

# Internal Correction CC-02 — the du03 commutant is metric-dependent

**Class:** CORRECTION + STRENGTHENED RESULT
**Raised:** 2026-08-06, from the Hodge/Satake analysis of H^1(X_0(143))
**Scope:** internal (du03 session only; no published paper affected)
**Regenerate:** `py studies/du03_metric.py`

---

## What was wrong, twice

**First error (du03_commutant.py).** The symmetric commutant of the
transported Hecke clock T_h on the 26-dim harmonic stage was computed as
5-dimensional, and a spin-factor Jordan-algebra reading was built on it,
together with a 1|12|12|1 block structure. The computation was correct;
the *interpretation* was not, because "symmetric" was taken with respect
to the Euclidean (edge-counting) inner product, and T_h is not
self-adjoint there.

**Second error (over-correction).** On noticing this I claimed the 5 was
an artifact and the true answer was 43. That was also wrong — it swapped
one metric-specific answer for another.

## What is actually true

T_h has Euclidean non-normality ||[T_h,T_h^T]||/||T_h||^2 = 0.136, but is
self-adjoint to 5.1e-15 for a Hecke-invariant metric G (the polarization
class). Both counts are correct and answer different questions:

| interaction class | symmetric w.r.t. | dim commuting with T_h |
|---|---|---|
| graph-local (potentials, wells, incidence-built) | Euclidean edge metric | **5** |
| Hodge-natural (periods, intersection form, modular forms) | polarization G | **43** |

cond(G) = 112.3, so the two metrics are genuinely different, not
proportional.

## The strengthened result

The decisive question is whether any interaction is admissible in **both**
senses: S = S^T (Euclidean, so realizable as a graph operator), [S,G] = 0
(so also G-symmetric), and [S,T_h] = 0 (so a rate is well-posed).

**Answer: the joint space is 1-dimensional and contains only scalars.**

A scalar shifts all 26 modes equally and splits nothing. Verified
independent of the choice of invariant metric — G lives in a 43-dim cone,
and 8 randomly sampled positive-definite members (cond 201–425) all give
joint dim 1, non-scalar 0.

**So the du03 obstruction is metric-independent.** No interaction can be
simultaneously graph-local and Hodge-natural and lift the harmonic
degeneracy. This is stronger than the earlier commutator-gate result,
which was itself metric-dependent.

## The reframe this yields

The two-clock problem is a **metric mismatch**, not merely a rate
mismatch. The graph clock Delta_1 is self-adjoint for the edge metric;
the Hecke clock T_h is self-adjoint for the polarization. The geometric
well fails to commute with G (0.165), confirming a graph-local
interaction is not Hodge-natural. The clocks do not merely run at
different rates — they disagree about which vectors are orthogonal.

Anchor count stays at 2, now for a metric-independent reason.

## Retracted from the du03 record

- "the 5 form a non-associative Jordan algebra" — true only in the
  Euclidean metric; not an invariant statement
- "common blocks 1|12|12|1" — metric-specific. In the G metric the
  generic invariant operator has 26 distinct eigenvalues (blocks all 1)
- "43 replaces 5" — over-correction, withdrawn

## Surviving

- T_h integral on the cycle lattice (du02 C2)
- the transported-well no-go (du02 C3 consequence)
- 13 lines / 12 distinct eigenvalues accounting
- anchor count 2
