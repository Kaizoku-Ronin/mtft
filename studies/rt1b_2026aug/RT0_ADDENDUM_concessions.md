# Addendum to RT-0 — corrections accepted, and what they close

Filed 2026-08-18 after Sol's audit of the RT-0 note. Script:
`rt1_cm_and_ghost.py`.

## Conceded, no defence

**1. "Totally real Hecke fields, therefore CM is excluded" is false.**
Sol's counterexample settles it: a CM newform can have coefficient field
Q outright (32a1). CM is the existence of a quadratic chi with
f (x) chi = f — a statement about twisting, not about the field being
real. CC-01 settles reality and says nothing about CM. My sentence
conflated two different properties and should be struck from the RT-0
note.

**2. Sym^2 gives A_1 in a 3-dimensional representation, not A_2.**
Sym^2 : SL_2 -> GL(Sym^2 C^2) has image PGL_2 = SO_3; the Lie algebra is
sl_2. A Gelbart-Jacquet lift produces a GL_3 automorphic representation
whose symmetry originates from A_1. Rank 3 does not imply type A_2.

This is the same error I had warned against one paragraph earlier —
"8 + 3 + 1 = 12 establishes nothing" and then a degree-3 L-function
offered as a route to su(3). Dimension is not algebra, including when
the dimension is the rank.

**3. "The geometric shadow of exactly that lift" is too strong.**
The canonical multiplication mu_2 : Sym^2 H^0(K) -> H^0(2K) and the
Gelbart-Jacquet lift both use a symmetric-square functor; no functorial
bridge between them has been established. By AG-D5 that is a proximity
claim with no mechanism. Sol's phrasing is correct: it *motivates* the
Sym^2 experiment, it does not realise it.

**4. "They cannot carry a coherent representation" over-claimed.**
What was shown is V_ghost | <W11,W13> = chi_-+ ^(+)3, a multiplicity-three
isotypic piece. That kills "3 generations from 3 AL directions"
permanently, which was the point, but it does not preclude some other
commuting algebra acting on the multiplicity space M_3.

## RT-0 itself: upgraded, with a methodological note

Sol's argument is right and it lifts RT-0 from Cert to **Pr**. For an
integer matrix, rank_{F_p} <= rank_Q, so dim ker_Q <= dim ker_{F_p} = 1;
scalars give dim stab_Q >= 1; hence dim stab_Q(I_2) = 1 exactly, and
dim Lie Stab_{PGL_13}(X_0(143)) = 0.

Worth recording *why* this worked, because it generalises: a mod-p rank
computation becomes a characteristic-zero certificate exactly when it
meets a known bound from the other side. Here the mod-p kernel hit the
scalar floor, so the two inequalities closed. Had it returned 5, we would
have known only dim ker_Q <= 5. One prime sufficed for the mathematics;
the second prime and the rational-normal-curve control were guarding the
implementation, not the theorem — which is the right division of labour
but should be stated as such rather than presented as two-route
verification.

## What the corrections let me close now

**(a) f1 is not CM** — computed, not inferred. Reading a_p off the shipped
e1 = 72 f_143a1: of the 34 primes p <= 140, exactly **2** have a_p = 0
(p = 2 and p = 83), density 0.059. CM forces density 1/2.

**(b) The CM audit for f2, f3 has exactly two candidate fields.** A
weight-2 CM newform with CM by disc -D has level divisible by D. At
squarefree level 143 the divisors are 11, 13, 143; of these -13 = 3 mod 4
is not a fundamental discriminant (the field Q(sqrt(-13)) has disc -52,
which does not divide 143). So the only possibilities are

    Q(sqrt(-11))   and   Q(sqrt(-143)).

RT-1's endomorphism audit is therefore two explicit fields to rule out,
not an open-ended search. Note this makes CC-09 load-bearing: h(-143) = 10,
not 7, and the second candidate field is exactly that one.

**(c) RT-1B — Hecke does not act where the ghost 3-space lives.** The
obvious candidate for "what acts on M_3" is the Hecke algebra. It supplies
nothing:

| quantity | value |
|---|---|
| dim H^0(2K) | 36 |
| dim span(T_2 H^0(2K)) | 36 |
| dim span(H^0(2K) + T_2 H^0(2K)) | **40 = dim S_4(Gamma_0(143))** |

T_2 carries H^0(2K) to a *different* 36-dimensional subspace; together
they fill S_4. Hecke acts on S_4, not on the doubly-vanishing subspace
where the ghost directions live. So M_3 inherits no Hecke action by this
route. Sol's question survives, but it must be re-posed inside S_4, or
with a commuting algebra that is not Hecke.

## RT-2C: where A_2 would have to come from

Accepting that Sym^2 lifts stay A_1, the question is where type A_2 Hodge
structures arise at all. The classical source is **order-3 cyclic covers**
— Picard curves y^3 = f(x) and the Deligne-Mostow families, whose
Jacobians carry Z[zeta_3] actions and whose monodromy sits in U(2,1).
A_2 comes from a Z/3 structure, not from rank-2 data raised in rank.

The arc has already closed that route on this curve, twice over:

* **X0(143) is not trigonal** — certified by `gate_generation`, rank 395
  excluding it via Enriques-Babbage. No degree-3 map to P^1.
* **Aut(X0(143)) is the Atkin-Lehner group, order 4** — no element of
  order 3, so no cyclic triple-cover structure.

So A_2 is not hiding in X_0(143) waiting for a cleverer construction. It
requires changing the ambient group — U(2,1), Picard modular surfaces —
which is an architectural move for the programme, not another wave on the
present object. That is the honest cost estimate, and it is worth having
before the next arc rather than after.

## Revised RT ordering (adopting Sol's, with RT-1B inserted)

| item | content | status |
|---|---|---|
| RT-0 | stab(I_2) = G_m; zero infinitesimal projective symmetry | **Pr** |
| RT-1 | endomorphism audit; two candidate CM fields; A_1^12 prediction | scoped |
| RT-1B | what acts on M_3 — Hecke ruled out; re-pose in S_4 | **partial negative** |
| RT-2A | Sym^2 f_i as certified GL_3 lifts | pre-register |
| RT-2B | verify their monodromy is A_1 in the 3-dim rep (positive control) | pre-register |
| RT-2C | genuinely A_2 data — requires leaving GL_2/modular curves | scoped above |
| RT-3 | hypercharge, chirality, anomalies — blocked until a gauge algebra exists | blocked |
