# mtft v0.24.0 — the Eisenstein/cuspidal release

Level-generic boundary arithmetic, the Wave-8 cuspidal machinery packaged,
and the package's first external literature anchors.

## Headline

`mtft.levels` lifts the modular-symbol core off the single curve X_0(143).
v0.23.0 hardwired `LEVEL = 143`, named 11 and 13 in the P^1 validity test,
and baked the triangle and cusp counts in as the literals 56 and 4, which
made *level-universality* — the program's largest standing risk since Wave 7
— untestable in code.  It is now a function call.

The immediate scientific return, **and the correction that matters most**:
the Wave-8 identity **C[2] = E_Eis** is *not* level-universal.

It holds at every squarefree level with two prime factors tested -- 33, 35 and
143 -- where the cusp set is a (Z/2)^2-torsor and E_Eis is 2-dimensional.  It
**fails at N = 105 = 3*5*7**, the first control with three prime factors and
eight cusps, where dim E_Eis = 8, dim C[2] = 6, and the two are distinct
(CC-17).  N = 105 carries the same genus 13 as 143, which isolates cusp
structure from genus as the operative variable.

So the honest fork called in advance resolves in the middle: the identity is
a theorem about the *cusp-torsor* regime, not about Eisenstein ideals in
general, and the Wave-8 SPLIT verdict inherits that conditionality.  What
singles 143 out is not the level but the shape of its boundary.

A second correction (CC-16): N = 15 has genus 1, so H_1 (x) F_2 is itself
only 2-dimensional and E_Eis is the whole space -- the identity holds there
for trivial reasons.  That row is **vacuous** and must not be counted as
support.  `cross_level_control` now flags every row with `informative`, and
the earlier draft of this changelog, which counted 15 as evidence and claimed
universality across "15, 33, 35 and 143", was wrong on both points.

## New modules

**`mtft.levels`** — level-generic Manin model for squarefree N with no
elliptic points.  `level_data`, `manin_model`, `hecke_matrix`,
`cuspidal_hecke`, `boundary_matrix`, `cusp_labels`, `supported_levels`.
Scope is enforced, not assumed: non-squarefree levels, levels with elliptic
points (N = 21 is the smallest), and genus-0 levels raise
`UnsupportedLevelError` with the specific reason rather than returning a
wrong answer.  `manin_model(143)` reproduces `hecke.model()` field by field,
and that equality is a shipped test — it is what licenses reading any
level-generic result back onto the certified 143 work.

**`mtft.cuspidal`** — cuspidal subgroups and Eisenstein torsion.
`cuspidal_group(N)` computes C(J_0(N)) by two independent routes (CRT
spectral projector; boundary-cokernel with no projector at all) and gates
their agreement, per protocol E2.  Also `eisenstein_subspace`,
`eisenstein_kernel_mod2`, `two_torsion_image`, `cross_level_control`,
`charpoly` (Faddeev-LeVerrier, exact), and `lambda2_congruence_scan`.

**`mtft.al_morphology`** — Atkin-Lehner read as maps of surfaces rather than
as operators.  Traces on H_1 give Lefschetz fixed-point counts, which give
quotient genera by Riemann-Hurwitz:

    g = 13  ->  g(X/W_11) = 7, g(X/W_13) = 6, g(X/W_143) = 2  ->  g(X*) = 1

with the full-quotient genus obtained two agreeing ways (Riemann-Hurwitz and
the character formula).  W_11 acts *freely* (trace +2, no fixed points), so
X_0(143) -> X_0(143)/W_11 is an unramified double cover.  All three
involutions act freely on the four cusps, so the cusps form a single free
(Z/2)^2-torsor and the boundary lattice is the augmentation ideal of Z[G] —
the structural reason the 26 -> 29 puncture can carry an obstruction at all.

## External anchors (new)

`test_mazur_prime_level_anchor` checks |C(J_0(p))| = numerator((p-1)/12) at
p = 11, 23, 47 and the computed groups are Z/5, Z/11, Z/23.  These are the
first tests in the package that confront published mathematics rather than
internal consistency.  Composite-level values computed here (N = 15: Z/2 x
Z/4; N = 33: Z/10 x Z/10; N = 35: Z/2 x Z/24; N = 143: Z/10 x Z/420) are
**not** yet checked against Ligozat / Chua-Ling / Yoo and are queued for
audit.

## Corrections (append-only)

**CC-11** `integral_lattice.saturate` — the denominator guard was dead: it
ran *after* `_as_obj` had already coerced with `int()`, which truncates
Fractions silently (`int(Fraction(1,2)) == 0`), so by the time the guard
looked, every denominator was 1.  The check now runs on the caller's matrix
before any coercion, and `_as_obj` raises rather than truncates.

**CC-12 / CC-13 / CC-14** — the same defect, found by auditing for the
family rather than the incident: `solve_in_lattice` truncated a rational
right-hand side, `rational_kernel` truncated a rational input matrix, and
`class_order` truncated while scaling.  All three silently returned wrong
answers on rational input.  Fixed; each ships a regression test that fails
against the old behaviour.  New `InexactInputError` and
`clear_denominators` make the correct workflow explicit.

**CC-15** cusp labelling.  The 2026-08-29 Wave-8 probe labelled the four
cusps of X_0(143) by model index as (1, 11, 13, 143).  The true divisor
labels in model order are (143, 13, 11, 1): index 0 is the cusp infinity,
which has width 1 and c = 0, hence gcd(c, N) = N.  The old names were
reversed by d -> N/d.  Because the Atkin-Lehner divisor-toggle commutes with
that reversal, **no Wave-8 gate is affected** — the group C = Z/420 x Z/10,
its order 4200, the two-route agreement, the SPLIT verdict and the selected
line all stand.  Only the *names* change:

    Z/420 generator   [C_143 - C_1]                  ->  [C_1 - C_inf]
    Z/10  generator   [C_13 - C_1] + 35[C_143 - C_1] ->  [C_11 - C_inf] + 35[C_1 - C_inf]
    AJ orders         70, 60, 420  now attach to  [C_13 - C_inf], [C_11 - C_inf], [C_1 - C_inf]

Found by the new tooling within an hour of its existing, which is the
argument for building tools before running more waves.

## Other changes

* `integral_lattice.snf_transform` — Smith normal form *with* transforms
  (U A V = S), needed for quotient-group generators rather than just orders.
  Cross-checked in tests against the existing independent HNF-ping-pong
  `smith_invariants`.
* `integral_lattice.int_kernel` — saturated integer kernel.
* `test_version_triple` now genuinely reads pyproject.toml and CITATION.cff
  instead of only `mtft.__version__`; through v0.23.0 a stale pyproject or
  CITATION could have shipped undetected.  Verified against a deliberately
  stale file.
* `lambda2_congruence_scan` ships the Wave-8 open conjecture
  a_ell = 1 + ell (mod lambda_2^2) as a standing falsifiable test, reporting
  `counterexamples` explicitly.  Confirmed for all 43 good odd primes
  ell <= 199; unproved in general, since Sturm bounds control congruences
  mod lambda_2, not mod lambda_2^2.

## Tests

717 passing (3 skipped), up from 705.  40 new in `tests/test_v0240.py`.

## Known open / queued for audit

* composite-level cuspidal groups vs Ligozat / Chua-Ling / Yoo (no LMFDB in
  container)
* CM decomposition of the 20 fixed points of W_143 into class numbers; the
  naive h(-143) + h(-572) split is in tension with genus theory
* the X_0(143)/<W_11, W_13> = genus 1 quotient and its relation to the
  level-11 old block is computed but not certified per-form
* elliptic-point levels (N = 21 and friends) need the torsion Manin
  relations before `levels` can cover them

## `mtft.kakeya` — Arf parity, direction sets, finite Kakeya

Leverage extracted from the theta-census structure, certified rather than
asserted.  On the G-invariant theta-characteristic space of X_0(143) (128
characteristics, affine F_2^7, parity radical R of dimension 5):

    Delta(Theta_odd)  = R,   31 = 2^5 - 1 nonzero directions
    Delta(Theta_even) = V,  127 = 2^7 - 1 nonzero directions

Odd spin states are directionally **confined to the radical**; even states are
directionally **saturated**.  The premise is not assumed: the observed 96/32
split forces q|_R = 0, because a nonzero q|_R would be surjective and give 64
odd characteristics rather than 32.  `arf_direction_theorem` checks every step
and reports the forcing condition explicitly.

The result is Kakeya-*adjacent*, not Kakeya: it concerns difference directions,
not the existence of a full line in every direction.  The return dict carries
`is_kakeya_theorem: False` and a test asserts it, to block promotion.

For the genuine article the module ships `besicovitch_set(p)` — tangent lines
to a parabola in F_p^2 — verified to be a Kakeya set and to sit exactly
(p-1)/2 above Dvir's bound binom(p+1, 2) for p = 5, 7, 11, 13.

`crt_direction_bridge(11, 13)` replaces an earlier loose claim.  Z/143 is a
product ring with zero divisors, so Dvir does not apply to it directly; the
correct statement splits first, and is verified bijectively:

    P^1(Z/143)  ~  P^1(F_11) x P^1(F_13),   12 * 14 = 168 = [PSL_2(Z):Gamma_0(143)]

**Deliberately not shipped:** the Atkin-Lehner permutation of P^1(Z/N).  W_Q
has determinant Q, so the naive row action mod N is not well defined on P^1;
a correct version needs the coset-level action.  Omitted rather than shipped
unverified — target for the next wave.

## Fixes in this pass

- Stale version pins in `tests/test_periods_v0210.py` and
  `tests/test_periods_v0220.py` still asserted `0.23.0`, so the suite was red
  on arrival.  Fixed; `test_version_triple` and these pins now agree.
- CC-16: `cross_level_control` gained a non-vacuity flag (`informative`).
- CC-17: `cross_level_control` gained N = 105 and the identity's failure there
  is now asserted by a test, so the negative cannot be silently lost.
