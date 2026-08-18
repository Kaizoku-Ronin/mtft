# CI-A and CI-B: the anomalous quadric and the ghost-only 3-space

**Session 2026-08-16. Both pre-registered before computation
(`PREREG_CI_AB.md`). All ranks exact over Q; Q* residual exactly 0.**

---

## Headline

Both anomalies dissolve into **one** structure, and it is not the one I
expected. The Atkin-Lehner sectors of S_2(Gamma_0(143)) are the spaces of
sections of **four line bundles on the genus-1 quotient E = X0(143)\***,
of degrees **(0, 6, 5, 1)**. Every rank in this study is h^0 of a product
of those bundles, capped by the number of monomials available. CI-B is
fully explained by it. CI-A turns out not to be about f2 at all.

## The descent structure

e1 is (+,+), so both involutions fix it and it descends to E; the f2 and
f3 forms pick up a global sign, which is invisible in projective space, so
their products descend too. Writing pi: X -> E for the degree-4 quotient,
pi_* K_X splits into eigen-line-bundles L_chi with

| chi | (+,+) | (+,-) | (-,+) | (-,-) |
|---|---|---|---|---|
| dim of sector | 1 | 6 | 5 | 1 |
| **deg L_chi** | **0** | **6** | **5** | **1** |

h^0 = deg for positive degree, h^0(O_E) = 1 for the trivial bundle, and
the degrees sum to 12 = chi(K_X). Predicted rank of each product space =
min(h^0 of the product bundle, number of monomials). Ten tests:

| product | monomials | deg | predicted | computed |
|---|---|---|---|---|
| L++ ^2 | 1 | 0 | 1 | 1 |
| L+- ^2 | 21 | 12 | 12 | 12 |
| L-+ ^2 | 15 | 10 | 10 | 10 |
| L-- ^2 | 1 | 2 | 1 | 1 |
| L++ · L+- | 6 | 6 | 6 | 6 |
| L++ · L-+ | 5 | 5 | 5 | 5 |
| L++ · L-- | 1 | 1 | 1 | 1 |
| L+- · L-+ | 30 | 11 | 11 | 11 |
| L+- · L-- | 6 | 7 | 6 | 6 |
| L-+ · L-- | 5 | 6 | 5 | 5 |

**Ten for ten.** And the class dims of H^0(2K) fall straight out:
12 = deg(L+-^2), 6 = deg(L++·L+-), 7 = deg(L+-·L--), 11 = deg(L+-·L-+),
summing to 36. The 12/6/7/11 grading found last session was this all along.

---

## CI-B — the ghost-only 3-space: **fully explained, all predictions hit**

| # | predicted | computed |
|---|---|---|
| B1 | deficiency confined to class (-,+) | yes |
| B2 | per-class deficiency (0, 0, 3, 0) | **(0, 0, 3, 0)** |
| B3 | newform rank in class (-,+) is exactly 4 | 4 |
| B4 | missing 3 supplied by oldspace monomials | yes |

Class (-,+) of H^0(2K) is H^0 of a degree-**7** bundle. Only two channels
reach it: L++·L-+ (degree 5) and L+-·L-- (degree **7**). The newspace has
**no (-,-) coordinate** — the single (-,-) line of S_2(Gamma_0(143)) is
oldspace, exactly as the package constant `AL_DECOMPOSITION` records with
(-1,-1): 0 — so the newspace can only use the degree-5 channel, and only
4 of its 5 dimensions (the 5th is the (-,+) oldform). Rank 4, deficiency 3.

**B5 answered decisively.** Which oldspace direction supplies the missing 3:

| added to e1·f2 (rank 4) | rank |
|---|---|
| + y1·y8 — the (-,+) ghost | 5 |
| + y13·f3 — the **(-,-) ghost** | **7** |
| + both | 7 |

The (-,-) ghost alone closes the gap completely; the (-,+) ghost recovers
only 1 of the 3. So the answer to "which 3 dimensions need the level-11
ghost" is: the three that only the degree-7 channel L+-·L-- can reach, and
L-- is the (-,-) line, which is entirely old. **The oldspace is
irreplaceable for one specific reason: the newforms have no (-,-) sector.**

## CI-A — the anomalous quadric: **predictions hit, interpretation refuted**

A1-A4 all confirmed. Q* lies in the (+,+) class (no y1·y_k terms), the
y1^2 coefficient is nonzero, and the quadric is the weight-4 identity

  a·e1^2 + q(e9, e10, e11, e12) = 0,   residual **exactly 0**,

with a = -2439613813. A4 confirmed too: under the grading the "excess 1"
is 11 products into a 12-dimensional target with rank 10, and the
class-(-,+) part has full rank 4.

**A5: negative.** The bilinear form of q is *not* Hecke self-adjoint —
T^t G - G T is nonzero. So q is not a trace form Tr_{K2/Q}(c x^2) and the
identity is not Hecke-canonical. G has rank 4 and is indefinite (diagonal
carries both signs).

**The decoy refutes the natural reading.** The tempting statement is
"e1^2 lies in Sym^2 of the f2 orbit" — a codimension-2 coincidence about
the 4-dimensional Galois orbit. It is not about f2. Replacing f2 by 20
random rational 4-dimensional subspaces of the 5-dimensional (-,+) sector:
**20 of 20 reproduce the relation.** The reason is in the table above:
Sym^2 of the whole (-,+) sector has rank exactly 10 = h^0(L-+^2), not 12,
and e1^2 lies in that 10-dimensional H^0(L-+^2) inside the 12-dimensional
H^0(M_++). Any 4-dimensional subspace spans the same 10-space, so any
choice produces the relation. The single true content is

  **e1^2 lies in H^0(L-+^2) subset H^0(M_++)** — equivalently, the
  degree-2 fixed divisor separating those two spaces on E is contained in
  the divisor of e1^2.

Since e1 = pi^*(omega_E) vanishes exactly on the ramification divisor and
the degree-2 divisor is ramification-supported, this is plausibly forced.
**Not established.** Registered open.

## An observation, flagged not claimed (AG-D5)

a = -2439613813 = **7^2 · 13 · 1957^2**, and 1957 = 19 · 103 is exactly
disc(K_f2), the discriminant of the quartic Hecke field of the f2 orbit.

a survives the invariance test: recomputed under six random unimodular
changes of the f2 lattice it does not move (up to sign), so it is a
genuine invariant of the integral structure and not basis noise.

That said, **no mechanism is established**, and disc^2 is exactly the kind
of factor that appears generically in Gram determinants of forms on rank-1
modules over an order. Per AG-D5 this is a proximity observation, not
evidence. Falsification route: compute the same invariant at other levels
N whose X0(N) has a 1-dimensional (+,+) sector and a Galois orbit filling
another sector. One instance is not a pattern; if disc^2 persists across
levels it is structural, and if it does not the observation dies.

## Status

EXACT: all ranks, the deficiency vector, Q*, the residual, the ten bundle
tests, the invariance of a. Cert: none upgraded beyond what was computed
here. No physics reading is claimed and none follows.

## Open, registered

1. Is e1^2 in H^0(L-+^2) forced by ramification, or a genuine condition?
2. The disc^2 observation, with the cross-level falsification above.
3. Identify L-- as an explicit degree-1 bundle on E, i.e. name the point.
   That would make the whole descent picture fully explicit.
