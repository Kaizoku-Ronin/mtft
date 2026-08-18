# CI-C: the descent made explicit, and CI-A closed

**Session 2026-08-16. Pre-registered in `PREREG_CI_C.md`. All predictions
confirmed except C3, which failed for an instructive reason and produced
a corpus correction.**

---

## 1. The ramification data

Riemann-Hurwitz on each degree-2 quotient, from genera already certified:

| involution | g(X/sigma) | #Fix (RH) | #Fix (class-number formula) |
|---|---|---|---|
| W11 | 7 | **0** | **0** |
| W13 | 6 | **4** | **4** |
| W143 | 2 | **20** | **20** |

Two independent routes, agreeing. The second uses
#Fix(W_Q) = h(-4Q)·prod_{p|M}(1+(-4Q\|p)) + [Q=3 mod 4] h(-Q)·prod(1+(-Q\|p)):

- **W11**: h(-44) = 3, h(-11) = 1, but kron(-44,13) = kron(-11,13) = -1,
  so both Euler factors vanish. **W11 acts freely on X0(143)** — no fixed
  points at all. Total 24 = deg K_X is carried entirely by W13 and W143.
- **W13**: h(-52) = 2, kron(-52,11) = +1, giving 2·2 = 4. The four fixed
  points are the CM points of discriminant **-52**.
- **W143**: h(-572) + h(-143) = 10 + 10 = 20.

Branch degrees on E are half the fixed-point counts:
**(d_W11, d_W13, d_W143) = (0, 2, 10)**.

## 2. The degree vector, re-derived independently (C4)

The bidouble relations 2L_chi = D_j + D_k, using only the branch degrees:

| chi | kernel | 2 deg L | deg L | sector route |
|---|---|---|---|---|
| (+,-) | <W11> | d_W13 + d_W143 = 12 | **6** | 6 |
| (-,+) | <W13> | d_W11 + d_W143 = 10 | **5** | 5 |
| (-,-) | <W143> | d_W11 + d_W13 = 2 | **1** | 1 |

The degree vector **(0, 6, 5, 1)** now has two derivations sharing no
step: one from h^0 of the Atkin-Lehner sectors, one from ramification.
That closes the E2 on the descent structure.

## 3. Naming the curve and the bundle

**E is the elliptic curve 143a1.** The (+,+) form is literally the 143a1
newform: e1 = **72 · f_143a1**, exact on all 30 coefficients checked.
143a1 has conductor 143, **rank 1**, trivial torsion, discriminant -1859.

**E has a distinguished rational point.** The four cusps, indexed by
d | 143 with width 143/d, are permuted by W11 as (1 11)(13 143) and by
W13 as (1 13)(11 143) — a **single free orbit**. So all four cusps map to
one point O in E(Q), the natural origin. (Consistent with W11 free.)

**L_-- named, three ways:**

1. **Its section, explicitly.** With g = 11a1 = eta(tau)^2 eta(11 tau)^2,
   the two oldspace eigendirections are

   e8  = g(q) + 13·g(q^13)   sector (-,+)
   **e13 = g(q) - 13·g(q^13)   sector (-,-)**

   so the unique section of L_-- is g(tau) - 13 g(13 tau). Both sit in
   W11 = -1 sectors, matching 11a1 having rank 0, root number +1, hence
   W11-eigenvalue -1 — an independent consistency check.

2. **As a covering datum.** L_-- is the line bundle defining the double
   cover C2 = X0(143)/W143 -> E, branched over the degree-2 divisor
   D_W13 = pi(Fix(W13)) — the image of the four discriminant -52 CM
   points. 2L_-- = D_W13, so L_-- = O_E(P) with 2P ~ D_W13.

3. **Through the genus-2 curve.** C2 has genus 2 and its two differentials
   are exactly the (+,+) and (-,-) sectors, so
   **Jac(X0(143)/W143) ~ 143a1 x 11a1**. C2 is bielliptic, covering 143a1
   via the newform line and 11a1 via the ghost line.

Still open, as pre-registered: **which** of the four square roots P is.
That needs the two branch points located in E(Qbar).

## 4. CI-A is closed: the relation is FORCED

Vanishing rule: at a fixed point of an involution sigma, a form in sector
chi vanishes iff chi(sigma) = +1 (locally sigma: z -> -z; chi(sigma) = +1
forces the coefficient function odd, hence a zero).

e1 is (+,+), so chi(W13) = +1, so **e1 vanishes at all four fixed points
of W13**. The divisor separating H^0(L_-+^2) from H^0(L_+-^2) is exactly
D_W13, degree 2. e1^2 vanishes to order 2 upstairs at each fixed point,
which pushes to order 1 downstairs through the ramification — exactly
what is needed. **So e1^2 lies in H^0(L_-+^2) necessarily, and the
"anomalous quadric" of CI-A is not an anomaly at all.**

This also explains the decoy result cleanly: nothing about f2 or about any
4-dimensional subspace enters the argument, which is why 20 of 20 random
subspaces reproduced it.

**Companion prediction C9, confirmed retroactively.** The same mechanism
says e1·e13 must be dependent on the 30 f3 x (-,+) products, because both
e1 and e13 have chi(W143) = +1 and so vanish on the 20 W143-fixed points.
Last session's table already shows it: class (-,-) has 31 monomials of
rank 11, and the 30 cross products alone have rank 11. The extra monomial
y1·y13 contributes nothing. Predicted after the fact was impossible — the
number was already on the page.

C10 also holds: div(e1) = Fix(W13) + Fix(W143), degree 4 + 20 = 24 =
deg K_X, so e1 = pi^*(omega_E) vanishes exactly on the ramification.

## 5. CC-09 (PROPOSED) — the class number of Q(sqrt(-143))

**C3 failed**, and tracing why produced a corpus error.

Arithmetica Generale **Pr 7.8.1** states that the class number of
Q(sqrt(-143)) is 7, that j(tau_0) is an algebraic integer of **degree 7**
over Q, and supports this with "the Minkowski bound is sqrt(143)/pi =
3.81, so we check primes 2 and 3" and "the discriminant -4 x 143 = -572
has 7 reduced forms".

**Correct values: h(-143) = 10 and h(-572) = 10.** The class group of
discriminant -143 is cyclic of order 10.

**Diagnosed mechanism.** The Minkowski bound for an imaginary quadratic
field is (2/pi)·sqrt(|d|), not sqrt(|d|)/pi — a factor of 2. The text used
3.81 where the correct bound is **7.61**, so it tested only p = 2, 3 and
missed the classes above 3.81, undercounting 10 as 7.

**Downstream.** The degree of j(tau_0) over Q in Pr 7.8.1 is **10**, not
7. Any claim keyed to a degree-7 minimal polynomial at the CM point needs
re-checking. Note the fixed-point count #Fix(W143) = h(-572) + h(-143) =
20 is unaffected, because it happens to use the corrected values and
still lands on the Riemann-Hurwitz answer — the arithmetic here is
self-checking.

**AG-D5 follow-through.** The leading coefficient of the CI-A quadric
factored as a = 7^2 · 13 · 1957^2, and 7 = h(-143) was an available
reading. That reading is now dead: h(-143) = 10. The disc^2 = 1957^2
factor stands as flagged; the 7^2 and 13 remain unexplained and are not
claimed.

## 6. Ledger

| item | status |
|---|---|
| W11 acts freely on X0(143) | Pr (two routes) |
| branch degrees (0, 2, 10) | EXACT |
| deg L = (0, 6, 5, 1) | Cert (E2: sectors and ramification) |
| E = X0(143)* = 143a1, rank 1 | Cert (q-expansion match, 30 coefficients) |
| cusps = one free AL orbit, origin O in E(Q) | EXACT |
| e13 = g(q) - 13 g(q^13), g = 11a1 | EXACT |
| Jac(X0(143)/W143) ~ 143a1 x 11a1 | Pr |
| CI-A forced | Pr (vanishing rule) |
| C3 (h(-572) = 13) | **FALSE** — honest miss, recorded |
| CC-09 h(-143) = 10 not 7 | PROPOSED |

## 7. Open

1. Locate the two branch points in E(Qbar) and pin P among the four
   square roots of D_W13.
2. CC-09 adoption, and an audit of anything in the corpus keyed to
   h(-143) = 7 or to a degree-7 CM minimal polynomial.
3. The disc^2 observation from CI-A, cross-level falsification as
   registered.
