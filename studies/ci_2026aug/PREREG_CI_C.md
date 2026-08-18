# Pre-registration — CI-C: naming the descent bundles

Filed 2026-08-16 before computation. Successor to PREREG_CI_AB.md.

Goal: make the descent picture on E = X0(143)* explicit, and settle the
CI-A open item (is e1^2 in H^0(L_-+^2) forced, or a coincidence?).

## Derivation done in advance (all from data already certified)

Riemann-Hurwitz on X -> X/sigma, degree 2, using the genera already
computed (7, 6, 2 for W11, W13, W143) and g(X) = 13:

| sigma | g(X/sigma) | #Fix(sigma) = 24 - 2(2g-2) |
|---|---|---|
| W11 | 7 | **0** |
| W13 | 6 | **4** |
| W143 | 2 | **20** |

Total 24 = deg K_X, and inertia on a smooth curve is cyclic, so all
ramification of pi: X -> E is simple. Branch degrees on E are half the
fixed-point counts: (d_W11, d_W13, d_W143) = **(0, 2, 10)**.

| # | prediction |
|---|---|
| C1 | **W11 acts freely on X0(143)** — no fixed points |
| C2 | #Fix(W13) = 4, #Fix(W143) = 20, checked independently by the class-number formula: #Fix(W_Q) = h(-4Q)·prod_{p\|M}(1+(-4Q\|p)) + [Q=3 mod 4] h(-Q)·prod(1+(-Q\|p)) |
| C3 | this forces h(-572) = 13 (= the genus), given h(-143) = 7 |
| C4 | bidouble relations 2L_chi = D_j + D_k reproduce deg L = **(0, 6, 5, 1)** — an independent route to the degree vector, which was previously read off h^0 of the sectors |
| C5 | the four cusps form a **single free orbit** under <W11,W13>, giving E a distinguished rational point O |
| C6 | the (+,+) form e1 is the newform of **143a1**, so E is that elliptic curve |
| C7 | the (-,-) and (-,+) forms are combinations of g(tau), g(13 tau) with g = **11a1** = eta(tau)^2 eta(11 tau)^2 |

## The CI-A resolution, predicted

Vanishing rule at a fixed point of an involution sigma: a form in sector
chi vanishes there iff chi(sigma) = +1 (locally sigma: z -> -z, and
sigma^*omega = chi(sigma) omega forces f odd exactly when chi = +1).

e1 is (+,+), so chi(W13) = +1, so **e1 vanishes at all 4 fixed points of
W13**. The divisor separating H^0(L_-+^2) from H^0(L_+-^2) is exactly
D_W13 = pi(Fix(W13)), degree 2. Order 2 upstairs pushes to order 1 down.

| # | prediction |
|---|---|
| C8 | **CI-A is FORCED**, not a coincidence |
| C9 | companion, same mechanism: e1·e13 must lie in the span of the 30 f3 x (-,+) products, because both e1 and e13 have chi(W143) = +1 and so vanish on Fix(W143) |
| C10 | div(e1) = Fix(W13) + Fix(W143) exactly (4 + 20 = 24 = deg K_X) |

## Naming L_--

C11: L_-- = O_E(P) is the line bundle defining the double cover
C2 = X0(143)/W143 -> E, branched over D_W13 = the image of the four CM
points of discriminant -52 fixed by W13. C2 has genus 2 and
Jac(C2) ~ 143a1 x 11a1, since its two differentials are the (+,+) and
(-,-) sectors and the latter is level-11 oldspace. 2P ~ D_W13.

Open, no prediction: **which** of the four square roots P is; that needs
the branch points located in E(Qbar).

## Decision rule

If C1 or C2 fails, the genera from the AL traces are wrong and everything
from the last two sessions is void. If C8's mechanism is right, C9 must
already hold in the data computed last session — that is a check I can
make without new freedom.
