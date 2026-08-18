# Pre-registration — CI-D: locating P

Filed 2026-08-16 before computation.

L_-- = O_E(P) with 2P ~ D_W13 = pi(Fix(W13)). Four square roots exist;
which is P?

## Predictions

| # | prediction | reasoning |
|---|---|---|
| D1 | **P lies in E(Q)**, and since E = 143a1 has rank 1 with trivial torsion, **P = nG** for an integer n and a generator G | L_-- is a Q-rational degree-1 line bundle on a curve with a rational point, so its unique effective divisor is a rational point |
| D2 | modular degree of 143a1 is **4** | X0(143) -> X0(143)* has degree 4 and X0(143)* = E if 143a1 is the optimal curve |
| D3 | the four W13-fixed points map to exactly **2** distinct points of E | W11 acts freely and commutes with W13, so it pairs them |
| D4 | Q1 + Q2 = 2P = 2nG lies in E(Q) | Galois-stable set |

No prediction on the value of n.

## Method

Fixed points of W13 solve 143c tau^2 - 26a tau - b = 0 with a + d = 0 and
13ad - 11bc = 1, giving tau = (13a +- sqrt(-13))/(143c) — discriminant
-52, as the class-number route already said. Images computed by the
modular parametrisation z(tau) = sum_n a_n q^n / n, q = exp(2 pi i tau),
then ellztopoint, Manin constant assumed 1 (optimal curve).

## Decision rule

If Q1 + Q2 does not land numerically on a rational point of E, the
normalisation or the Manin constant is wrong. **Report and stop — do not
tune anything to force a fit.** If D3 fails (four distinct images), the
claim that phi factors through the full quotient is wrong and CI-C needs
revisiting.
