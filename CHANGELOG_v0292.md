# mtft v0.29.2 — surface.arithspin: arithmetic spin structures and the three-family flux (2026-09-13)

New GP-free module `surface.arithspin` (ARITH-SPIN-01 / GAL-01), all results EXACT unless noted:
- eta-quotient modular units on Gamma_0(N) (Newman conditions, Ligozat cusp orders, vectorised), the
  cuspidal relation lattice (Hermite basis) and cuspidal class group (Smith form): for X0(143)
  Z/420 x Z/10, order 4200; orders of (0)-(oo), (1/11)-(oo), (1/13)-(oo) = 420, 60, 70.
- K ~ 6 D_cusps; the four cuspidal theta characteristics (cuspidal 2-torsion (Z/2)^2), all AL-invariant;
  h^0 by Serre duality from AL-signed cusp orders of the frozen weight-2 eigenbasis: O(3D) has h^0 = 5 (odd;
  matches KK-06's cuspidal spin Dirac (5,5)), O(6(0)+6(1/11)) has h^0 = 2 (even); control h^0(K) = 13.
- the unit u = eta(13t)eta(143t)/(eta(t)eta(11t)) spanning H^0(S0) with the constant; u∘W13 = -(1/13)/u,
  u∘W11 = -u (numerical to 1e-8 at interior points; forced by divisors up to the constants).
- W13-fixed CM points (discriminant -52), u(P) = ±i/sqrt 13 with two of each sign, hence the purity theorem:
  L3 = S0 (x) O(P1+P2+P3) has exactly three zero modes (h^1 = 0) for every triple — a W13-symmetric,
  W11-breaking, Q(i, sqrt 13)-rational three-family flux; the odd exceptional spin O(3D) is impure.
Context recorded: the W11 parity theorem (v0.29.1) and the non-rationality of any AL-invariant odd-degree
flux (no cuspidal torsion twist makes O(3 oo) AL-invariant; CM fields have even degree).
Tests: `tests/test_surface_arithspin.py` (2, ~11 s).  Base: live 0.29.0 + unpushed 0.29.1.  Pin four-way.
