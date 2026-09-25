# mtft v0.22.0 — involutions, the oldspace abelian surface, the quadratic
# Hamiltonian layer, and Bergman harmonic channels

Promotes the verified v0.22-arc studies (Sol + Manus exploration + Claude
battery, 2026-08) into the wheel.  Four new `mtft.periods` modules; the
period surface now declares `__all__`.

## periods.involutions
- Frozen exact integer W_11, W_13 on the promoted 26D Hecke basis
  (`al_matrix`), certified at call time: W^2=I, [W_11,W_13]=0, det +1,
  W^T E_H W = E_H, [W,T_p]=0 for p=2,3,5,7 (EXACT), eigenspace dims
  (14,12,4) = 2 g(X/W), U_11 = -W_11, rank im[U_13,W_13] = 4 = oldspace.
- `al_signs` / `sector_census`: homology-side AL sign decode
  eps(W11)=(+,-,+), eps(W13)=(+,+,-) on (143a1,f2,f3); reconstructs the
  canonical differential census (1,6,5,1) and resolves the 5-sector as
  f2(4) + old_plus(1).
- `route2_fixed_intersections`: fix(W_143) meets (ell,ghost,q4,q6) in
  (2,2,0,0) — Route 2 on homology.
- `oldspace_projector`: P_old = (U13^2-I)(I-2U13)/90, EXACT idempotent.
- `star_symplectic` / `star_charge_orbit`: star transported to symplectic
  charge coordinates is integral, star^2=I, and ANTI-SYMPLECTIC
  (star^T J star = -J, EXACT), so charge energy is star-invariant; the
  minimal-shell degeneracy E(n3)=E(n5-n6)=0.881330420747955 is
  symmetry-protected (one 4-element orbit).

## periods.oldtorus
Call-time kernel-route re-derivation of the oldspace abelian surface:
Smith(E|L_old)=(2,2,18,18) (type (2,18)); J_arith=(U13-2I)/3 with
[L9:L_old]=9 and mod-3 rank 2 (9=3^2 structural); E/2 on L9 integral
unimodular alternating (principal); J_arith does NOT preserve E/2;
charpoly((J_arith W13)|L9)=(x^2-4x-1)^2; entropy 2 log(2+sqrt5).

## periods.hamiltonian
Hodge adjoint and Hermitian split (A,K); Hamiltonian J-split (A_+,A_-);
`pairing_stability` rho (width 0.1234286299, distance 0.4248813827, both
< 1); 13 Williamson `symplectic_frequencies` of F=JA (degree: 3 x13
exactly); `oldspace_routing` honest negative (ghost carries only
12.3%/14.8% of pairing power); `hecke_block_routing` (width ~75.4%
intra-block, distance ~54.5% inter-block).  Degree is the EXACT null
control: V_degree = 3I over Q, hence A_- = 0 for any complex structure.

## periods.channels
Bergman bilinear coefficients B_{n,m} = a_n^dag (Im tau)^-1 a_m; channel
series C_k(y); `channel_density` = exact FFT-free rearrangement of
`bergman_density` (gate residual ~7e-40 at full k-range); `mode_crossover`
solves |C_4|=|C_1| at y*/y0 = 2.302140221833918907.  The k=1 anomaly is
destructive interference (B_{2,1}=-0.3002 vs B_{3,2}=+2.3015), not a
vanishing leading coefficient.

## Gates, CLI, tests
`python -m mtft.periods verify` now runs 10 gates (six v0.21 + four new).
New CLI subcommands: involutions, oldtorus, hamiltonian [--potential],
crossover.  New test tier `tests/test_periods_v0220.py`.

## Fixes
- `mtft.periods.__init__` declares `__all__` (36 v0.21 names + v0.22).
- `cp_channel_report`: adds `commutator_star_odd_abs`; the relative ratio
  returns NaN below a norm floor instead of an ill-conditioned value
  (exposed by the exactly-zero degree commutator).
- `integral_lattice.saturate` raises TypeError on non-integer input
  (Fraction vectors from `rational_kernel` must be denominator-cleared).
- `tests/test_x0143_particle_box.py` skips cleanly without SciPy
  (importorskip) under a core-only install.

## Open register
Legend registration for the periods surface; interaction-algebra study
(ker/im of V -> V_-); certified shortest-vector search on (Z^26, G);
Kimi audit of the AL decode, L9 chain, star orbit, and Bergman channels.
