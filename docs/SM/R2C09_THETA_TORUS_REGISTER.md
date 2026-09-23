# R2C-09 — the torus factor on 143a1; CC-32 (2026-09-19)

## CC-32 (my error in R2C-08)
The R2C-08 selection rule used Kuenneth parity alone.  The Yukawa is a (0,2)-form on S and Lambda^{0,2}(S) = Lambda^{0,1}(X) (x) Lambda^{0,1}(E),
so the BIDEGREES (q_X, q_E) of the two families and the (0,1) Higgs must sum to (1,1).  The R2C-08 "solutions" (e.g. Q (3,-1), u^c (-1,-3))
sum to (2,2) and vanish.  Corrected scan (`product_surface.scan_family_pairs_full`): exactly two mixed-origin solutions (mirror pairs):
   S1: Q (3, 1) from the curve [3 pure modes at the CM points], u^c (1, -3) from the torus [index -3, with 3 vector-like pairs from h^0(S0(P)) = 2],
       Higgs (-4, 2), bidegree (1,0), 32 modes (h^1(X,-4) = 16 x h^0(E,2) = 2), slope-free at A_X = 2 A_E.
   S2: Q (-3, 1) from the curve [3 modes via h^1(S0(-sum P)) = 3, Serre-dual to the family space], u^c (1, 3) from the torus [6 vs 3],
       Higgs (2, -4), bidegree (0,1), 4 modes (h^0(X, O(sum P - P)) = 1 x h^1(E, -4) = 4), slope-free at A_X = A_E/2.
   Every curve x curve pair fails (Higgs with a zero degree on one factor or no modes of the required bidegree).

## The torus factor (theta functions on 143a1) — `research.theta_torus`
143a1: c4 = 64, Delta = -1859 = -11 * 13^2, j = -262144/1859 (exact); tau = 1/2 + 1.0232745927...i on Re tau = 1/2.  Validation: the standard
theta bases are orthogonal with equal norms (off-diagonal 2e-26) and theta_1 theta_2 closes in the degree-3 span (residual 1e-25) — the
classical facts reproduced by the quadrature before any new number.
   S1: B = <theta3_j, theta1 theta2_alpha>, 3 x 2 (orthonormal bases): singular values (1.103, 0.859), rank 2; rows 2 and 3 equal (z -> -z).
   S2: B' = <theta4_alpha, theta1 theta3_j>, 4 x 3: singular values (1.124, 0.979, 0.843), rank 3.

## Structural consequence (EXACT)
In both solutions the u^c curve factor is H^0(S0(P)), h^0 = 2 (h^1 = 1 for any single CM point: a + b u vanishing at P is one condition).  The up
mass matrix M_{i,(a,j)} = sum_v A_{i a ...} B_{... j} therefore has rank <= 2 at leading order: the lightest up-type quark is massless at
leading order, from a Riemann–Roch number.  In S2 the curve factor A is the Serre pairing of the family space with its 2-dimensional subspace
of sections vanishing at two CM points (products sigma_a s with s the section of O(sum P - P)); the curve Higgs is a single mode, so m_c/m_t
depends only on the four torus Higgs directions and the fixed A — the next computation, and a far more constrained prediction than M1's
eighteen-dimensional anarchy.
OPEN: 8D anomalies; moduli stabilisation at A_X = 2 A_E or A_E/2; the three vector-like u^c pairs; the down and lepton triangles on the surface.


## Addendum (2026-09-23, v0.33.0): the structural consequence, qualified; second route to the torus factor
rank M <= min(3, 2 * kunneth_rank(v)).  In S2 (curve Higgs factor one-dimensional) rank <= 2 holds for every VEV; in S1 a Kuenneth-rank-two VEV
reaches rank 3, so "m_u = 0 at leading order" holds there for Kuenneth-rank-one VEVs only (`theta_torus.up_mass_rank_bound`, `up_mass_rank_examples`).
The torus factors are reproduced without quadrature by the theta multiplication formula, B_{j,b} = c_j(b) (2k/((k+1) Im tau))^{1/4}, to 1e-25
(`theta_torus.torus_factor_closed_form`); the singular values are a numerical certificate on two routes (tau is transcendental).  h^0(O(sum P - P)) = 1
for P in {P1, P2, P3} is proved (gonality >= 4); for P = P4 the S2 Higgs factor has no sections.  `V0330_COMPENDIUM_NOTES.md` items 4–6.
