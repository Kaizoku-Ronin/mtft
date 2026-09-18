# SM-02 — step 5 machinery: Yukawa sections with CM poles (2026-09-13)

## Two exact facts that make step 5 constructive
1. div(f1 dz) = all 24 Atkin–Lehner fixed points.  The (+,+) differential (the 143a1 newform) is even under
   W13 and W143; a weight-2 form with eigenvalue +1 must vanish at each fixed point of that involution
   (the automorphy factor there is −1), so f1 vanishes at the 4 W13-fixed and the 20 W143-fixed CM points,
   which is its whole degree.  Verified numerically: 4e-10 (W13 points), 1.6e-4 (W143 points, limited by the
   130-term expansions at height 0.01); W13/W143-odd sectors do not vanish.  The same argument shows every
   W13-even differential vanishes at the W13 points and every W143-even one at the W143 points.
   (Weight-1 CM forms do NOT vanish at the fixed points — the naive weight-1 argument fails; recorded.)
2. Cuspidal degree-3 fluxes are impure with S0 (one vector-like pair, because either the constant or u
   survives the twist), while CM-point fluxes are pure: purity forces the CM points.

## Construction of H^0(S0 (x) O(P1+P2+P3)), dim 3 (DIAGNOSTIC precision, exact structure)
Sections phi = G / (f1 · f_K) with f_K = g1^2 w (cuspidal divisor 2A, w = eta quotient (2,0,−2,0)) and G in the
36-dimensional canonical-ring space H^0(K^2) (quadratic products of the 13 eigenforms) satisfying:
ord_q(G) >= 8 at cusps 0 and 1/11 (12 exact conditions via AL signs), G = 0 at the 20 W143 points
(evaluated after height-maximising reduction; W11-partners generated explicitly), G(P4) = 0.
36 − 33 = 3.  The sections are 1, u and a third with simple poles at exactly the two points of the triple that
share a sign of u — precisely as the purity theorem's bookkeeping predicts (L(A + same-sign pair) is already
3-dimensional).  Data: /home/claude/yukawa_sections.npz (basis indices, nullspace).
Precision caveat: 130 coefficients give 1e-4 at the lowest reduced heights; rank certification of the 33
conditions needs longer expansions (a PARI dump to q^400) — the STRUCTURE (which points carry poles) is the
check that passed.

## Remaining for the up-type Yukawa of M1
Q_i u^c_j H_k = the multiplication H^0(S0(SigmaP)) (x) H^0(S0(SigmaP)) -> H^0(K(2 SigmaP)) (dim 18), paired with
the Higgs wavefunctions in H^1(O(−2 SigmaP)) = H^0(K(2 SigmaP))^*.  Basis of H^0(K(2 SigmaP)): G'/f1^2 with G' in
H^0(K^3) (60-dim, cubic products) vanishing to order 2 at the 20 W143 points and at P4, cusp orders fixed
— 60 − 40 − 2 = 18.  Then the 3 x 3 x 18 tensor is the change of basis of phi_i phi_j f1^2 f_K^2 in that space:
the same solve as `yukawa.yukawa_tensor`, with point conditions.  Normalisation by the compact metric
(uniformising factor) is the last ingredient.  Tool to ship: `surface.rrspace` (Riemann–Roch spaces from
canonical-ring elements with CM-point conditions), after a q^400 basis dump makes the point conditions
certifiable.
