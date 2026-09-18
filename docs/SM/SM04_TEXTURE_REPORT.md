# SM-04 — W13 texture of the Yukawa tensor, the degeneracy theorem, and the first falsifiable relation (2026-09-13)

## Exact structure
W13 acts on the three-family space by T(phi) = u (phi o W13) (T^2 = −1/13; normalised, an involution) and on
the Higgs space by the weight-6 slash.  Grades: sections (−, +, +); Higgs directions 8 even, 10 odd.  The
tensor obeys the selection rule  Y_ij^k = 0 unless grade_i grade_j grade_k = −1  (forbidden entries 7e-8,
allowed O(1)).  Consequences:
- a W13-EVEN Higgs VEV couples odd to even only: the 3 x 3 mass matrix has rank 2 — one massless family;
- a W13-ODD Higgs VEV couples odd–odd and even–even: rank 3, block-diagonal (1 + 2) in the graded basis.
The Z2 surviving the parity theorem is a flavour symmetry with a definite texture.  Which VEV direction is
realised is a potential question (open).

## Degeneracy theorem (untwisted M1) and its test
With L_b = L_a and L_d = L_c = O, the up (Q u^c H_u), down (Q d^c H_d) and charged-lepton (L e^c H_d) couplings
all involve the same bundles S0(SigmaP), S0(SigmaP) and O(−2 SigmaP): the same tensor.  Since down and
lepton share H_d, untwisted M1 PREDICTS m_d : m_s : m_b = m_e : m_mu : m_tau at the compactification scale.
This is the SU(5)-type relation known to fail by factors of ~3 in the light generations (Georgi–Jarlskog).
So untwisted M1 is falsified at leading order — the program's first genuine contact with data — unless the
twist that separates L_b from L_a supplies the missing structure: a cuspidal 2-torsion twist is a quadratic
character of Gamma_0(143), and the twisted down/lepton wavefunctions are modular forms with character.
That is the next computation (twisted Riemann–Roch spaces via PARI character spaces, Gamma_0-only reduction).

## On the older arithmetic values
Nothing computed before this program (alpha^-1, Koide, weight-derived ratios) enters or is explained by the
construction.  What the construction can deliver are dimensionless ratios once the sections are normalised
(compact metric) and a Higgs direction is fixed; the inputs then are R, the 6D gauge coupling and that
direction — not "the Planck mass and the electron".  The first comparison with data is the relation above.
Tool: `rrspace.w13_grading_and_texture`; test in `test_surface_rrspace.py`.
