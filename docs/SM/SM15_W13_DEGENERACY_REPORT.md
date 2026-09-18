# SM-15 — the W13-symmetric vacuum: degeneracy theorem and epsilon admixture (2026-09-14)

## Setting (post CC-26)
With the correct kinetic normalisation the leading-order mass matrices of M1/M2 are anarchic (all ratios O(1)).  The
one exact structure available is the W13 grading of the untwisted up sector: family grades (+, -, -), Higgs space
8 even + 10 odd, selection rule grade_i grade_j grade_k = -1.  The HYM metric is W13-invariant, so the Gram
matrices are block-diagonal in the grading (off-grade entries 5e-3 family / 5e-4 Higgs, mesh level) and the
normalisation can be done grade by grade (`hym.w13_graded_normalisation`).

## Degeneracy theorem (EXACT; verified in the corrected normalisation)
For every W13-even Higgs direction the mass matrix has the form [[0,a,b],[a,0,0],[b,0,0]] in the graded basis, hence
singular values (r, r, 0): one massless family and an exactly degenerate heavy pair.  Numerically (1, 1, 1e-9) on
random even directions.  Odd directions give rank 3 with O(1) ratios.

## Epsilon admixture (DIAGNOSTIC, 300 random pairs per eps)   v = v_even + eps v_odd
   eps      1e-4      1e-3      1e-2      1e-1      0.3       1.0
   m1/m3    7.2e-5    7.2e-4    7.1e-3    6.9e-2    0.18      0.25
   m2/m3    0.9999    0.9991    0.991     0.906     0.76      0.59
The light family scales linearly, m1/m3 ~ 0.7 eps; the heavy pair stays degenerate to O(eps).

## Consequence (EXACT)
A W13-symmetric (or nearly symmetric) Higgs vacuum implies m_c = m_t (and, in M2, degenerate heavy Dirac neutrinos).
The observed up spectrum (m_c/m_t ~ 3e-3, m_u/m_t ~ 6e-6) therefore requires O(1) breaking of W13 in the Higgs
direction, at which point the spectrum is anarchic again.  The W13 grading cannot be the origin of the quark
hierarchy in this class.  The observed hierarchy needs a second small parameter acting inside the heavy pair —
nothing in the present construction supplies one.

## Corrected CKM (replaces SM-14)
Random direction pairs in the corrected M2: |V_us|, |V_cb|, |V_ub| all 0.54 at the median (5–95%: 0.15–0.88):
Haar-like anarchy.  The observed small mixing is not a typicality statement of this class at leading order.

## Lessons kept in the register
1. CC-26: a normalisation convention must be tested by an identity (A N A^dag = I) and by basis invariance, never
   by whether a desired hierarchy appears.
2. The exact textures (rank 2, degeneracy) survived the correction because they are algebra; the hierarchies did
   not because they were numerics without an invariance gate.  Exact and diagnostic must be kept apart in the text
   as well as in the code.
3. A symmetric vacuum is a hypothesis with consequences (here: m_c = m_t); test the consequence before building on
   the hypothesis.
Tools: `hym.w13_graded_normalisation`, `hym.w13_epsilon_scan`; test `test_surface_w13_texture.py` (degeneracy theorem,
tolerance 1e-6 over a measured 1e-9).
