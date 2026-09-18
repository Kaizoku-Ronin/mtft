# SM-06 — the chi_13 twist: down-type tensor, and the down–lepton transposition theorem (2026-09-13)

## Fixed and validated
The twisted sector failed its product test in SM-05 (residual 8e-2).  Cause: the cusp-order conditions taken
from PARI's `mfslashexpansion` for character forms were inconsistent.  Second route: the expansions of a
chi_13-form at the cusps 1/13, 1/11 and 0 are, up to constants, the q-expansions of F|W_11, F|W_13, F|W_143,
and PARI's `mfatkininit` gives those operators as RATIONAL matrices on S_k(143, chi_13) (W_11^2 = −1,
W_13^2 = W_143^2 = +1, columns = images).  With the conditions rebuilt from these matrices the twisted family
space (40 − 37 = 3) and twisted Higgs space (68 − 50 = 18) satisfy the product relation
G_i G^t_j = f_K sum_k Y G'^t_k to 4e-9 (the transposed orientation gives 0.4, the slash route 8e-2).  The same
E2 that rejected the first construction accepts this one.  Frozen: `_data/x0143_chi13_spaces.npz`
(S_2, S_4, S_6 to q^1200 and the six AL matrices, 473 KB).  Tool: `rrspace.twisted_spaces_chi13`,
`rrspace.down_yukawa_M1`, `rrspace.reduce_gamma0`, `rrspace.chi13_values`.

## Result
Model M1 with L_b = L_a (x) t3 (chi_13 twist on the b stack): the down-type tensor Q_i d^c_j H_d^k is
3 x 3 x 18 with residual 4e-9, generic rank 3.  It is a different tensor from the up-type one (untwisted x
untwisted vs untwisted x twisted, different Higgs spaces): the up/down degeneracy of untwisted M1 is lifted
structurally.  The twisted sector is Brill–Noether pure (h^0(S0 (x) t3) = 0).

## Transposition theorem (EXACT; the down–lepton degeneracy survives 2-torsion twists)
In the five-stack class, Q in (c, L-bar), d^c in (b, c-bar), L in (d, L-bar), e^c in (b, d-bar), with the shared
Higgs H_d in (L, b-bar).  Write L_d = L_c (x) t.  Then L lives in V_L (x) t and e^c in V_b (x) t^{-1}, where V_L, V_b
are the Q and d^c spaces, and both couplings land in the same Higgs target.  If t is 2-torsion (t = t^{-1}),
the lepton tensor is obtained from the down tensor by exchanging the roles of the two family spaces:
   t trivial:   Y_lep = Y_down;      t = t3:   Y_lep(i,j,k) = Y_down(j,i,k)   (with L_b = L_a (x) t3),
so M_e = M_d or M_e = M_d^T (kinetic normalisations transpose along), and the singular values — the mass
ratios m_e : m_mu : m_tau and m_d : m_s : m_b — coincide exactly at the compactification scale.  Hence the
Georgi–Jarlskog failure of untwisted M1 is NOT cured by any quadratic twist.  Lifting it requires
L_d L_c^{-1} of order >= 3: a cubic (or higher) cuspidal torsion class — the order-3 cuspidal classes correspond
to the cubic characters of (Z/13)^*, which are even, so their twisted spaces are computable by the same
route — or a continuous Wilson line on the singlet stack.  Purity then has to be re-checked for the order-3
twists (the 2-torsion twists t1, t2 are impure because of S1 = O(3D)).

## Status
EXACT: the theorem, the identification of the twists, the dimension counts.  DIAGNOSTIC (1e-9): the tensor
entries.  Not yet: normalisations (compact metric), the cubic-character sector, the Higgs directions.
The first falsifiable relation of the program has therefore been sharpened: in this class, quarks and leptons
of the down sector can only be split by a non-quadratic twist of the singlet stack.
