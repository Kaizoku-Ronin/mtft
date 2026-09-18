# SM-07 — the cubic twist of the singlet stack: the lepton tensor and the lifted degeneracy (2026-09-13)

## Construction (validated by the same product test as SM-06)
The cuspidal group's 3-torsion (140((0)−(oo)) has order 3) corresponds to the cubic characters of (Z/13)^*,
which are even.  Conrey labels mod 143: chi_133 (cubic), chi_100 (its conjugate), chi_56 = chi_13 chi_133-bar
(sextic), chi_23.  Character convention fixed by a discriminating test — direct q-series evaluation versus the
reducer through a Gamma_0(143) element with d = 2 — residual 1e-14 for chi, 1.7 for chi-bar.
With L_d = L_c (x) t_cubic:  L in the chi_133-twisted family space, e^c in the chi_56-twisted one, H_d in the
chi_13 Higgs space of SM-06.  Cusp conditions via the AL matrices (W_13, W_143 into the conjugate spaces):
both family spaces certify at 37 of 40 conditions (dim 3, pure).  The lepton tensor L_i e^c_j H_d^k lands in
f_K x (chi_13 Higgs space) with residual 8e-9 — the character identity chi_133 * chi_56 = chi_13 checked by
the E2.  Generic rank 3.  Frozen: `_data/x0143_cubic_spaces.npz` (654 KB).  Tools: `rrspace.load_cubic`,
`conrey_chi`, `cubic_values`, `cubic_family_space`, `lepton_yukawa_M1`, `family_equivalence_residual`.

## Result: the down–lepton degeneracy is lifted
The transposition theorem (SM-06) forces m_e : m_mu : m_tau = m_d : m_s : m_b for every quadratic twist.  With
the cubic twist the lepton and down tensors are built from different family spaces and are NOT related by
any GL(3) x GL(3) change of family bases: alternating least-squares fit residual 2e-2 to Y_down and 2e-2 to
its family transpose, against 6e-16 for the self-fit control.  So the equality of ratios is no longer forced.
Its size — whether the arithmetic produces the Georgi–Jarlskog factors — needs the kinetic normalisations
(the Hermitian–Yang–Mills metrics on the flux bundles, i.e. Green's functions of the hyperbolic Laplacian at
the CM points, computable with the SPEC finite-element stiffness matrix) and a Higgs direction.

## Vacuum structure (exact, recorded)
The extension classes of the recombined bundles are flat directions of the leading-order potential: every
non-split extension of the given degrees is stable and has the same Atiyah–Bott energy, so the Higgs
direction is a modulus (P^17 for each Higgs sector) at this order.  Consequently the program's predictions
are relations that hold for ALL Higgs directions — such as the 2-torsion transposition theorem — while
individual mass ratios are functions on the moduli space until something lifts it.

## Status
EXACT: existence and purity of the cubic-twisted sectors, the character identity, the non-equivalence.
DIAGNOSTIC (1e-8): the tensor entries.  Open: normalisations (HYM metrics), the moduli problem, the up-sector
counterpart (a twist between the a and c stacks would split u from d further), neutrino sector (nu^c in the
chi_100-twisted space couples to L and H_u).
