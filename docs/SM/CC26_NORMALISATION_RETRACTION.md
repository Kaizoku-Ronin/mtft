# CC-26 — transposed family contraction in the HYM normalisation: retraction of SM-08…14 mass-ratio claims
(2026-09-14; found by Astra, Atlas v0.31.1 Study Edition, audit item V0311-C01; verified and corrected here)

## The error
`hym.up_sector_mass_ratios` (v0.30.6–0.31.1) and every session computation behind SM-08…14 normalised the family
indices with `einsum('ai,bj,abk->ijk', L_A^{-1}, L_B^{-1}, Y)`, i.e. with (L^{-1})^T instead of L^{-1}.  For N = L L^dag
the orthonormal sections are A phi with A = L^{-1} (A N A^dag = I); the transpose is not a kinetic normalisation.
Astra's two-by-two rational example detects it by hand (A^T N A has diagonal 37/9, 13/9).  Verified on the saved
Grams: A N A^dag = I holds for L^{-1} and fails for (L^{-1})^T; the corrected normalisation is invariant under
re-expression of the family bases (5.7e-11 on M1; ~1e-8 on M2 through the regularised Higgs Gram) and the transposed
one is not.

## What changes (corrected numbers, same meshes, same directions)
   M1 up (SM-08/09):   released (1.4e-8, 1.7e-4)  ->  correct  m1/m3 = 0.26 (0.07–0.53), m2/m3 = 0.62 (0.34–0.87)
   M2 up (SM-13):      released (1e-5, 7e-3)      ->  correct  0.19 (0.05–0.39), 0.56 (0.34–0.81)
   M2 down (SM-12/13): released (2.8e-3, 7.2e-2)  ->  correct  0.25 (0.07–0.53), 0.61 (0.35–0.87)
Astra's independent replay of M1 (128 directions): medians (0.236, 0.625) — agreement.

## Retracted
- All mass-ratio distributions, hierarchies and their comparison with observation in SM-08, 09, 10, 12, 13, 14.
- The "hierarchy classes" (strong x strong / mild x strong / mild x mild) and "twist character -> hierarchy class".
- The "spin-structure invariant 7.8e-5" as a physical statement: it is the Gram spectrum in the SVD-nullspace
  coefficient basis — an invariant of that coefficient metric, not of the bundle with its HYM metric.  The numbers
  themselves (and their two-mesh agreement) stand only as properties of that basis.
- The CKM distributions and the joint fit of SM-14; the "up hierarchy reproduced" headline of SM-13.
- The frozen `x0143_m2_tensors_h02.npz` of v0.31.0/0.31.1 (regenerated with the correct contraction; the raw tensors
  and Gram matrices are now shipped so the normalisation can be re-audited).

## Unaffected (exact, algebraic)
Yukawa tensors and their product tests; the W13 selection rules and textures (rank-2 for even VEVs, exact); the
chi_13, cubic and sextic twisted spaces and the character-neutrality algebra; the transposition theorem; the
down = lepton identity of M2; nu Dirac = up in M2; the four-stack no-go, W11 parity, purity; sin^2 theta_W = 3/8;
the Liouville metric, Gauss–Bonnet gate and Green's functions; the Gram matrices themselves.

## What the corrected picture says
At leading order the kinetic normalisation does NOT generate fermion mass hierarchies in this class: all three
families come out O(1).  Any hierarchy is reachable by choosing a Higgs direction near the rank-deficient loci
(det M(v) = 0 is a cubic hypersurface; minimisation reaches m1/m3 -> 0 and m2/m3 ~ 1e-14 in both M1 and M2), so
masses are unconstrained until a mechanism selects the direction.  The one structured mechanism in hand is exact:
the W13 texture gives rank 2 at even VEVs, so a small odd admixture epsilon gives m1/m3 ~ epsilon — a
Froggatt–Nielsen-like origin of hierarchy from an exact selection rule plus one small symmetry-breaking parameter.
That is the next experiment, and it is a different claim from the retracted one.

## Process
The test that would have caught this — normalisation invariant under re-expression of the family basis — now ships
(`hym.normalisation_is_basis_invariant`, `hym.normalise_yukawa` as the single implementation).  The gate tolerance
follows the v0.31.1 margin rule (measured 1e-8, tolerance 1e-6).  Credit: Astra (finding), Kimi (CI margin rule).
