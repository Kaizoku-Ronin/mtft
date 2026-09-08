# Correction to the previous Higgs audit wording

7 September 2026. Supersedes the interpretation of **SM-AUDIT-01** in the earlier Standard Model assessment and electron–photon report. Their historical numerical results and source variants are preserved.

I previously described `HIGGS.lambda_quartic = 0.5179175466708401` as an implementation defect because it was four times the coupling used by `HosotaniMTFT.higgs_self_coupling()`. That description assumed that both API names denoted the same potential convention. The factor-four comparison was correct; calling the larger number intrinsically wrong was too strong.

MTFT 0.26.2 makes the conventions explicit:

- The existing `lambda_quartic` is retained, with `m_H² = lambda_quartic * v² / 2`, corresponding to `V = (lambda_quartic/4)(H†H-v²/2)²`.
- The new `lambda_quartic_pdg = lambda_quartic/4 = 0.12947938666771003` uses `m_H² = 2*lambda_quartic_pdg*v²`. It agrees with the default Hosotani value.

These are two names and normalizations for the same potential coefficient. Our earlier local source variant chose the second convention for the original public accessor; the release instead preserves the original API and supplies a conversion accessor. Use the release's explicit accessor for comparisons. The earlier local patch should not be applied to 0.26.2.

The separate **SM-AUDIT-02** determinant defect was a mathematical implementation error. Release **CC-21** adopts the full-determinant row correction and retains `random_su_n_defective_v0261` for reproducing old draws.

The electron–photon charge, covariance, Maxwell, Coulomb, and Pauli experiments did not depend on the numerical value of this Higgs quartic. Their conclusions are unaffected by this convention clarification.
