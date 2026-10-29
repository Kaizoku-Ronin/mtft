# Unreleased — after v1.0.0

Everything after v0.34.0 shipped in v1.0.0 ([release notes](CHANGELOG_v100.md)). Since then:

## CC-71 (2026-10-29): the Tano metric of the diatomic gases — S°(X₂) over the mass of X₂

- `mtft.tano_metric` divided the standard entropy of H₂, N₂ and O₂, which is tabulated per mole of molecules, by the molar mass of one
  atom. That doubled their Tano values: H 129.64 → 64.82, N 13.68 → 6.84, O 12.82 → 6.41 J g⁻¹ K⁻¹. `Element` gains `atoms` (atoms per
  standard-state formula unit; default 1, and 2 for H, N, O), and `tano_metric` now divides by `atoms × molar_mass`. The values before the
  correction stay exposed as `tano_metric.TANO_METRIC_RETRACTED_CC71` (append-only). `seebeck_diagnostic` computes its own T_a and now
  uses the same normalisation; nothing changes for the metals it is meant for.
- Second route: the three gases counted from first principles reproduce the stored S° per mole of molecules: H₂ 130.682 against 130.68
  (ortho and para hydrogen; counting every J once instead gives 5.71 more, close to R ln 2), N₂ 191.564 against 191.61, O₂ 205.128 against
  205.14 (only odd rotational levels exist, because ¹⁶O has spin 0; R ln 3 for the electrons' triplet). So the division by 2M is the right one.
- Downstream, the default screen (`materials_screening()`, `python -m mtft screen`) moves with it: LaH₁₀ λ_eff 2.092 → 1.444 and T_c
  930.1 → 750.5 K, H₃S 743.1 → 598.7 K, and the O/Ba proxy for YBCO 125.6 → 109.8 K; MgB₂, Al and Pb are unchanged. The screen's other weakness, the form
  Θ e^(−1/λ) without the (1 + λ) renormalisation or μ*, is a model choice rather than part of this correction and stays as it is (DIAGNOSTIC).
- Found on 2026-10-28 in a review of Roger's Tano repository workbook, where the same doubling affects seven rows outside the package
  (H, N, O, F, Cl, Br, I). `tests/test_cc71_tano_diatomic.py` (11 tests).
