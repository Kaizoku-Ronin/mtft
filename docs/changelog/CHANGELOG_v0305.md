# mtft v0.30.5 — SM-07: cubic-character twist and the charged-lepton tensor (2026-09-13)

`surface.rrspace`: frozen cubic/sextic Conrey-character spaces (133, 56 and conjugates, q^1200) with AL matrices;
`conrey_chi` (convention fixed by a direct-vs-reduced test, 1e-14 vs 1.7), `cubic_values`, `cubic_family_space`
(dim 3, pure, 37/40 conditions), `lepton_yukawa_M1` (3 x 3 x 18, residual 8e-9, character identity
chi_133 chi_56 = chi_13 verified by the product test), `family_equivalence_residual` (alternating GL(3)xGL(3) fit).
Result: with a cubic twist on the singlet stack the lepton tensor is inequivalent to the down tensor and its
family transpose (2e-2 vs 6e-16 control): the down–lepton mass-ratio degeneracy forced by quadratic twists is
lifted.  Recorded: the Higgs directions are flat moduli at leading order; predictions are direction-independent
relations.  Full fast suite before build; clean tree; MANIFEST covers *.npz.  Base: PyPI 0.30.0 + 0.30.1–0.30.4.
