# mtft v0.30.7 — SM-10: twisted-sector normalisations (2026-09-14)

`surface.hym`: `reduce_gamma0_vec` — vectorised Gamma_0(143)-only reduction tracking the composite matrix modulo 13
exactly (character evaluation; matches the per-point reducer to 1e-15; test added) and `twisted_section_values`
for the chi_13 (d^c, H_d), cubic (L) and sextic (e^c) sections at arbitrary point sets, with the low-height mask.
Results recorded (DIAGNOSTIC): HYM norm spectra (7.8e-5,6.9e-4,1) Q, (1.9e-4,1.1e-3,1) d^c, (2.6e-5,3e-4,1) e^c,
(8.2e-3,0.11,1) L — the twist character sets the hierarchy strength; quark sectors over-hierarchical, lepton
sector realistic in spread; Georgi–Jarlskog over-corrected (65, 290 vs 3, 1/3).  Base: PyPI 0.30.0 + 0.30.1–0.30.6.
