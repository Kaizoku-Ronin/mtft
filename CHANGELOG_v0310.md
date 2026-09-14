# mtft v0.31.0 — Standard Model campaign SM-10…14 consolidated (2026-09-14)

Since 0.30.6: vectorised character reduction and twisted-section evaluation (`hym.reduce_gamma0_vec`,
`twisted_section_values`); flux-divisor scan (`hym.family_gram_for_divisor`) with the divisor-independent
spin-structure invariant 7.8e-5; certified twisted Gram spectra by twist character; cubic Higgs spaces and mixed
Yukawa tensors with a character-violating control (`rrspace.cubic_higgs_space`, `cubic_family_space_100`,
`mixed_yukawa`; frozen `x0143_cubic_higgs.npz`); model M2 (up = cubic x untwisted, down = lepton = cubic x cubic)
with frozen normalised tensors (`hym.load_m2_tensors`, DIAGNOSTIC h = 0.2), `mass_matrix`, `mass_ratios`, `ckm`,
`fit_higgs_direction`; `smflux.hypercharge_normalisation` (sin^2 theta_W = 3/8 for M1/M2).
Results: up-sector hierarchy reproduced (m_c/m_t ~ 7e-3, m_u/m_t ~ 1e-5 at the median vs 3e-3, 6e-6); down within
x3; leptons on the naive SU(5) relations (Georgi–Jarlskog factors missing); nu Dirac tensor = up tensor; CKM
generically large — observed mixing realisable but atypical (0.2%); all reports SM-10…14 in outputs/SM.
Base: PyPI 0.30.0 + the 0.30.1–0.30.9 stack.
