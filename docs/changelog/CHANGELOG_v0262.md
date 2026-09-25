# mtft v0.26.2 — CC-21 and CC-22 from Astra's electron–photon audit (2026-09-07)

## CC-21 — `lattice.random_su_n` produced no SU(N) matrices
`Q[0,:] /= det**(1/N)` rescales the determinant by det^(-1/N); every draw
failed det = 1 (60/60 at N = 2, 3, 5).  Fixed by dividing the single row by the
full determinant; `random_su_n_defective_v0261` retained for reproducing
earlier runs.  Consequence: any hot-start result produced through
`LatticeGauge` with random links was not an SU(N) simulation and must be
rerun before interpretation; cold starts are unaffected.  Found by Astra
(SM-AUDIT-02), verified here.

## CC-22 — Higgs quartic convention (documentation; no numeric change)
`HIGGS.lambda_quartic = (γ/Ω)²/2` satisfies m_H² = λv²/2 and its 0.074% is
against 2m_H²/v² in that convention.  Astra's SM-AUDIT-01 rewrote it as
(γ/Ω)²/8 (the m_H² = 2λv² convention); that is a convention change and would
silently alter a stated formula, so it is NOT applied.  Added
`lambda_quartic_pdg` = λ/4 with the convention documented.

Tests: `tests/test_cc21_cc22.py` (2).  Release pin bumped four-way.
