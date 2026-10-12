# LEP-01 — the charged-lepton texture (2026-10-02)

Module `mtft.matter.leptons`; tests `tests/test_lep01_leptons.py`; command `python -m mtft matter leptons`; correction CC-39.

1. **Georgi–Jarlskog at the matching scale (NUMERICAL).** At 1.1e6 GeV: m_mu/m_s = 3.06, m_e/m_d = 0.295, m_tau = m_b. A common rescaling of eps cannot fix both
   light families (1.76 needed for the second, 0.74 for the first).
2. **Y_e = Y_d is exact with SU(4) kept (EXACT).** Leptons and down quarks share the Higgs entry (L, B) and its heavy modes; with equal gamma-charges, curve degrees and
   exceptional charges on colour and A1 their jets coincide too. Checked across 2,037 SU(4)-preserving scan backgrounds.
3. **The Georgi–Jarlskog-friendly order (STRUCTURAL).** SU(4)-breaking gamma-flux reshapes the lepton jets at order eps; the third family moves only at relative order
   eps^2. b–tau survives, the light-family ratios can move by order one.
4. **How different leptons must be (MODEL).** With equal overlaps (Yukawa entries = autocorrelations of the condensate's jet amplitudes s): up (1, 0.025, -0.003),
   down (1, 0.221, 0.029), leptons (1, 3.26, -0.31). Quarks are charge-0 dominated; leptons need the charge-1 jet to dominate.
5. **Koide (NUMERICAL).** Q = 0.666664 (leptons, pole), 0.746 (down), 0.888 (up). The texture does not single out 2/3. The relation s1 s2 = -1, noticed after fitting,
   predicts m_e 5.5 per cent low and Q = 0.6672: chosen after looking, so numerology.
6. **CC-39.** Vertex gamma-charges sum to zero, so the up and down vertices differ in two charges (up (1, 0, 1), down (1, 2, 3)) and the Higgs entries in curve degree.
   RT-1's 'exactly 3' becomes 'near 3'; FS-4's plane-wave phase is a gauge artefact (the gauge-invariant effect is wavefunction reshaping).

OPEN: the lepton dressing by the condensate (it points in the lepton direction of the 4); the effective potential for s — the Legendre transform of the generating
functional — which would fix s, eps and every texture from first principles.


**CC-61 (2026-10-08).** Item 4's "leptons need the charge-1 jet" (LEP-4) is retracted: the equal-overlap texture has several exact solutions per sector, and the leptons also fit with the charge-0 jet dominant, (1, 0.339, 0.031) (LEP-6). The earlier reading was the optimizer's path, exposed when a SciPy upgrade changed which solution it reached first. `fit_alignment` now returns a stated canonical representative (the rho branch nearest 1/2, then the smallest coefficients) and `all_alignments` returns the full set.
