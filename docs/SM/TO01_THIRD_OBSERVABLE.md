# TO-01 — a pre-registered third observable (2026-10-05)

Module `mtft.matter.third_observable`; tests `tests/test_to01_third_observable.py`.

**Pre-registration (fixed before computing).** Observable: V_cb predicted from the sectors' own mass ratios using only the design's structure — no fitted parameter, independent
of where the top sits. Pass window: [0.030, 0.055].

1. **The second family is induced (EXACT within the design).** Its direct entry costs three units of curve-direction breaking, at most 2.4e-4 (r = 18.6) to 1.8e-5 (r = 105),
   far below m_c/m_t = 0.0036 and m_s/m_b = 0.019. So m2/m3 = theta_L theta_R, and with left–right symmetric tails theta_L = theta_R = sqrt(m2/m3).
2. **Prediction (NUMERICAL).** V_cb = |sqrt(m_s/m_b) - e^{i phi} sqrt(m_c/m_t)| lies between 0.078 and 0.199 for every phase: outside the window by at least a factor 1.9. Fail.
3. **What it corrects (CC-50).** A geometric angle common to both sectors cancels in the CKM matrix; V_cb here is fixed by mass ratios, not by the top's distance. SP-2, DD-4 and
   CA-2's r-from-V_cb become caveats; CA-3's Fricke candidate loses its premise and is superseded.
4. **What the design needs.** Left–right asymmetric tails: theta_L^d ~ 0.05, theta_R^d ~ 0.38, a ratio near 7.7. [Corrected by CC-51: theta_L is set by the right-handed block's tail, so Q's tail must be the longer one; VS-3's ordering points the wrong way.]

Open: tail ranges as a function of the blocks' curve degrees; the up–down asymmetry from the SU(2)_R-breaking fluxes.
