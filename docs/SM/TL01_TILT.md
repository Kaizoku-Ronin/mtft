# TL-01 — the gamma tilt is a rotation (2026-10-07)

Module `mtft.matter.tilt`; tests `tests/test_tl01_tilt.py`.

1. **The tilt is a rotation (EXACT).** The gamma bundle restricts to the flat torus bundles N_{c phi(x)} — a Poincare-type pullback whose curvature mixes the curve and
   torus directions (odd under w -> -w). At the jet point it makes the x- and w-jets non-orthogonal: the physical jets are a rotated pair (LV-01's pure lever).
2. **Block rotations cancel (EXACT).** The quark doublet block is shared by up and down, so its tilt cancels in U_u^dagger U_d; right-handed tilts drop out of
   Y Y^dagger; V_us = V_ub = 0 (checked to 1e-17). The tilt must act inside the vertex and differ between up and down (EPS-01's eps_d/eps_u split).
3. **V_ub measures the down quarks' share of V_cb (EXACT at leading order, 3% checked).** With sector-dependent aligned tilts, V_ub = V_us |s23_d|, so c = |s23_d| / V_cb.
4. **The programme's numbers (NUMERICAL).** s23_u/s23_d = (m_c/m_t)/(z m_s/m_b) = 0.19/z, so c >= 1/(1 + 0.19/z). At DP-01's z = 0.25: c >= 0.565, outside the
   pre-registered window [0.35, 0.45] — **fail**. The window needs z <= 0.157 with nearly opposite 2–3 phases.
5. **Verdict (STRUCTURAL).** c = 0.40 needs LV-01's local arm or stronger up-sector dressing; both need the condensate's SU(2)_R split at the jet point.

CC-58 caveats CA-1 (and so the V_us window in FH-3's concordance).
