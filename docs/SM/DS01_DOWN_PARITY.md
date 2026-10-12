# DS-01 — the down sector is parity-symmetric (2026-10-07)

Module `mtft.matter.down_parity`; tests `tests/test_ds01_down_parity.py`.

1. **What the condensate touches (EXACT).** At first order the nu^c condensate's D-term reaches u^c (+1), e^c (+1), nu^c (+2), L (-1) and H_u (-1), and none of Q, d^c,
   H_d — CD-1's pattern. The down sector has no first-order parity breaking; the Cabibbo angle is an up-sector quantity.
2. **The floor (EXACT).** m_s/m_b = tan theta_L^d tan theta_R exactly, and V_ub = V_us sin theta_L^d (TL-3), so c >= (m_s/m_b)/(V_cb tan theta_R) and
   |V_ub| >= V_us (m_s/m_b)/tan theta_R. For theta_R <= 45 deg: c >= 0.47 (0.42–0.52), |V_ub| >= 0.0043 (0.0039–0.0048) against 0.0037 (exclusive)–0.0042 (inclusive).
   The pre-registered window [0.35, 0.45] is missed at central inputs, grazed within uncertainties.
3. **What the data say (NUMERICAL).** tan theta_R = 1.17 (49.5 deg), sin theta_L^u = 0.050, 2–3 phase 48 deg, z = 0.063.
4. **The tails cannot do it (ESTIMATE).** T(Q)/T(d^c) = 71 needed; curve flux gives at most 13.8 (FL-01), so c >= 0.91 and |V_ub| >= 0.0084: 2.3 times too large.
5. **Verdict (STRUCTURAL).** The design needs first-order parity breaking in the down sector or a different home for the left-handed top. The fifth open failure.

CC-59 supersedes TL-5.
