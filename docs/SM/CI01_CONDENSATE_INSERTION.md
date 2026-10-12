# CI-01 — the condensate insertion (2026-10-06)

Module `mtft.matter.condensate_insertion`; tests `tests/test_ci01_condensate_insertion.py`.

1. **The projector is rank one (EXACT).** The condensate's internal direction is fixed by its bundle, so Sigma Sigma^dagger = |Sigma|^2 P_lepton everywhere. An entry's
   lepton/quark ratio is 1 + rho w_ij: rho the universal projector-to-singlet coefficient ratio, w_ij the share of the quark entry carried by the condensate channel.
2. **Tails are universal (EXACT).** Tail-carried entries involve only zero-mode profiles fixed by flux that commutes with SU(4): lepton factor 1 (times dressing <= 1).
   The 2–3 entries are tail-carried (TO-01, FL-01), so the muon's 3.09 forces y_tau/y_b <= 0.57 at matching against 1.40–1.53 — factor >= 2.46 (1.76 with b–tau unification).
3. **The entry pattern (NUMERICAL, running-independent).** 2–3 entries x1.758 relative to the 3–3 entry, 1–2 entries x0.959: lepton tails 1.76 times quark tails with
   local entries untouched. Predicts f_1 f_2 = f_3^2; data 0.919 (within m_d/m_s's 8%).
4. **The condensate channel is too short-armed (ESTIMATE).** Gaussian jets (l_1 = 0.268), holomorphic 1-jet condensate at P1, singlet coefficient zero (most favourable):
   the channel's 2–3 excess is bounded by the tail lever arm 2d/(n l_1). Needed 5.9 (f_3 = 1.53) to 7.7 (1.40), impossible below f_3 = 1.043; the n = 3 tail gives 2.4–3.0.
5. **Verdict (STRUCTURAL).** Neither dressing (QL-01) nor insertion gives the lepton pattern. It needs SU(4) breaking in the zero modes — e.g. a discrete B–L Wilson line.

CC-53 supersedes QL-3. Open: the B–L Wilson line; the overall b–tau normalisation (GC-01's threshold problem).
