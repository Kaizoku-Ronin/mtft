# GC-01 — what the gauge couplings fix (2026-10-05)

Module `mtft.matter.gauge_couplings`; tests `tests/test_gc01_gauge_couplings.py`.

1. **Shape versus volume (EXACT).** J^2/2 = A_X A_E - sum eps^2/2 and Re f_i = J^2/2 - kappa F_i^2/2, F_i^2 = 2ab - 8c^2 - sum e^2. With no torus flux, A_X enters only through
   the volume: coupling ratios fix V/kappa, never A_X/A_E. V_cb remains the only measure of the shape.
2. **Matching (NUMERICAL).** alpha_4 = alpha_3, alpha_L = alpha_2, alpha_Y^-1 = alpha_R^-1 + (2/3) alpha_4^-1. The flux-only prediction R = (sum e_R^2 - sum e_4^2)/(sum e_L^2 -
   sum e_4^2) is independent of V/kappa. Data: R = 3.99 at m_Z, 3.678 at the b-tau crossing, 1 at 5.1e13 GeV; alpha_2 = alpha_3 at 9.7e16 GeV.
3. **Quantised (EXACT given f).** Two-point designs: R = 1/9, 3/7, 1, 7/3, 9; single point: 1/5, 1, 5. None reaches 3.678. Matching scales 8.9e11, 5.1e13, 1.4e14, 2.1e14 GeV.
4. **Tension (ESTIMATE).** At those scales y_b/y_tau = 0.66-0.71: b-tau would need 30-35 per cent corrections. The 1.1e6 GeV scale used since PS-01 is not where the gauge
   couplings put Pati–Salam (CC-48 caveats PS1-4); results evaluated there inherit the caveat.
5. **Which home (EXACT by symmetry).** The torus involution fixes P1, so the top's torus coordinate is critical at the four 2-torsion points; no automorphism fixes a cusp,
   so the curve coordinate has no symmetry selection.

Open: threshold corrections (vector-like pairs, adjoints, KK towers) that could reconcile the scales; the shape enters the couplings only there.
