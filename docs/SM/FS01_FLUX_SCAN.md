# FS-01 — scanning SU(2)_R-breaking fluxes at the single blown-up cusp (2026-10-02)

Module `mtft.matter.flux_scan`; tests `tests/test_fs01_flux_scan.py`; correction CC-38. Stacks c, A1, L, B, A2 keep their Matter-04 exceptional charges, so indices
and fermion jets are fixed; per stack the curve degree a_s (chirality-neutral under the twist) and the gamma-charge c_s are scanned.

1. **Curve degrees suffice (EXACT).** With no gamma-flux, a = (0, 1, 1, -1, -2) at eps = A_E/3 (or the SU(4)-preserving (1, 1, 1, -1, -5)) is D-flat with one
   Standard-Model-singlet condensate: Pati–Salam breaking with 288 zero modes and every Higgs entry global. **CC-38:** PS-01's 'the gamma direction is needed' assumed
   Matter-04's curve degrees and is superseded.
2. **TK-1 lifted (EXACT).** With curve degrees free, nu^c gamma-charge +-1, SU(4) kept, a global (L, B) Higgs and a D-flat single condensate coexist.
3. **A gamma-localised Higgs sets its sector's eps (STRUCTURAL).** It enters the Yukawa linearly as a plane wave e^{i k(x) w}, k ~ Delta phi', so
   eps_sector ~ (2.1-6.7)|Delta|/sqrt(r); with SU(2)_R broken the up Yukawa uses the (L, A2) Higgs and the down Yukawa the (L, B) Higgs: eps_d/eps_u = |Delta_d/Delta_u|.
4. **The mod-4 rule (EXACT).** With SU(4) kept, Delta_u + Delta_d = 0 mod 4. Near the data's 2.30: 3 (charges 1, 3), 5/3 (3, 5), 11/5 (5, 11); 2 needs (4, 8).
5. **Families (cheapest representatives, |a| <= 3, |c| <= 3):**

| family | curve degrees a | gamma-charges c | eps / A_E | kappa | nu^c | Higgs (up, down) | zero modes |
|---|---|---|---|---|---|---|---|
| minimal | (0, -1, 1, -1, 0) | 0 | 1/3 | 0 | 0 | (0, 0) | 288 |
| ratio 3 | (-1, -1, 1, 0, 2) | (0, 0, -1, 2, 0) | 2/3 | 1/8 | 0 | (-1, -3) | 672 |
| neutrino-optimal, both localised | (-1, -1, 2, 0, 0) | (0, 0, 0, -1, 1) | 2/3 | -1/4 | -1 | (-1, 1) | 416 |
| neutrino-optimal, up global | (-1, -1, 1, 0, 2) | (0, 0, -1, 3, -1) | 1/2 | 1/16 | 1 | (0, -4) | 1056 |

6. **The separation (EXACT).** Ratio 3 together with a charged nu^c needs Higgs charges (2, 6) of opposite sign — e.g. c = (1, 1, 0, -6, 2), 3104 zero modes. The cheap
   options separate: the up/down asymmetry (ratio 3) or the neutrino suppression (charged nu^c), not both.
OPEN: the order-one factors in eps ~ |Delta|; ampleness beyond the test curves (exceptional curve, fibre and X0(143) x O through the point, graphs of +-phi — all
positive for r > 1); the Higgs-sector bosonic spectra at the new curve degrees.
