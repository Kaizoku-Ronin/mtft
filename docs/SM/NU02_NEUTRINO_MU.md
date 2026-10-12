# NU-02 — how small can mu be? The arithmetic geometry of the partners (2026-10-02)

Module `mtft.matter.neutrino_mu`; tests `tests/test_nu02_neutrino_mu.py`.

1. **The modular parametrisation from scratch (NUMERICAL to 1e-12).** 143a's coefficients from the 2-division polynomial (a_11 = a_13 = -1, the (+,+)
   newform); phi(tau) = sum a_n q^n / n has integral periods in the basis (Omega+, Omega+ tau_E), Omega+ = 1.9699231646, tau_E = 1/2 + 1.0232746 i
   (Manin constant 1); phi o W_143 = phi to machine precision.
2. **The partners' points (NUMERICAL, certified by (j(tau), j(143 tau)) pairs).** For each nonzero 2-torsion point P_eta the fibre phi^-1(P_eta) has exactly
   four points modulo Gamma_0(143); every Atkin–Lehner involution acts on it without fixed points, so it carries the regular representation, like the cusps.
   Real spin structure (P = 1/2): tau ~ 0.0789 + 0.0276 i, 0.4629 + 0.0194 i, 0.9211 + 0.0276 i, 0.5371 + 0.0194 i — depth rho = 3.589 below the unit
   horocycles of the cusps, stretch s = 0.1965. Complex pair: rho = 3.397, s = 0.2049. Unit-area separations sigma = pi d^2 = 0.768 and 0.996.
3. **Suppression (ESTIMATE).** mu/M ~ exp(-beta sqrt(r)), r = A_C/A_E, with beta from the bulk geometry: 0.589 - 0.768 (real), 0.732 - 0.996 (pair) for constant
   or uniform stretch. The arithmetic separation alone is about one width; the size hierarchy r supplies the exponent; the cusp region adds a further,
   metric-dependent factor. Requirements from the bulk alone: second family (5.7e-5) r ~ 160 - 280 (real) or 96 - 178 (pair); third family (3.8e-9)
   r ~ 640 - 1080 (real) or 380 - 700 (pair). At r = 56 pi (Matter-05's QCD window) the pair gives 1.8e-6 - 6.0e-5.
4. **Flux (EXACT).** The exponent scales as 1/|c|; the cheap Pati–Salam flux (c = 2 on nu^c) halves it; (1, 1, -2, 0, 0) is D-flat with the (4,1,2) condensate
   alone and has c = 1 on nu^c, at 384 localised extras.
OPEN: the smooth metric near the cusps (Arakelov) and exact overlaps; the Dirac pairing M; family structure from Atkin–Lehner sectors in a symmetric background.
