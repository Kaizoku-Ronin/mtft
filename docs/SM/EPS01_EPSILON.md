# EPS-01 — what sets epsilon (2026-10-02)

Module `mtft.matter.epsilon`; tests `tests/test_eps01_epsilon.py`. Local model at the blown-up point: x along X0(143), w along 143a1; families x^2, x w, w^2 with
torus charges 0, 1, 2; the Yukawa entry (q, q') reads the torus-charge (q + q') component of the D-term the Higgs mixes through (PB-01).

1. **Rank one (EXACT).** A torus-constant D-term couples only q = q' = 0.
2. **The gamma-localised condensate keeps it (EXACT).** One band of the fibre family N_{c phi(x)} is a plane wave along the torus: its density is torus-constant.
3. **Rotating the families does not help (EXACT).** Mixing x^2, x w, w^2 among themselves keeps an outer product.
4. **Rank three (EXACT structure).** Y_{qq'} ~ D_{q+q'}; a condensate spread over jets of degree <= 2 has D-charges in [-2, 2]: texture [[1, eps, eps^2], [eps, eps^2, 0], [eps^2, 0, 0]],
   determinant -eps^6, masses (1, eps^2, eps^4) for generic coefficients, mixing (eps, eps, eps^2). Epsilon = |D_1|/|D_0| is the condensate's torus-charge misalignment
   at the point — a vacuum-alignment quantity. Against data (eps_u = 0.060, eps_d = 0.139, order-one coefficients): m_u/m_t 9.1e-6 (7.4e-6), m_c/m_t 5.1e-3 (3.6e-3),
   m_d/m_b 2.4e-4 (9.4e-4), m_s/m_b 2.7e-2 (1.9e-2), |V_us| 0.11 (0.225), |V_cb| 0.15 (0.041), |V_ub| 0.020 (0.0037).
5. **Up versus down (EXACT structure).** A rank-one aligned condensate gives every D^a one profile: eps_u = eps_d. The data need eps_d/eps_u = 2.30.
6. **Parity (EXACT).** The eps-linear entries need w -> -w broken; the gamma-flux breaks it.
7. **Size (ESTIMATE).** Arakelov density at the cusp from the frozen periods: rho = K(inf, inf)/(13 * 4 pi^2) = 0.02434, so phi' = 3.217/sqrt(r) and
   eps_gamma ~ (2.1-6.7)|c|/sqrt(r); eps_d = 0.14 needs r ~ 230-2300 (c = 1), the large-r window the neutrino sector also prefers.

**Toolkit finding (TK-1, EXACT).** A Pati–Salam-symmetric traceless gamma-flux that keeps the bidoublet free of gamma-flux (a global Higgs, as PB-01 needs) has
c_L = c_R, hence gamma-charge 2 c_4 on nu^c — never +-1. NU-02's optimum and PB-01's Higgs are incompatible within Pati–Salam-symmetric fluxes; an SU(2)_R-breaking
gamma-flux (B != A2) is the escape — the same ingredient item 5 needs for eps_d != eps_u.
