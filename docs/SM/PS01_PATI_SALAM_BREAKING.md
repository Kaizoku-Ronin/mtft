# PS-01 — the Pati–Salam branch: breaking, anomalies, and the first comparison with data (2026-10-02)

Module `mtft.matter.ps_breaking`; tests `tests/test_ps01_breaking.py`. Index-level and tree-level statements unless marked.

1. **Hypercharge (EXACT).** T_3R(B) = +1/2, T_3R(A2) = -1/2, B-L = 1/3 on colour and -1 on A1 give stack hypercharges (1/6, -1/2, 0, 1/2, -1/2) for
   (c, A1, L, B, A2) — the PV-02 assignment — and the Pati–Salam branching of all three blocks.
2. **Anomalies (EXACT).** SU(4)^3, SU(4)^2-Q4 and both Witten counts (24 doublets) cancel. The stack U(1)s are anomalous (SU(2)_L^2 QL = -12,
   SU(2)_R^2 QR = 12, SU(2)_L^2 Q4 = 6, ...); the only anomaly-free combination is Q4 + QL + QR, which no field carries. Both physical extra U(1)s must
   therefore get Green–Schwarz (Stueckelberg) masses; Y sits inside SU(4) x SU(2)_R and is anomaly-free. The axion completion is the open photon gate.
3. **Wilson lines on 143a1 (EXACT).** Flat bundles keep indices, slopes and rank: N_A1 != N_c gives SU(3) x U(1)_{B-L}, N_A2 != N_B gives U(1)_{3R};
   the Matter-03 Higgs pair survives iff N_L = N_B. Rank 5 remains — Wilson lines alone cannot reach the Standard Model.
4. **The Kahler FI drives electroweak breaking (EXACT D-term algebra).** With delta = A_E - 3 eps the stack FI terms are (0, 0, 2 delta, -delta, -delta)
   for (c, A1, L, B, A2); the (1,2,2) block alone cancels them.
5. **The gamma direction drives Pati–Salam breaking (EXACT on the product).** gamma (gamma^2 = -8; primitive part of the graph of the degree-4 modular
   parametrisation phi) is orthogonal to K, so gamma-flux leaves every chirality unchanged; it enters slopes as -8 c kappa. The Pati–Salam-symmetric,
   traceless flux (c, A1, L, B, A2) = (1, 1, 0, -2, -2) at delta = 0 gives group FI (4: -32 kappa, R: +32 kappa), cancelled by the (4,1,2) block alone —
   a Pati–Salam-breaking condensate set by kappa, independent of the electroweak delta. On each fibre {x} x 143a1, O(gamma) is the flat bundle of phi(x),
   so a block with gamma-charge c and exceptional charge e has chi = -4c^2 - e(e-1)/2; the bosonic 4c^2 modes localise on phi^-1(E[c]) — for c = 1 exactly
   the four cusps of X0(143) (all rational; Manin–Drinfeld; E(Q)_tors = 0 from point counts; deg phi = 4). Under the Scherk–Schwarz twist the fermionic
   ones sit on phi^-1(P_eta). Cost: 8c^2 vector-like modes per rank pair. OPEN: the blow-up point lies over the cusp at infinity; nonabelian D-terms of the
   condensate; F-terms; one-loop competition.
6. **b–tau against data (NUMERICAL, SM one loop).** The flux keeps SU(4) unbroken, so Y_d = Y_e and Y_u = Y_nu at the matching scale. SM running gives
   y_b = y_tau at 1.1 x 10^6 GeV (0.56 - 2.3 x 10^6 for m_b(m_Z) = 2.79 - 2.93 GeV) — the PeV scale where K_L -> mu e already bounds the Pati–Salam vector
   leptoquark (2100 TeV for maximal mixing, arXiv:1210.3697). At that scale m_s/m_mu = 0.33 and m_d/m_e = 3.4 (the Georgi–Jarlskog pattern, not equality),
   and Y_nu = Y_u with M_R <~ 10^6 GeV gives a type-I seesaw mass ~ 0.03 GeV. Those two failures are the next gates (an extra Yukawa source for the light
   generations; an inverse seesaw using the extra neutral vector-like states, or suppressed Dirac couplings).

**CC-38 (2026-10-02).** Item 5's necessity claim assumed Matter-04's curve degrees; with curve degrees free (chirality-neutral under the twist) Pati–Salam breaking by one singlet is D-flat without gamma-flux (FS-01).
