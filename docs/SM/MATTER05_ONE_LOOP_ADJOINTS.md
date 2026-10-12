# Matter-05 — one-loop masses of the bosonic adjoints, and a first look at Pati–Salam (2026-10-02)

Modules `mtft.matter.one_loop_adjoints`, `mtft.matter.pati_salam`; tests `tests/test_matter05_one_loop_adjoints.py`.

## Result: the adjoint question closes positively at one loop
Under the Matter-04 Scherk–Schwarz twist (fermions antiperiodic on an even spin structure of 143a1), twisted supersymmetry still pairs the boson and fermion
spectra on X0(143) for every flat connection and holomorphic Phi, so the one-loop potential factorises as V = -(n/32 pi^2) int dt t^-3 Tr_C(t) P(t) with
P(t) = (A_E/2 pi t) sum_{eta(w) = -1} e^{-|w|^2/4t} > 0.
- **Sign theorem (EXACT on the product, any product metric, gauge-multiplet loops).** All 27 bosonic adjoints per generator have strictly positive mass^2.
  Curve Wilson lines: method of images with the positive scalar heat kernel on the universal cover (every species reduces to it on a curve, up to
  a-independent index constants); the Hessian is a sum over closed geodesics of X0(143) weighted by (alpha(a).[gamma])^2. Phi: the commutator potential
  |alpha(phi)|^2 |omega|^2 >= 0 raises every level (min-max). Torus Wilson line: only odd windings survive, each with a positive weight.
- **Matter loops:** positive for Phi and for the torus Wilson line (without background gauge Wilson lines); for the curve Wilson lines the flux blocks'
  image sums carry magnetic phases and the sign is OPEN.
- **Closed forms** (leading order in A_E/R_C^2): m_Phi^2 = (n/2 pi^4) g_4^2 (A_C/A_E^2) S_eta Sum|alpha|^2;  M_theta^2 = (6n/pi^4) g_4^2 (A_C/A_E^2) T_eta Sum|alpha|^2.
- **143a1's three even spin structures split 1 + 2**, as its 2-torsion does over Q: antiperiodic on the tau-cycle (real point 1/2): S = 2.107219,
  T = (0.473605, 1.633614); the complex-conjugate pair (tau/2, (1+tau)/2): S = 3.238510, T = (0.770478, 2.468032) — Phi masses 1.537 times heavier.
  Torus Wilson lines are 2.7 to 9.3 times heavier (in mass^2) than Phi.
- **Curve Wilson lines are the lightest:** their image sums need a closed geodesic; the shortest in Gamma_0(143) has trace 4, length 2 arccosh 2 = 2.6339
  (e.g. [[19,-2],[143,-15]]), area 56 pi. QCD with six flavours: b0 = 7 below every adjoint mass, -6 while only the 26 real curve-Wilson-line octets are
  light, -20 with all 54 — so the window above the lightest adjoint must be short: 143a1 should not be much smaller than X0(143)'s curvature radius
  (r = A_C/A_E of order 56 pi; the blow-up correction 1/(3r) is then about 0.2 per cent).

## Pati–Salam first look (index counts only)
At the Matter-04 minimum the stack classes coincide (A1 = c, A2 = B). Without separating Wilson lines the gauge group enhances to U(4) x U(2) x U(2):
4 = colour + A1, 2_L = L, 2_R = B + A2. Inter-group blocks under the twist: (4,2,1) — 72 Weyl, three net families of (q, l); (4bar,1,2) — 72 Weyl, three net
families of (u^c, d^c, e^c, nu^c) with (c,B) = d^c, (c,A2) = u^c, (A1,B) = e^c, (A1,A2) = nu^c; (1,2,2) — 144 Weyl in a real representation, whose bosons hold
the Matter-03 Higgs pair. Total 288. Separating Wilson lines are flat and leave these indices unchanged: they are Wilson-line breaking to
SU(3) x U(1)_{B-L} x U(1)_{3R}. OPEN, next: breaking to the Standard Model, the extra U(1)s, anomalies, Yukawa structure and mass relations, the vector-like pairs.
