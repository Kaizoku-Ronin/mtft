# VE-01 — the vacuum energy of the twisted compactification: minus a square, and the point where it vanishes (2026-10-17)

Module `mtft.matter.vacuum_energy`; tests `tests/test_ve01_vacuum_energy.py`; every number from `report()`. Addresses puzzle I.6 (a vacuum
energy from the parent action) with what the programme's own terms give, and corrects SC-4 (CC-70). The question Roger put was the vacuum
energy of "our anti-de Sitter space"; the programme has three: the modular-time plane H (Euclidean AdS_2, the space of ratios of two clocks),
the curve X0(143) = Gamma_0(143) \ H that carries its metric, and the four-dimensional vacuum. The first two are arenas; the third is the one
with an energy.

## 1. The perfect square (EXACT, one loop, product approximation, through a_1)

For the U(8) multiplet with torus Wilson lines U (holonomy U^w around the winding w of 143a1), the Scherk–Schwarz supertrace factorises block
by block as (1 - eta(w)) e^{i (theta_i - theta_j) . w} (Matter-05), so summing the blocks gives

    V(U) = -c sum_{eta(w) = -1} K(w) |Tr U^w|^2,      c > 0,   K(w) = |w|^-8 + (2 pi/3r) |w|^-6 + ... > 0  (summed over the lattice: S0 (1 + 2.78/r)).

MS-5 makes the block coefficient flux-independent through the a_1 term of the curve heat kernel, so the square holds to relative order 1/r;
flux enters at a_2 ((B t)^2). Consequences: V <= 0 everywhere; V = 0 exactly where Tr U^w = 0 on every odd winding; that locus is the
global **maximum**. Checked against the block-by-block sum to 1e-9 at random Wilson lines.

## 2. The Z/2 table (EXACT)

For the Z/2 Wilson lines of the three Pati–Salam groups, with `Background.twist` deciding which blocks have light bosons and which light fermions:

| W on (4, L, R) | bosons untwisted | fermions untwisted | (Tr W)^2 | V / (c S0) | PA-01's verdict |
|---|---|---|---|---|---|
| (+, +, +) jet vacuum | 64 | 0 | 64 | -64 | minimum |
| one SU(2) flipped (4 cases) | 40 | 24 | 16 | -16 | flat at one loop |
| (-, +, +) or (+, -, -) | 32 | 32 | 0 | 0 | maximum |

N(bosons) - N(fermions) = (Tr W)^2 in all eight cases. The 0 row is PA-01's matter-untwisted background: the Bose–Fermi-balanced point of the
compactification — the "super no-scale" condition n_F = n_B of Kounnas–Partouche, met here by the Pati–Salam ranks 4 = 2 + 2 — where the
one-loop vacuum energy vanishes. PG-6's saddles sit at V = -32 c S0, half the jet value, by 143a1's real structure (Re tau = 1/2).

## 3. The balanced point is the top (EXACT)

| | Hessian signature (16 Cartan Wilson lines) | negative directions | flat directions |
|---|---|---|---|
| jet vacuum | (0 negative, 2 zero, 14 positive) | — | the overall U(1), exactly |
| balanced point | (2 negative, 14 zero, 0 positive) | the relative U(1) T = diag(1_4, -1_2, -1_2) on both torus cycles | overall U(1) exactly; the other 12 are maxima at quartic order (V = -c K u^4, fitted power 3.96) |

The two tachyonic curvatures are exactly minus the jet vacuum's Wilson-line curvatures (PG-4), with PG-6's theta-coordinate ratio 166.935/28.893.
Downhill from 0 lie the 16-points and then the jet vacuum. So the balanced point is a vacuum only in a construction with **no continuous torus
Wilson lines at all**. What the programme has for that:

- T is anomalous — SU(2)_L^2 T = 12, SU(2)_R^2 T = -12, SU(4)^2 T = grav^2 T = 0 (PS-01's table) — so a Green–Schwarz mass M_St for it (I.2, open)
  lifts the two tachyons iff M_St exceeds the jet vacuum's Wilson-line masses, 1.85–2.2 and 3.4–4.1 M_KK on the measured background (SC-2). It does
  nothing for the twelve non-abelian directions.
- Atkin–Lehner left–right parity (the `two-cusps-dparity` background: theta_4 -> -theta_4, theta_L <-> theta_R) freezes theta_4 at a half-period and
  leaves theta_L = theta_R = b continuous; that direction's projection on T is -4 b, so it stays tachyonic.
- A magnetised torus removes the moduli but costs 24 chiral families per flux unit (TP-1).

The history is a caution in the same direction: Moore's Atkin–Lehner proposal (1987) for a vanishing one-loop cosmological constant met no-go
theorems; Kachru–Kumar–Silverstein's cancellation at one and two loops acquires a non-perturbative value in the dual (Harvey 1998).

## 4. The curve corrections — CC-70 (EXACT at a_1; a_2 and geodesics ESTIMATE/NUMERICAL)

SC-4 said the closed forms of Matter-05 are "O(1) uncertain" because 56 pi/r = 0.63–0.77 is not small. The torus images weight the curve's heat
kernel at proper time <t> = 0.11 A_E, so the curve's expansion parameter is <t>/R_C^2 = 0.11 x 56 pi/r = 0.07–0.085, and term by term, relative
to the leading closed form at r = 230–280:

| term | V | Wilson-line mass^2 | standing |
|---|---|---|---|
| a_1 = MS-5's -chi/12 per block (topological) | +(pi/3)(a_1/N)(S1/S0)/r = +2.78/r = +1.0 to +1.2% | +1.0 to +1.7% (two directions) | EXACT given MS-5 |
| a_2 (scalar-Laplacian value (1/60 pi) int K^2 dA per d.o.f. as the size) | 5–8 x 10^-4 | — | ESTIMATE |
| shortest closed geodesics of Gamma_0(143), length 2 arccosh 2 = 2.634 | 3–5 x 10^-6 per primitive class | — | NUMERICAL |

The closed forms are good to about a per cent on the measured background; the geodesic sum SC-4 asked for would change nothing visible. The
O(1) uncertainty that remains in SC-2 is lambda/16 = 4.5–6.5 — higher loops — and no heat kernel fixes that. CC-70 also corrects the exponent in
MS-01's definition of S_k (|w|^-(8-2k): S0 multiplies the A_X term, S1 = Matter-05's S the a_1 term, as the quoted values already had it).

## 5. The curve's own rung (EXACT in the cusped metric; the first discrete mode DIAGNOSTIC)

The hyperbolic Laplacian of Gamma_0(143) has no eigenvalue in (0, 1/4): Selberg's eigenvalue conjecture is a theorem at this level (Booker, Lee
and Strömbergsson 2018, Gamma_1(N) for N <= 880; Gamma_0(143)'s spectrum lies in Gamma_1(143)'s), and 1/4 is the Breitenlohner–Freedman bound of
the AdS_2 the curve is a quotient of. So in the cusped metric the curve's KK tower starts at

    1/(2 R_C) = sqrt(56 pi / r) / (4 pi) M_KK = 0.063–0.070 M_KK,

7.1–7.8 times below M_SS — the ladder's lowest rung, missing from SC-01 — with the first Maass form at sqrt(0.39)/R_C = 0.079–0.087 M_KK (SPEC-01,
DIAGNOSTIC). In the smooth blown-up metric the lowest eigenvalue is not computed. GC-01's thresholds should be checked against this rung; and
because the curve's tower starts below M_SS, the one-loop count weighs whole 8D multiplets, not 4D zero modes — which is why section 2 counts
8D degrees of freedom.

## 6. The ledger (STRUCTURAL; numbers ESTIMATES under convention C1)

On the measured background (r = 230–280, M_KK = 8.9e11–2.1e14 GeV; rho_Lambda = 2.5e-47 GeV^4):

| term | sign | size | / rho_Lambda |
|---|---|---|---|
| classical internal curvature, +M_P^2/R_C^2 | + | 7e58–5e63 GeV^4 | 10^105–10^110 |
| one-loop twisted Casimir at the jet vacuum, 32 complex directions | − | 3.7–4.5 M_KK^4 = 2e48–9e57 GeV^4 | 10^95–10^105 |

Negative internal curvature is positive four-dimensional energy (Douglas–Kallosh; Silverstein's de Sitter constructions use it), and Douglas–
Kallosh's no-go — everywhere-negative internal scalar curvature forbids a maximally symmetric four-dimensional vacuum without large warping or
large corrections — is PV-01's and PG-4's runaway seen from the other side. The classical term beats the loop by 10^5–10^10. The C3X AdS control
(`compactification.c3x_ads_control`) reaches AdS_4 only with an added negative bulk term and then has ell_4 < sqrt(3) R_C = 2.0–2.2 sqrt(A_E): no
scale separation, as Lüst–Palti–Vafa's distance conjecture expects.

Nothing here is a vacuum; I.6 stays open. What is new is that the anti-de Sitter ingredient the programme has — the twisted Casimir sink — is
cancelled exactly at one point of the Wilson-line torus, and that point is the top of the potential. A vanishing one-loop energy at a vacuum
needs the Z/2 twist realised rigidly (PA-01's open alternative), now motivated twice: tree-level Yukawas (YK-01) and Tr W = 0.

## Package changes

New `mtft.matter.vacuum_energy` (potential, block_sum, z2_table, hessian, signature, quartic_profile, t_anomalies, dparity_projection,
stueckelberg_criterion, curve_corrections, curve_kk_rung, energy_ledger, report) and `tests/test_ve01_vacuum_energy.py` (6 tests); ledger VE-1..VE-6
and CC-70 (294 claims); SC-4 carries CC-70's caveat; puzzle I.6 updated (open); wave VE-01 registered under Piece IV; MS-01's S_k exponent fixed.

Sources: Kounnas and Partouche, arXiv:1607.01767; Moore, Nucl. Phys. B 293 (1987) 139, via Gannon and Lam, hep-th/9201028; Kachru, Kumar and
Silverstein, hep-th/9807076; Harvey, hep-th/9807213; Booker, Lee and Strömbergsson, arXiv:1803.06016; Douglas and Kallosh, arXiv:1001.4008;
Silverstein, arXiv:0712.1196; Lüst, Palti and Vafa, arXiv:1906.05225.
