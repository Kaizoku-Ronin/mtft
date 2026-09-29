# PV-01 — parent actions and vacua: the 8D anomaly gate, the B1 finite search, and the flux-energy landscape (29 September 2026)

Tools built for the parent-action ledger (`PARENT_ACTION_REQUIREMENTS.md`, class B1 and section C), and what they returned on first use.
Every statement below is EXACT unless it says otherwise; the finite scans are exact within their stated boxes and are not theorems.

## 1. `research/anomaly8d.py` — the ten-form anomaly of chiral 8D content
- Per positive-chirality Weyl fermion I_10 = [A-roof(T) ch(F)]_10; representations by lambda-ring rules on the truncated Chern character
  (Adams operation psi^2 for Lambda^2 / Sym^2).  Identities reproduced: tr_{Lambda^2} F^5 = (N - 16) tr F^5 + 5 tr F tr F^4 + 10 tr F^2 tr F^3
  (Sym^2: N + 16); tr_{Lambda^2} F^3 = (N - 4) tr F^3 + 3 tr F tr F^2 (the 4D rule); adjoint and every V (+) V*: zero.
- SU(N >= 5), content (n_F, n_A, n_S): the irreducible quintic vanishes iff n_F + (N - 16) n_A + (N + 16) n_S = 0, and the entire remainder is
  s_3 (alpha s_2 + beta p_1): ONE 2-form with X_6 = tr F^3 completes the Green–Schwarz cancellation.  The cubic Casimir that obstructs the
  6D lift (R2C-04, `hom_type_obstruction`) is the Green–Schwarz factor in 8D.  Exact cancellation without Green–Schwarz forces n_F = n_A = n_S = 0.
- U(N): the c_1 terms need a second field (axion, X_2 X_8): `gs_decomposition`; U(6) with Lambda^2 + 10 fundamentals verified in explicit
  Chern roots (`specialise`).  There is no pure gravitational anomaly in 8D.
- Split-stack parity lemma: on a split background, blocks of Lambda^2 / Sym^2 of one stack are L_i^2 with index 4 a_i b_i (even); adjoint
  blocks are vector-like; only L_i^{+-1} and mixed L_i L_j blocks can be odd.  With the trivial colour stack of S1/S2 the coloured families
  come only from mixed blocks, so u^c and d^c need two singlet stacks (ledger row 3 made concrete).

## 2. The B1 finite search (`b1_scan`, `hypercharge_assignment`, `sm_net_content`, `b1_trilemma`)
Content Lambda^2(N) + (16 - N) N (quintic-free), N = 5 + k, on stacks c (U(3)), L (U(2)) and k singlet stacks with bidegrees in a box;
acceptance = net left-handed content EXACTLY three Standard-Model families plus hypercharge-neutral singlets and vector-like pairs, for
some rational stack hypercharges (exact Gaussian elimination; the free directions are pinned exhaustively at SM values).
- k = 1 (U(6)), k = 2 (U(7)), box 3, and the Sym^2 class: NO Standard-Model assignment (the lepton-doublet count cannot be 3).
- k = 3 (U(8) = U(3) x U(2) x U(1)^3, Lambda^2(8) + 8 x 8): box 3 gives 60 assignments, box 4 gives 156, each with exactly three families;
  one class has a single vector-like Higgs-type doublet pair and only integer-charged singlets.  These are the first chiral, quintic-free 8D
  contents that land the Standard Model on X0(143) x 143a1 as net indices.  (Net indices only; vector-like pairs beyond the ones listed
  are not excluded; 4D anomalies are automatic for this content once the extra U(1)s are Stueckelberg-massive.)
- The trilemma follows the move.  With (Y) the bidegree Yukawa rule, (H) a slope-free Higgs block L_L L_s^-1 for both u and d, (C) a
  colour-matchable slope of the rank-5 rest:  box 3 — Y+H (6) and Y+C (6) are disjoint;  box 4 — H and C can both hold but only at
  DIFFERENT shapes (r_H = 2, r_c = 5/2 and the mirror); no assignment in either box has Y, H and C at one ratio.  Catalogue exemplars:
  `B1_U8_LAMBDA2_YH` (Higgs at r = 2, colour never), `B1_U8_LAMBDA2_YC` (colour at r = 3/2, Higgs never), `B1_U8_LAMBDA2_HC` (2 vs 5/2).
- Not decided here: the mu-stability of any non-split bundle with these Chern classes (the rank-2 lemma was proved for Hom-type blocks; its
  tensor-block analogue is open), global anomalies (Omega_9^Spin), the Green–Schwarz completion itself, and every scale.

## 3. `research/vacuum_energy.py` — the flux landscape (section C of the ledger)
- int_S |F|^2 = 4 pi^2 (a^2 / r + b^2 r) for a line bundle (a, b) at r = A_X/A_E: scale-invariant, >= 8 pi^2 |a b| with equality at r = |a/b|
  (AM >= GM on the two legs = (anti-)self-duality).  Split background: E(r) = 4 pi^2 (A/r + B r), minimiser r* = sqrt(A/B), E* = 8 pi^2 sqrt(AB),
  excess over the topological floor = the Cauchy–Schwarz defect (zero iff every |a_i/b_i| coincides).
- S1 and S2 have the same energy function (A = 19, B = 11, r* = sqrt(209)/11); neither reaches the Hermitian–Yang–Mills floor of its
  class at any r (excess polynomials 41 r^2 - 14 r + 89 and 13 r^2 - 10 r + 13, both positive).  The floors: Delta = 14 with c_1 = (-5, -5)
  (S1, minimum 32/3 at r = 1) and Delta = 50 with c_1 = (7, 1) (S2, minimum 32/3 at r = 7): mirrors swap Delta and 2 |A_1 B_1|.  Neither
  minimum sits at a wall (r = 2, 1/2).  The M4 triangle is polystable at r = 3 and touches its floor there (control).
- Einstein frame: flux ~ E(r) Vol^-2, hyperbolic curvature ~ +48 pi M_8^6 r^-1/2 Vol^-3/2 — both positive, both decreasing in the volume: no
  volume stationary point (decompactification runaway, the surface analogue of `compactification.einstein_frame_potential`); the shape r is
  stabilised at fixed volume.  No cosmological constant is selected by the classical terms; delta^-6 e^{-2/alpha} does not appear (CC-34).

## 4. `mtft.parent` — candidates as data, gates as verdicts
`ParentCandidate` (dimension, group, content, supersymmetry flag, split stacks, family and Higgs blocks); `gates()` returns True/False/None
with witnesses for ledger rows 1-10 and section C; `ledger_table()` and `python -m mtft parent [ID]` print the scoreboard.  The three-stack
column reproduces the ledger's status column exactly.  `test_pv01.py` ties every number above to the code; the box-3 and box-4 scans are
slow tests (45 s, 120 s).

## Status
EXACT: the anomaly identities, the parity lemma, the energy functions, minimisers and floors, the scoreboard.  EXACT (finite): the scans within
their boxes.  OPEN: supersymmetry of the parent (undeclared — it decides whether the U(8) class needs a supersymmetric completion or only its
two Green–Schwarz fields), non-split stable bundles for the B1 classes, the Green–Schwarz construction, strong CP after the lift, all scales.
