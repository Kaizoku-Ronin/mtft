# Unreleased — after v0.34.0

## PG-01 corrections (2026-10-16): capture needs hypotheses (CC-67); the divisor series' second pole and summatory remainder (CC-68)

From Astra's review of `ca289dad` and Codex's RESPONSE-02 preflight, re-derived independently before filing; append-only.

- CC-67 caveats PG-5: dE/dt = -3H phi'^2 <= 0 holds for any H >= 0, but capture of the shape does not follow "for any H(t) > 0" — an integrable
  H leaves E(t) >= E(0) exp(-6 int H) (Gronwall; H = 0.1 e^-t keeps 77% in the toy run) and a fast-growing H freezes the field off the floor
  (Codex's exact solution). What holds, with proof: a convex coercive single well and H >= 0 of bounded variation with int H = infinity give
  capture — constant H and the radiation and matter eras (energies ~ t^-1.5 and t^-2); the hypotheses are sufficient, not necessary. New `vacuum_atlas.capture_conditions()`; the docstrings
  of the module, `toy_run` and `energy_argument` and the IV.7 answer now state the hypotheses.
- CC-68 caveats CC-65 and corrects `mtft.divisor_gas`, the Legend (`two_primon_energies`, `mean_weight_glaisher`) and the report: Z_D has a
  double pole at s = 0, -1/(2 s^2) - (ln 2 pi)/(2 s) + c0 with c0 = (zeta''(0) - gamma_1)/2 = -0.9667703052, besides the pole at s = 1; and
  sum_{n<=x} w_n has no O(1) remainder (it jumps by the unbounded w_n). What holds: -(1/2) log^2 x - log x - O(1) <= S(x) + zeta'(2) x <= 0,
  so S = -zeta'(2) x + O(log^2 x), with the s = 0 terms as its mean (proved through the first Riesz mean) and exact term by term under a smooth
  weight. New `divisor_gas.pole_at_zero()` (two routes), `summatory_bounds()`, `smoothed_sum()`.
- Ledger 282 claims (221 active, 36 caveat, 19 superseded, 6 retracted): PG-5 and CC-65 now caveat with notes, statements unchanged; CC-67 and
  CC-68 added (IV). Generated MATTER_LEDGER.md, MATTER_PUZZLE_MAP.md and docs/CAPABILITIES.md refreshed. Report: section 6 of
  `docs/SM/PG01_ENERGY_TWO_CLOCKS.md`.
- `tests/test_pg01_corrections.py` (8 tests); `tests/test_pg01_energy.py` updated for the two caveats and IV.7's claims. README's correction range
  now reads CC-36 to CC-68.

## PG-01 addendum (2026-10-14): the Wilson-line torus has two saddles (CC-66, PG-6); two visual pages

- `mtft.matter.vacuum_atlas.wilson_line_torus(eps)`: every critical point of the relative-U(1) Wilson-line potential on the whole torus R^2/Z^2
  (Newton from a grid of starts), typed by the eigenvalues of G^-1 H in the kinetic metric `kinetic_metric()` (143a1's Hodge metric, det 1), with
  the Morse count #min - #saddle + #max = chi(T^2) = 0 and the jet-vacuum curvatures in theta-coordinates and in the metric; `saddle_pairing()`
  guards the argument for why both saddles of the real spin structure sit at V = 0 exactly ((m, n) -> (-m - n, n) keeps |m + n tau| because
  Re tau = 1/2). Four critical points for each even spin structure, all at half-periods and nondegenerate: sink, source and two saddles (PG-6;
  exact in PA-01's product approximation that the half-periods are critical, numerical by two routes that there are no others). `atlas()` now
  lists the two saddles (nine rows, `saddles == 2`); `python -m mtft matter atlas` prints the census and both sets of curvatures.
- CC-66 caveats PG-4 (append-only; PG-4's statement is unchanged and now carries the note): its "No saddle appears anywhere" held only along PA-01's
  straight path, and conflicted when filed with BL-3, which had already found saddles at these half-periods for the B–L Wilson line. Also
  recorded: the Hessian eigenvalues 28.9 and 166.9 are theta-coordinate values; the Wilson-line masses follow the metric curvatures 37.394 and
  128.985 (ratio 1.857, not 2.404). PA1-1 and the rest of PG-4 stand.
- `tests/test_pg01_wilson_torus.py` (9 tests, two routes: the package's census on PA-01's box lattice, and in the test an independent disc sum
  with a winding-number count, a Lipschitz exclusion and a contraction bound for the census, and its own metric curvatures);
  `tests/test_pg01_energy.py` updated for the nine-row atlas, PG-6 and CC-66. Report: section 5 of `docs/SM/PG01_ENERGY_TWO_CLOCKS.md`.
- Puzzle map: IV.7 cites PG-6 and CC-66 and its answer names the saddles; piece IV's verdict no longer says "no saddle". Ledger 280 claims
  (221 active, 34 caveat, 19 superseded, 6 retracted); generated MATTER_LEDGER.md and MATTER_PUZZLE_MAP.md and docs/CAPABILITIES.md refreshed.
- Visuals: `viz/wilson_line_torus.html` (the force −G^-1 grad V on the holonomy torus, drawn in its own metric, three spin structures) and
  `viz/field_line_islands.html` (the standard map as the local model of a field-line Poincaré map, Ford circles, Greene's K_c); both self-contained,
  light and dark, checked in headless Chromium.

## PG-01 (2026-10-12): energy and the two clocks — the weights as the primon gas's energy, the vacuum as a dynamical system, the simulation gate (corrected the same day, CC-65)

- `mtft.divisor_gas` (five Legend entries, all EXACT): w_n = sum_{d|n} (log d)/d is the energy numerator of n's divisor gas at the Hagedorn
  temperature, w_n(beta) = -d/dbeta sigma_{-beta}(n) at beta = 1 (checked against the sieve); the marked ensemble's internal energy <log n>_D beside the
  unmarked gas's -zeta'/zeta = sum Lambda(n) n^-beta; the residue of Z_D at the wall is the mean weight -zeta'(2) = (pi^2/6)(12 ln A - gamma - ln 2 pi)
  (Glaisher); the alphabet at the negative integers W(-k) = zeta(-k)(ln A_{k-1} - B_k H_{k-1}/k) (Bendersky–Glaisher), which reads PV-02's
  W(-1) = -(ln 2 pi)/24; chi(SL(2,Z)) = -1/12, the rational Euler characteristic of the amalgam, equals zeta(-1) (Harder, K = Q) and the genus
  formula at 143 is orbifold Gauss–Bonnet (that |SL(2,Z)^ab| = 12 = 1/|chi| is a coincidence of size); the exact decomposition
  U_D = -zeta'/zeta - zeta''(beta+1)/zeta'(beta+1) (`marked_energy_decomposition`, `residue_at_first_zero`: the marked energy keeps the poles at the
  zeros); Bost–Connes quoted, no map to a physical temperature claimed. `python -m mtft.divisor_gas`.
- `mtft.matter.vacuum_atlas` (PG-01, 5 claims, piece IV): the critical-point atlas in the eight-types language — the shape modulus in a cosh well
  E/E* = cosh(phi - phi*) with floor at the self-dual point and curvature exactly 1 (PG-1, a sink); the volume with no finite critical point when both
  Einstein-frame coefficients are positive (runaway, N22); the condensate's non-hyperbolic flat direction lifted to exactly 4 m^2 by a mass term with
  0 < m^2 < xi; the Wilson line's 2D sink at the jet vacuum and 2D source at the untwisted point (full Hessians, all three spin structures, exact in
  PA-01's product approximation); adjoints sinks, A_X no critical point; no saddle anywhere (PG-4). Hubble friction: spiral iff m > 3H/2; capture of
  the shape guaranteed by the energy argument dE/dt = -3H phi'^2 <= 0 (`energy_argument`), a DIAGNOSTIC toy run illustrating it beside the volume's
  runaway (log Vol ~ (2/3) log t in the toy); the simulation gate (PG-5; sub-puzzle I.7: I.1, I.6/GS-1, GC-01, a cosmological background, then
  standard physics). `python -m mtft matter atlas [--json]`; figure `scripts/pg01/` -> `docs/SM/figs/pg01_energy_clocks.png`; report
  `docs/SM/PG01_ENERGY_TWO_CLOCKS.md`; `tests/test_pg01_energy.py` (9 tests).
- CC-65 (same day, after independent review, which confirmed every number of the wave): the first filing's sentence that marking by w_n moves the
  energy off the zeros of zeta was false (the marked energy keeps every pole of the von Mangoldt series and adds poles at the shifted zeros of zeta';
  RH-blindness belongs to Z_D and to the summatory function, whose error is O(log^2 x/x)); "the twelve in front is zeta(-1)" tightened to
  -1/zeta(-1) = 2 pi^2/zeta(2); the Bendersky identity named as a definition rearranged; the Wall–Serre/Harder attribution and the coincidence-of-size
  reading of |SL(2,Z)^ab| = 12; the Bost–Connes paragraph now says what the marked gas is and is not. PG-2 and PG-3 superseded by PG-4 and PG-5
  (three results stated without their hypotheses; a run-length-dependent crossing count presented as a result). Append-only, as CC-63 and CC-64.
- Puzzle map: wave PG-01 on piece IV (IV.7 answered, claims PG-1, PG-4, PG-5, CC-65), I.6's settle text notes the runaway, I.7 the gate (open; PG-5),
  N22 (PG-4). Ledger 278 claims (220 active, 33 caveat, 19 superseded, 6 retracted); 73 sub-puzzles, 29 open.

## GF-01 (2026-10-11): forces first — the eight dots, the sealed touch and charge-quantisation rules, the parent columns, the clock of twelve

- `mtft.matter.forces` (GF-01, 9 claims): a rank-8 parent as eight dots; grouping the dots chooses the forces (carriers = circles + lines inside)
  and the lines between groups are the particles — difference lines (U(8) adjoint), sum lines (Lambda^2(8), PV-01) and sheets (the spinor of SO(16),
  E8). The Standard-Model fields are the difference lines of U(3) x U(1) x U(2) x U(1)^2 with PS-01's hypercharges and the four Yukawa couplings are
  the four triangles of dots whose third side is the Higgs (H_u on (L, A2), H_d on (L, B), as DS-01). Reproduces PS-01, DS-01 and PV-01's 28 blocks.
  Sealed rules: PR-S3, DS-01's touch criterion on the three supersymmetric block menus (mixed: U(8) FAIL — only the nu^c line, DS-1 restated;
  SO(16) and E8 PASS through sum lines and sheets; PV-01's non-supersymmetric U(8) has the sum lines as fermions only and was not scanned); PR-S5,
  charge quantisation through the Z6 of the Standard Model (PASS on every line; holds for every odd T_3R coefficient, so necessary and not what
  fixes hypercharge; first sealed as PR-S4, whose statement read the global form backwards and was superseded the same day, CC-64). Hypercharge by two routes: anomalies on the stack U(1)s leave <T_3R, B - L>, Y(Q) = 1/6
  fixes B - L and a neutral Higgs component fixes T_3R to +-1 (the u <-> d relabelling); the Pati–Salam route agrees. Couplings as the shape of the
  force torus (1/e^2 = 1/g^2 + 1/g'^2 exact, theta_W the diagonal's tilt; GC-01's matching is Pythagoras on the Pati–Salam torus). The three parent
  columns of I.1 with PV-02's twist: U(8) recovers PV-02's three conditions and six doublet pairs; SO(16) and E8 through U(8) are closed — the sum
  line e_c + e_L makes SO(16)'s quark-doublet count 2 phi_c, E8's three quark doublets come with three chiral exotic (3, 2)_{-5/6} doublets, and in
  both the lepton rows force net(L) + net(e^c) = 0; inconsistent at every hypercharge shift (the lattice (1/48) Z on [-2, 2] contains every special
  shift, and a generic shift fails too). So within the supersymmetric menus the blocks that reach the down sector live in columns that do not give
  the exact Standard Model: the wall of piece IX is parent-level there.
  The clock of twelve: SL(2,Z)^ab = Z/12 (S order 4, ST order 6), G2's twelve directions in the colour plane, twelve carriers; the golden loop (trace 3,
  4 log phi, through i) and the twelve loop (trace 4, 2 log(2 + sqrt 3), through rho) cross the mirrors only at CM points; the golden loop is reversed
  by S, the twelve loop only by a mirror; X(5) is the icosahedron. The systole of X0(143) is the twelve loop because (5/13) = -1 — logged as trial
  X17 (selects 118, 141, 142, 143; tried after looking). CLI `python -m mtft matter forces [--touch P] [--column P] [--box N]`; figures
  `scripts/gf01/` -> `docs/SM/figs/gf01_*.png`; report `docs/SM/GF01_FORCES_FIRST.md`; `tests/test_gf01_forces.py` (16 tests).
- CC-64 (same day, after independent review, which confirmed every number from its own construction of the SO(16) and E8 root systems): GF-2, GF-3,
  GF-7 and GF-8 superseded by GF-13, GF-10, GF-11 and GF-12, and seal PR-S4 by PR-S5 — scope of the touch verdict, direction of the global-form
  reading, the E8 wording, the twelve loop's second self-crossing and the turn angle at rho.
- Puzzle map: wave GF-01 on piece V (V.7 charge quantisation and hypercharge, V.8 which parent's lines carry three families — both answered);
  GF-2 on IX (IX.10 answered: no consistent parent reaches the down sector at first order), GF-5 on VI (VI.6 answered: the couplings as a shape),
  GF-8 on II (II.11 answered: the clock and where phi lives), GF-9 on IV (IV.6 cites X17); I.1's settle text is now a two-way choice (GF-6, GF-7) and
  N21 records the closure. Pre-registered tests T7 (PR-S5, pass) and T8 (PR-S3, mixed). Ledger 272 claims (216 active, 33 caveat, 17
  superseded, 6 retracted); 71 sub-puzzles, 28 open. Piece I still carries no ledger claim of its own.

## TR-01 and PR-01 (2026-10-10): time reversal, sealed pre-registrations, the trials counter

- `mtft.matter.time_reversal` (TR-01, 7 claims): time reversal tau -> -conj(tau) on a gauge coupling has two mirror classes under duality,
  theta = 0 (Farey edges, fixed form of content 2) and theta = pi (Farey medians, content 1), the two conjugacy classes of reflections of
  PGL(2, Z); two exact routes agree on all 212 reflections with entries <= 12. The oloid and triangle mirrors and 143a1's period lattice are of
  the theta = pi class. The thirteen star quotients of SQ-01 are identified as curves (the curve of modular degree 4 in each class; 65a2, not the
  optimal 65a1): nine with one real circle, four with two. The root number is the t -> -t parity of Lambda (SQ-2 restated). CM points on a
  mirror are the ambiguous classes, read off the reduced form; W_143's four real fixed points are two of each class. A sealed strong-CP test for
  I.4 (TR-8, seal PR-S2, puzzle test T6, outcome pending; the first seal PR-S1/TR-6 rejected theta = pi points outright and was superseded the
  same day after independent review, CC-63, which also corrected PR-01's trial X12 and TR-5's count of W_13's fixed points). CLI `python -m mtft matter timerev`; provenance `scripts/time_reversal_curves.py`;
  figure `docs/SM/figs/tr01_mirror_classes.png`.
- `mtft.matter.prereg` (PR-01, 3 claims): rules sealed by SHA-256 of statement, code, parameters and window; `verify()` recomputes from the code
  as it stands, `evaluate()` refuses a changed rule; T0–T5 listed as unsealed legacy. The trials counter logs every property tried as a reason
  for 143 on the thirteen torus levels (16, none pre-registered, 4 single it out, all post hoc) with the chance baseline that applies only to
  sealed rules. The first correction, CC-63, is the registry working as designed: the superseded seal stays in the record with its hash. CLI `python -m mtft matter prereg [--verify]`. The test suite refuses an unlogged distinction.
- Puzzle map: I.4 cites TR-8 with pre-registered test T6 (outcome "pending", a new outcome value); II.10 answered (the thirteen curves and their
  time-reversal classes); IV.6 cites PR-3. Ledger 258 claims (CC-63 included). Piece I still carries no ledger claim of its own.

## v0.34.0 release-review corrections (reconciled 2026-10-01)

- Search scope: `anomaly8d.b1_scan` imposes `a_c*b_c = 0` and two individual mixed antitriplet blocks of index -3. Its 60/156
  assignments and the finite trilemma are statements about that imposed subclass, not an exhaustive B1 classification; an explicit
  U(8) assignment with `a_c*b_c = 1` has exactly `SM_NET(3)` plus 63 neutral singlet indices (`tests/test_parent_edge_cases_review.py`).
- `parent.gates`: unsupported content (e.g. `u1(q)`) makes the 8D anomaly verdict OPEN with the content named, instead of being dropped;
  identically zero slope conditions are reported as the whole positive-r locus.
- `surface.blowup`: equal-stack slopes give the whole ample cone in `polystable_locus`; `collinear_family(a_L)` has its D-flat locus at
  eps/A_E = 1/3 and ampleness r > 1/3 for every nonzero a_L (the released docstring's a_L/3 scaling was wrong), nine doublet pairs at a_L = 1.
- `anomaly8d.sm_net_content` keeps weak-triplet labels.
- `product_surface.s1_colour_recombination_dims`: gluing dimension 32 - 15*pairing_rank; independence alone does not fix the rank;
  the historical 17, 2, 2 values are the generic counts (`tests/test_s1_recombination_rank.py`).
- CC-35 wording: mu_N(y) is a defined arithmetic stiffness diagnostic, NOT the Hessian of the implemented single-trace Polyakov action
  (at the SU(3) centre-symmetric point, kappa = 1, y = 0.18174, the tangent gradient is 0.01928985 and the actual Hessian eigenvalues are
  3.60e-8 and 0.0984541; at y = 0.3 one is -0.00242941). The CC-35 lemma has one exception: SU(2), n = 2 (Tr U^2 = chi_1 - 1).
- Release navigation: `CHANGELOG_v0340.md` added (the shipped notes, verbatim); README points at the v0.34.x line; MANIFEST.in now ships
  `viz/*.html`, `viz/*.png`, `docs/**/*.png` and `scripts/cosmology/*.py` — the cause of the three artifacts missing from the v0.34.0
  sdist (they were excluded by the manifest). `docs/SM/figs/cc35_polyakov_kernel.png` and `viz/arithmetic_susy_picture.html` restored
  from the PV/CC-35 handoff; `scripts/cosmology/cc34_cosmology_audit.py` must come from the repository (not in any sdist).
- Record: `docs/SM/RECONCILIATION_V0340_REVIEW_MATTER03.md`.

## Matter-01 (2026-10-01): U(1) reference engine and conditional bridge

- New `mtft.matter` CLI and massive explicit-spinor tree-level QED amplitude
  for electron-positron annihilation into two photons; full-sphere event
  factor, Ward/Bose/gauge-shift/boost checks and independent scalar oracle.
- Frozen PDG reference mass/coupling inputs with explicit override provenance.
  These are not MTFT predictions; no quantum-state time evolution is added.
- Conditional same-charge Gram normalization, electroweak matching and uniform
  Einstein-frame breathing gauge prefactor; C3X and B3 remain distinct.
- Exact B3 flux-class audit and RR cohomology lower bounds. The archived
  minimum-index-extras witness can have much larger cohomology bounds than
  the simpler collinear control. Actual spectra, normalized Yukawas, exotic
  masses and the photon axion/Stueckelberg kernel remain OPEN.
- Offline interactive lab, scientific overview, atlas evidence ledger without
  mapping promotions, focused runner and [report](../SM/MATTER01_REPORT.md).
- Existing results and release version are preserved. Full release tests and
  long archived studies are separate from this milestone gate.

### Matter-02 collinear representative and holomorphic coupling

- Specify the cusp-at-infinity blow-up and trivial elliptic flat factor.
- Add exact cohomology dimensions, rational cusp-jet certificates and a
  globally valid rank-three lepton–Higgs holomorphic residue witness.
- Append the opposite archived/BHV chirality conventions without changing
  released index functions or source data.
- Audit hypercharge traces and explicitly INPUT axion/Higgs mass kernels;
  physical photon and normalized particle masses remain open.
- Add the local interactive mode lab, scientific overview, independent
  certificates and a focused reproducibility runner.

### Matter-03 Higgs branch, deformed light modes and kinetic inputs

- Derive the arbitrary-Hermitian-metric obstruction to a nonzero single
  directed Higgs on the unchanged full B3 wall background.
- Certify an explicit reverse field and odd exceptional spectral valuation
  seven; canonical-twisted surface stability supplies classical paired
  F/D-solution existence. No physical metric or Higgs scale is computed.
- Add a charged deformation-complex census: six aa ranks total 41 and leave
  237 degree-one classes per sign for the fixed pair; all aligned pairs in
  these spaces leave at least 160. Observed family selection remains open.
- Add lower-index field-coordinate kinetic normalization with explicit Grams,
  canonical Higgs direction, vev and parent factor, including complex-basis
  covariance tests. Historical section normalization remains unchanged.
- Add derivations, source APIs, focused runner, offline lab and scientific
  overview. Standard Model metastability is kept distinct from this branch.
- Preserve v0.34.0 release pins, original arithmetic data and prior archives.

### Matter-04 the adjoint question and spin structures (2026-10-02)

- `mtft.matter.spin_twists`: the four spin structures of 143a1 (2-division polynomial 4x^3 - 4x^2 - 4x - 7, discriminant 16 Delta, S3 over Q(sqrt(-11)));
  adjoint fermion zero modes by twist — (1, 14, 13) for BHV, (0, 12, 12) for a curve twist, none for any even spin structure on 143a1 (Scherk–Schwarz);
  charged fermion blocks reduced to exceptional jets (|Q| = 1 split census 279 -> 153 per sign); flat-circle one-loop sign; QCD b0 (-78 with the 27 adjoints,
  7 below the Scherk–Schwarz scale); the flux minimum under the twist (288 against 1152 for the Matter-02 collinear choice; pairwise coinciding classes).
- Report `docs/SM/MATTER04_ADJOINT_SPINORS.md`; `tests/test_matter04_spin_twists.py` (5 tests).

### Matter-05 one-loop adjoint masses; Pati–Salam first look (2026-10-02)

- `mtft.matter.one_loop_adjoints`: sign theorem for the 27 bosonic adjoints under the Scherk–Schwarz twist (positive one-loop mass^2 from gauge loops; Phi and
  torus Wilson line positive with matter loops; curve Wilson lines' matter-loop sign open), closed forms at leading order, 143a1 lattice sums for the three even
  spin structures (split 1 + 2: S = 2.107219 versus 3.238510), the systole of X0(143) (trace 4), QCD windows (b0 = 7, -6, -20).
- `mtft.matter.pati_salam`: the Matter-04 minimum regrouped as U(4) x U(2) x U(2) — three net families each of (4,2,1) and (4bar,1,2) with nu^c, a (1,2,2) block;
  index counts only.
- Report `docs/SM/MATTER05_ONE_LOOP_ADJOINTS.md`; `tests/test_matter05_one_loop_adjoints.py` (4 tests).

### PS-01 Pati–Salam breaking, anomalies and the b–tau check (2026-10-02)

- `mtft.matter.ps_breaking`: Pati–Salam hypercharge and branching from the stacks (= PV-02); exact anomaly table (nonabelian and Witten cancel; both
  physical stack U(1)s anomalous, needing Green–Schwarz masses); Wilson lines on 143a1 (rank-preserving; Higgs pair iff N_L = N_B); the Kahler FI drives the
  (1,2,2) bidoublet (electroweak), the gamma direction with a Pati–Salam-symmetric gamma-flux drives the (4,1,2) condensate (Pati–Salam breaking); gamma-flux
  is chirality-neutral with 4c^2 localised modes, bosonic ones at phi^-1(E[c]) (the four cusps for c = 1); SM one-loop b–tau equality at 1.1 x 10^6 GeV,
  Georgi–Jarlskog light-generation ratios, seesaw estimate.
- Report `docs/SM/PS01_PATI_SALAM_BREAKING.md`; `tests/test_ps01_breaking.py` (6 tests).

### NU-01 the neutrino gate (2026-10-02)

- `mtft.matter.neutrino_gate`: the cheaper Pati–Salam-breaking gamma-flux (1, 1, -1, -1, -1) (D-flat, bidoublet untouched, 256 against 384 localised extras);
  the neutral-fermion inventory (net 3 nu_L and nu^c; 19 vector-like (nu^c, N) pairs per entry; no B-L-neutral singlets); disjoint bosonic (phi^-1(E[c])) and
  fermionic (phi^-1(cQ = P_eta)) localisation and the resulting holomorphic selection rule; the Atkin–Lehner regular representation on the four cusp-localised
  modes (one in the (-,-) sector); B-L selection rules; Dirac masses at M_PS and the seesaw requirements (type I 2.9 x 10^14 GeV; inverse mu(3) = 4.3 MeV).
- Report `docs/SM/NU01_NEUTRINO_GATE.md`; `tests/test_nu01_neutrino_gate.py` (3 tests).

### NU-02 how small can mu be (2026-10-02)

- `mtft.matter.neutrino_mu`: the 143a newform from the 2-division polynomial, the modular parametrisation with an integrality check of its periods (Manin
  constant 1) and phi o W_143 = phi; the four points over each nonzero 2-torsion point (certified by j-invariant pairs), free Atkin–Lehner orbits (regular
  representation), depth (3.589 / 3.397) and stretch (0.1965 / 0.2049); the bulk suppression mu/M ~ exp(-beta sqrt r) with beta 0.59 - 1.0 and the shape ratios
  it requires; a D-flat flux with c = 1 on nu^c.
- Report `docs/SM/NU02_NEUTRINO_MU.md`; `tests/test_nu02_neutrino_mu.py` (4 tests).

### AL-01 Atkin–Lehner symmetry and one family per cusp (2026-10-02)

- `mtft.matter.cusp_families`: lattice of multi-point blow-ups at cusps; exhaustive check that full Atkin–Lehner symmetry (all sixteen gauge actions) forbids
  three families; the two-cusp left–right-parity background (216 zero modes, tree-level Yukawa zero); the one-family-per-cusp background on three cusps
  (exactly three families with no vector-like quarks or leptons, 96 zero modes, D-flat, ample iff A_X > A_E, unique up to Atkin–Lehner, tree-level Yukawas
  diagonal in the cusp basis so CKM = 1 at tree level).
- Report `docs/SM/AL01_CUSP_FAMILIES.md`; `tests/test_al01_cusp_families.py` (4 tests).

### YK-01 Yukawa couplings of jet families; CC-36 (2026-10-02)

- `mtft.matter.yukawa_jets`: tree-level Yukawas between exceptional-jet fermions vanish (cup product of jet classes is zero); classification of Z/2 Wilson lines —
  tree-level Yukawas exist iff the quark–lepton group carries the Wilson line equal to the spin structure; the resulting Standard-Model twist pattern (matter
  fermions, Higgs bosons and gauge bosons keep their zero modes; sfermions, Higgsinos, gauginos are jets or removed); its price on the three-cusp background (896
  against 96 fermion zero modes).
- CC-36 (`docs/SM/CC36_JET_YUKAWA_RETRACTION.md`): retracts AL-01's tree-level diagonal-Yukawa statement and the tree-level origin assumed in PS-01 / NU-01 (the SU(4)
  relations themselves stand).
- Report `docs/SM/YK01_YUKAWA_JETS.md`; `tests/test_yk01_yukawa_jets.py` (4 tests).

### PA-01 Path A: the matter-untwisted point is a one-loop maximum (2026-10-02)

- `mtft.matter.path_a`: origin-by-origin content of a matter entry (25 vector-like pairs + net 3); the sources a vector-like mass would need (stack-U(1) moduli
  vevs); the one-loop potential of the relative-U(1) Wilson line, which only matter feels — the untwisting point is its maximum (Hessian exactly opposite to the
  jet vacuum's), so the background rolls back to the jet vacuum; there the 25 pairs become messengers at the Scherk–Schwarz scale (M sqrt(A_E) = 3.106 / 3.537),
  so Path A converges on Path B (Yukawas from messenger exchange, ~ v_R / M).
- Report `docs/SM/PA01_PATH_A.md`; `tests/test_pa01_path_a.py` (3 tests).

### PB-01 messenger-exchange Yukawas (2026-10-02)

- `mtft.matter.messenger_yukawa`: SU(n) invariant counts by Weyl integration; the holomorphic messenger route closes (no single condensate insertion, Sigma^n
  invariants vanish for identical bosons, the D-term vacuum keeps one orientation); the non-holomorphic route through Higgs–heavy-bidoublet mixing driven by the
  residual nonabelian D-term; the cusp Bergman kernel of weight-2 forms as a Klein-group circulant with sector weights (2.485, 4.438, 4.502, 1.066).
- Report `docs/SM/PB01_MESSENGER_YUKAWA.md`; `tests/test_pb01_messenger_yukawa.py` (3 tests).

### CK-01 quark mixing and the cusps; CC-37 (2026-10-02)

- `mtft.matter.ckm_cusps`: exact horoball distances between the cusps of X0(143) (log 143, log(143^2/13), log(143^2/11)); the 1/143 bound on mixing between families at
  different cusps, which excludes every split assignment (CC-37: AL-01's one-family-per-cusp background retired; the Matter-04 single point is canonical); the family
  triplet as the degree-2 jets x^2, x w, w^2 with torus charges 0, 1, 2; the rank-one leading Yukawa; the Froggatt–Nielsen texture Y_{qq'} ~ eps^{q+q'} with a
  Monte-Carlo comparison to the quark masses and the CKM angles.
- Report `docs/SM/CK01_CKM_CUSPS.md`; `tests/test_ck01_ckm_cusps.py` (3 tests).

### EPS-01 and the matter-programme toolkit (2026-10-02)

- `mtft.matter.epsilon` (EPS-01): exact selection rules at the blown-up point — torus-constant D-term gives rank one; a single-band gamma-localised condensate keeps
  it; family rotations keep it; Y_qq' follows the D-term's torus charge, so epsilon is the condensate's misalignment; a rank-one aligned condensate gives
  eps_u = eps_d (data need 2.30); the gamma-flux estimate eps ~ (2.1-6.7)|c|/sqrt(r) with the Arakelov cusp density 0.02434 from the frozen periods.
- `mtft.matter.background`: a `Background` object and seven presets assembling indices, zero modes, twists, D-flatness, ampleness, anomalies, the Yukawa verdict,
  Wilson-line stability, messenger scale, CKM admissibility, family structure and gamma-localisation; `ps_symmetric_gamma_tradeoff` (TK-1: a global Higgs forces
  gamma-charge 2 c_4 on nu^c).
- `mtft.matter.ledger`: 37 entries Matter-04 ... EPS-01 with class, status, module and test; `check()`; generated `docs/SM/MATTER_LEDGER.md`.
- Command line: `python -m mtft matter presets | background NAME | ledger | epsilon` (and `python -m mtft.matter ...`).
- CK-01 qualification recorded (rank one assumes a torus-constant D-term). Docs: `docs/SM/EPS01_EPSILON.md`, `docs/SM/MATTER_PROGRAM_MAP.md`.
- Tests: `tests/test_eps01_epsilon.py` (3), `tests/test_toolkit_matter_program.py` (4).

### FS-01 SU(2)_R-breaking flux scan; CC-38 (2026-10-02)

- `mtft.matter.flux_scan`: stack-level backgrounds (curve degrees and gamma-charges per stack) with exact Kahler solving for D-flat Pati–Salam breaking by one
  Standard-Model singlet; the scan; curve degrees alone suffice (CC-38 supersedes PS-01's gamma necessity); TK-1 lifted; a gamma-localised Higgs sets its sector's
  epsilon; the mod-4 Higgs-charge rule with SU(4) kept (eps_d/eps_u = 3 at the smallest charges); four background families; ratio 3 with a charged nu^c needs
  charges (2, 6).  Command line: `python -m mtft matter scan`.
- Report `docs/SM/FS01_FLUX_SCAN.md`; `tests/test_fs01_flux_scan.py` (4 tests); ledger 43 entries (`docs/SM/MATTER_LEDGER.md` regenerated).

### RT-01 the ratio-3 background (2026-10-02)

- `mtft.matter.ratio_three`: leading eps_d/eps_u = 3 (the Higgs entries differ only by gamma-charge); the finite-width correction gives 2.84-2.96 where the eps
  magnitudes fit (25 per cent above 2.30); jet-sector neutrinos in the one-orientation vacuum behave as a type-I seesaw at v_R (27 MeV, independent of the
  D-term suppression).  Command `python -m mtft matter ratio3`.  Report `docs/SM/RT01_RATIO_THREE.md`; `tests/test_rt01_ratio_three.py` (2 tests); ledger 46 entries.

### LEP-01 the charged-lepton texture; CC-39 (2026-10-02)

- `mtft.matter.leptons`: Georgi–Jarlskog quantified at the b–tau scale; Y_e = Y_d exact in every SU(4)-preserving background; the Georgi–Jarlskog-friendly order of
  SU(4)-breaking gamma-flux; equal-overlap alignment fits (quarks charge-0 dominated, leptons charge-1 dominated); Koide values and a recorded numerology example.
  Command `python -m mtft matter leptons`.  CC-39 downgrades RT-1 and FS-4 to caveats.  Report `docs/SM/LEP01_LEPTONS.md`; `tests/test_lep01_leptons.py` (3 tests);
  ledger 52 entries.

### CP-01 spin-1 phases of the condensate alignment (2026-10-03)

- `mtft.matter.condensate_phases`: degree-2 jets as spin 1 of the local SU(2); the two quartic invariants and the coherent/polar phases; the tree-level D-term quartic
  (1/5)|c|^4 - (1/15)|theta|^2 prefers polar; textures of each phase; the lepton sign-flip test (m_tau = m_b, m_mu m_e = m_s m_d; 3.57 and 0.282 against 3.06 and 0.295).
  Command `python -m mtft matter phases`.  Sector phases of the equal-overlap fits (quarks coherent, leptons polar).  Report `docs/SM/CP01_CONDENSATE_PHASES.md`; `tests/test_cp01_condensate_phases.py` (4 tests); ledger 57 entries.

### CP-02 zero flips, parity, and CC-40 (2026-10-03)

- `mtft.matter.condensate_phases`: `zero_flip`, `rho_branches`, `landau_ground`, `geometric_test`.  CC-40 retracts CP-5 (zero-flip artefact); CP-3 and CP-4 carry
  branch caveats.  New ledger entries CP-6 to CP-9 (62 entries); `tests/test_cp01_condensate_phases.py` now 6 tests.

### CP-03 the parity-odd term as a transverse field (2026-10-03)

- `mtft.matter.condensate_phases`: `zeeman_ground`, `texture_invariants`, `LAMBDA_CRITICAL`, `rho_departures`.  CP-8 now a caveat (necessary, not sufficient).
  Ledger entries CP-10 to CP-14 (67 entries); `tests/test_cp01_condensate_phases.py` now 8 tests.

### XQ-01 the overlaps X_qq' (2026-10-03)

- `mtft.matter.overlaps`: cusp Bergman jets of X_0(143) from the period-normalised q-expansions, correlations, charge-conserving weights, the family-normalisation
  convention and the data's rho targets.  CC-41 downgrades CP-12's targets.  Report `docs/SM/XQ01_OVERLAPS.md`; `tests/test_xq01_overlaps.py` (3 tests); ledger 72 entries.

### HJ-01 the family vertex of the jet classes (2026-10-04)

- `mtft.matter.jet_vertex`: the curve-charge selection rule, the round-blow-up texture (2 : 1 : 1), the 143a1 and cusp-peaked Higgs profiles, the Taylor texture.
  CC-42 retracts XQ-4; XQ-1 is a caveat.  Report `docs/SM/HJ01_JET_VERTEX.md`; `tests/test_hj01_jet_vertex.py` (3 tests); ledger 77 entries.

### HP-01 the light Higgs of the ratio-3 background (2026-10-04)

- `mtft.matter.higgs_profile`: gamma-localisation over phi^-1(E[Delta]), the graph Taylor data of 143a1, the A-type verdict, the superpotential rule, the coupling
  Higgs's exceptional charge.  CC-43 supersedes FS-4 and RT-1; RT-2 and FS-3 become caveats.  Report `docs/SM/HP01_HIGGS_PROFILE.md`;
  `tests/test_hp01_higgs_profile.py` (3 tests); ledger 82 entries.

### YA-01 character audit of the Yukawa amplitude (2026-10-04)

- `mtft.matter.yukawa_audit`: overall-U(1) characters, the holomorphic-Higgs obstruction, the scale ratio, the curve-charge-only entries and texture.  CC-44
  caveats PB1-2, EPS1-2, HJ-2, HJ-4.  Report `docs/SM/YA01_YUKAWA_AUDIT.md`; `tests/test_ya01_yukawa_audit.py` (3 tests); ledger 87 entries.

### LB-01 lattice breaking at the blow-up scale (2026-10-04)

- `mtft.matter.lattice_breaking`: 143a1's lattice multipoles, the mod-2 parity selection, the six-unit breaking at the blow-up scale, the jet-Yukawa bound and the
  top Yukawa at the b–tau scale.  Report `docs/SM/LB01_LATTICE_BREAKING.md`; `tests/test_lb01_lattice_breaking.py` (3 tests); ledger 91 entries.

### TF-01 can any jet configuration host the top? (2026-10-05)

- `mtft.matter.top_family`: the minimal-breaking theorem u = |e_Q| + 1, the two-fold anisotropy of 143a1's lattice, the two-plus-one flux-family design.  CK1-2
  becomes a caveat.  Report `docs/SM/TF01_TOP_FAMILY.md`; `tests/test_tf01_top_family.py` (3 tests); ledger 94 entries.

### TP-01 the flux family fails, the constant-jet top survives (2026-10-05)

- `mtft.matter.top_point`: flux chirality 24 b, the global-condensate cost, the two-point jet counts, the unique-Higgs count, the Q-block D-flat relation.  CC-45
  retracts TF-3; CC-46 caveats TF-2 and LB-4.  Report `docs/SM/TP01_TOP_POINT.md`; `tests/test_tp01_top_point.py` (3 tests); ledger 100 entries.

### TA-01 auditing the two-point design (2026-10-05)

- `mtft.matter.two_point_audit`: fibre geometry of 143a1, the blow-up constraints, jet-tail mixing and the required separation.  TP-3 becomes a caveat.
  Report `docs/SM/TA01_TWO_POINT_AUDIT.md`; `tests/test_ta01_two_point_audit.py` (3 tests); ledger 104 entries.

### SP-01 where the top can live (2026-10-05)

- `mtft.matter.special_points`: canonical-metric path lengths on X_0(143), the curve size fixed by V_cb, cusp-to-cusp mixing.  CC-47 caveats CK1-1 and CC-37.
  Report `docs/SM/SP01_SPECIAL_POINTS.md`; `tests/test_sp01_special_points.py` (2 tests); ledger 109 entries.

### GC-01 what the gauge couplings fix (2026-10-05)

- `mtft.matter.gauge_couplings`: the shape/volume statement, the Pati–Salam matching ratio, quantised flux ratios and their scales, b–tau at those scales.
  CC-48 caveats PS1-4 as the Pati–Salam scale.  Report `docs/SM/GC01_GAUGE_COUPLINGS.md`; `tests/test_gc01_gauge_couplings.py` (3 tests); ledger 115 entries.

### TH-01 thresholds and the scale tension (2026-10-05)

- `mtft.matter.thresholds`: vector-like pair beta coefficients and the left–right cancellation, the required asymmetry, adjoint route, pair counts per design,
  the DKL factor of 143a1.  Report `docs/SM/TH01_THRESHOLDS.md`; `tests/test_th01_thresholds.py` (3 tests); ledger 119 entries.

### VM-01 what gives the vector-like pairs their masses (2026-10-05)

- `mtft.matter.vector_like_masses`: SU(4) operator selection, vacuum sources, the tree-level verdict, the cutoff ratio a conjugate condensate would imply.  TH-3
  becomes a caveat.  Report `docs/SM/VM01_VECTOR_LIKE_MASSES.md`; `tests/test_vm01_vector_like_masses.py` (2 tests); ledger 123 entries.

### CO-01 is there a conjugate condensate? (2026-10-05)

- `mtft.matter.conjugate_orientation`: the closed-form D-term vacuum (8/13), the degenerate conjugate family and its lifting, the right-handed scalar masses,
  the gauge-loop asymmetry.  VM-4 superseded.  Report `docs/SM/CO01_CONJUGATE_ORIENTATION.md`; `tests/test_co01_conjugate_orientation.py` (2 tests); ledger 127 entries.

### VS-01 the condensate's feedback on the (4,2,1) scalars (2026-10-05)

- `mtft.matter.vacuum_stability`: squark and slepton masses after condensation, the stability bound xi_4L >= (11/13) xi_R4, the curve-degree ordering, the
  radiative asymmetry scaling.  CC-49 caveats TP-4 and FS-1.  Report `docs/SM/VS01_VACUUM_STABILITY.md`; `tests/test_vs01_vacuum_stability.py` (2 tests); ledger 132 entries.

### GS-01 Green–Schwarz relaxation, the slepton direction, the vector-like symmetry (2026-10-05)

- `mtft.matter.green_schwarz`: the relaxed D-term minimum, the irreducible SU(4) slepton shift, U(1)_V charges and invariance.  PS1-3 becomes a caveat.
  Report `docs/SM/GS01_GREEN_SCHWARZ.md`; `tests/test_gs01_green_schwarz.py` (2 tests); ledger 135 entries.

### MS-01 the Casimir energy and the curve's size (2026-10-05)

- `mtft.matter.moduli_casimir`: odd-winding lattice sums (reproducing Matter-05), the minimum condition and location, the heat-kernel sign count.
  Report `docs/SM/MS01_MODULI_CASIMIR.md`; `tests/test_ms01_moduli_casimir.py` (2 tests); ledger 139 entries.

### MS-02 the full twisted spectrum (2026-10-05)

- `mtft.matter.moduli_casimir`: Dolbeault-pair curvature coefficients, the per-block a_1/N, the total-curvature term, the no-minimum verdict.  MS-4 superseded.
  `tests/test_ms01_moduli_casimir.py` now 3 tests; ledger 142 entries.

### CS-01 pre-registered canonical-shape test (2026-10-05)

- `mtft.matter.canonical_shape`: the pre-registered window and rule, four candidate families, the chance level.  No hit.  Report `docs/SM/CS01_CANONICAL_SHAPE.md`;
  `tests/test_cs01_canonical_shape.py` (2 tests); ledger 145 entries.

### DD-01 a parity selection rule for the CKM matrix (2026-10-05)

- `mtft.matter.double_duty`: parities under the torus involution, the decoupled-family theorem, the required breaking and coefficient, the refined separation.
  Report `docs/SM/DD01_DOUBLE_DUTY.md`; `tests/test_dd01_double_duty.py` (2 tests); ledger 149 entries.

### CA-01 the Cabibbo angle and two measurements of the curve's size (2026-10-05)

- `mtft.matter.cabibbo`: r from V_us via the gamma-tilt, r from V_cb per home, their gap, the post-hoc Fricke candidate, the V_ub bracket, the c = 1 vector-like scale.
  Report `docs/SM/CA01_CABIBBO.md`; `tests/test_ca01_cabibbo.py` (2 tests); ledger 154 entries.

### TO-01 a pre-registered third observable (2026-10-05)

- `mtft.matter.third_observable`: the induced-second-family bound, the V_cb range from mass ratios, the left–right asymmetry needed.  CC-50 caveats SP-2, DD-4, CA-2
  and supersedes CA-3.  Report `docs/SM/TO01_THIRD_OBSERVABLE.md`; `tests/test_to01_third_observable.py` (2 tests); ledger 158 entries.

### FL-01 which tail sets which angle (2026-10-05)

- `mtft.matter.flux_tails`: the angle–tail identification, the tail ordering the data demand, the curve-flux bound.  CC-51 corrects TO-01 item 4.
  Report `docs/SM/FL01_FLUX_TAILS.md`; `tests/test_fl01_flux_tails.py` (2 tests); ledger 163 entries.

### CD-01 the condensate's u^c-but-not-d^c dressing (2026-10-05)

- `mtft.matter.condensate_dressing`: one-step gauge selection, the dressing factor z the data need, the nearly symmetric tails at z = 0.25.  FL-2 caveat, FL-4 superseded.
  Report `docs/SM/CD01_CONDENSATE_DRESSING.md`; `tests/test_cd01_condensate_dressing.py` (2 tests); ledger 167 entries.

### DP-01 dressing strengths and the quark–lepton conflict (2026-10-05)

- `mtft.matter.dressing_profiles`: the massless light state and its weight, the dressing strengths z_u = 0.25 needs, the lepton requirement z_e = 3.09, the factor-12
  conflict.  CD-4 becomes a caveat.  Report `docs/SM/DP01_DRESSING_PROFILES.md`; `tests/test_dp01_dressing_profiles.py` (2 tests); ledger 171 entries.

### GT-01 the two gaugino towers (2026-10-05)

- `mtft.matter.gaugino_towers`: the W_R tower's density share with cusp-concentrated flux, the subspace strength bound, the reversal requirements.
  Report `docs/SM/GT01_GAUGINO_TOWERS.md`; `tests/test_gt01_gaugino_towers.py` (2 tests); ledger 175 entries.

### QL-01 the quark–lepton split (2026-10-06)

- `mtft.matter.quark_lepton_split`: the b–tau bound on dressing, the per-family lepton/quark targets, the vertex insertions, the Wilson-line limit.
  Report `docs/SM/QL01_QUARK_LEPTON_SPLIT.md`; `tests/test_ql01_quark_lepton_split.py` (2 tests); ledger 180 entries; CC-52 retracts CD-4.

### CI-01 the condensate insertion (2026-10-06)

- `mtft.matter.condensate_insertion`: the rank-one projector, the tail theorem, the running-independent entry pattern, the lever-arm bound.
  Report `docs/SM/CI01_CONDENSATE_INSERTION.md`; `tests/test_ci01_condensate_insertion.py` (2 tests); ledger 186 entries; CC-53 supersedes QL-3.

### BL-01 a B–L Wilson line on 143a1 (2026-10-06)

- `mtft.matter.bl_wilson_line`: twisted image sums over 2-torsion points, the lepton point against V_cb, the one-loop saddle structure, the image-pair factor.
  Report `docs/SM/BL01_BL_WILSON_LINE.md`; `tests/test_bl01_bl_wilson_line.py` (3 tests); ledger 192 entries; CC-54 caveats CI-5.

### FA-01 the flavour audit (2026-10-06)

- `mtft.matter.flavour_audit`: the scoreboard, image-summed V_cb and the corrected r, Gauss–Bonnet curvature radii, the arithmetic of X0(143)(Q).
  Report `docs/SM/FA01_FLAVOUR_AUDIT.md`; `tests/test_fa01_flavour_audit.py` (2 tests); ledger 198 entries; CC-55 supersedes CA-2.

### AG-01 the canonical Green's-function kernel (2026-10-07)

- `mtft.matter.arakelov_kernel`: the conformal kernel identity, the second-kind differential, frozen kernel values at the homes, the pre-registered r test.
  Provenance `scripts/arakelov_kernel_x0143.py` (cypari2). Report `docs/SM/AG01_ARAKELOV_KERNEL.md`; `tests/test_ag01_arakelov_kernel.py` (2 tests); ledger 204 entries; CC-56 caveats FA-2.

### FH-01 the top at the Fricke point (2026-10-07)

- `mtft.matter.fricke_home`: the Fricke point's image on 143a1, torus suppression of parity-odd tails, the four-way r concordance, the V_ub texture-zero relation.
  Report `docs/SM/FH01_FRICKE_HOME.md`; `tests/test_fh01_fricke_home.py` (2 tests); ledger 208 entries.

### LV-01 the lever at the jet point (2026-10-07)

- `mtft.matter.lever`: the two-arm texture formulas, the opposite-arm solution, the CP-angle test with real and twisted arms; CC-57 supersedes FH-4.
  Report `docs/SM/LV01_LEVER.md`; `tests/test_lv01_lever.py` (2 tests); ledger 213 entries.

### TL-01 the gamma tilt is a rotation (2026-10-07)

- `mtft.matter.tilt`: block-tilt cancellation, c as the down quarks' share of V_cb, the dressing bound on c, the pre-registered window test (fail at z = 0.25).
  Report `docs/SM/TL01_TILT.md`; `tests/test_tl01_tilt.py` (3 tests); ledger 219 entries; CC-58 caveats CA-1.

### DS-01 the down sector is parity-symmetric (2026-10-07)

- `mtft.matter.down_parity`: the condensate's reach by stack charges, the exact angle relation, the V_ub floor, the data solution, the tail cap.
  Report `docs/SM/DS01_DOWN_PARITY.md`; `tests/test_ds01_down_parity.py` (2 tests); ledger 225 entries; CC-59 supersedes TL-5.

### PA-01 the parity audit (2026-10-08)

- `mtft.matter.parity_audit`: the parity ledger of odd sources and their first-order reach, the odd/even split of the complex-structure moduli, the up-sector loophole, the verdict.
  Report `docs/SM/PA01_PARITY_AUDIT.md`; `tests/test_pa01_parity_audit.py` (2 tests); ledger 229 entries.

### OM-01 and YB-01 the two experiments (2026-10-08)

- `mtft.matter.odd_modes`: the Lefschetz count of forced odd modes under the torus flip, the zero-torus-degree design, the topological form of the parity verdict.
- `mtft.matter.bridged_top`: the bridged top Yukawa at the Fricke separation against the pre-registered y_t >= 0.5 (fail); closure of the two-point flavour design.
  Report `docs/SM/OM01_YB01_EXPERIMENTS.md`; `tests/test_om01_yb01_experiments.py` (3 tests); ledger 235 entries; CC-60 caveats PA-4.

### Programme map and consolidation (2026-10-08)

- `mtft.matter.puzzle`: the programme map. Ten pieces (I–X, numbered in dependency order), each a question the programme asks, with its status
  (stands, partial, open, closed), verdict, owned ledger waves and its own sub-puzzles (open, answered, failed, moot; what would settle each open one);
  the pre-registered tests and their outcomes; the honest negatives; the chronology of the waves; and the disposition of every report in docs/SM.
  Every ledger claim belongs to exactly one piece. `python -m mtft matter puzzle [--piece ID | --open | --json | --markdown | --dot]`.
  Generated `docs/SM/MATTER_PUZZLE_MAP.md`; `tests/test_puzzle_map.py` (9 tests) checks it against the ledger, the module docstrings and the reports.
- CC-61 retracts LEP-4: the equal-overlap texture has several exact solutions per sector, and leptons also fit with the charge-0 jet dominant (LEP-6).
  Exposed by the environment change below; `leptons.fit_alignment` now returns a stated canonical representative (the rho branch nearest 1/2, then
  the smallest coefficients) and `leptons.all_alignments` the full set. CP-01's zero-flip test and LEP-01's fit test no longer depend on the optimizer path.
- CC-62 (bookkeeping): the parity audit's wave is relabelled PAR-01 (PA-01 is Path A's); claims PA-1 to PA-4 unchanged.
- Environment: the suite now runs on Python 3.13, NumPy 2.5, SciPy 1.18. Two tests depended on which exact solution SciPy's least_squares reached
  first; both are now path-independent (above).
- Ledger 238 entries (187 active, 33 caveat, 12 superseded, 6 retracted).
  `docs/SM/MATTER_PROGRAM_MAP.md` (2026-10-02) is marked superseded by the generated map. Notes for the integration pass:
  `docs/SM/MATTER_PROGRAMME_INTEGRATION_NOTES.md`. Full routine suite: 1279 passed, 4 skipped (PARI/GP), 40 slow deselected.

### SQ-01 star quotients of supersingular semiprime levels: is X0(143) special? (2026-10-08)

- `mtft.matter.star_quotients`: the genus of X0(N)* = X0(N)/<W_p, W_q> at all 105 levels N = pq with p < q supersingular, by two routes —
  newform signs (PARI, frozen; 99 levels) and class-number fixed points with Eichler symbols and the Lefschetz formula (pure Python, live; all
  105) — which agree wherever both run, and at 75 levels with three or four primes up to 700. Sphere at 21 levels, elliptic curve at 13
  (57, 58, 65, 77, 82, 91, 118, 123, 141, 142, 143, 145, 155; all rank 1), genus >= 2 at 71. What singles out 143 is arithmetic only (twin
  primes; new dimension p with consecutive primes; orbit pattern (1, 4, 6)); 123 = 3 x 41 matches its index, genus, cusps and quotient type.
  Every sphere level is a Monster element order (GAP CTblLib). `python -m mtft matter levels [--level N] [--json]`.
- Ledger SQ-1 to SQ-8; programme map: wave SQ-01 on piece II (sub-puzzles II.7–II.9, negative N20), the open sub-puzzle IV.6 (why this level,
  with its test pre-registered as SQ-8), and SQ-7 on IX.9 (another level would not reopen quark flavour).
- Provenance: `scripts/star_quotient_scan.py` (needs cypari2) reproduces the frozen tables; `--check` compares them.
- `tests/test_sq01_star_quotients.py` (12 tests; one needs cypari2). Ledger 246 entries (195 active, 33 caveat, 12 superseded, 6 retracted).
  Full routine suite: 1295 passed, 40 slow deselected (PARI/GP 2.15.4 came with the GAP install, so the 4 PARI/GP tests no longer skip).
