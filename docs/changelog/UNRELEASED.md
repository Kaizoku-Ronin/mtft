# Unreleased — after v1.0.0

Everything after v0.34.0 shipped in v1.0.0 ([release notes](CHANGELOG_v100.md)). Since then:

## CC-71 (2026-10-29): the Tano metric of the diatomic gases — S°(X₂) over the mass of X₂

- `mtft.tano_metric` divided the standard entropy of H₂, N₂ and O₂, which is tabulated per mole of molecules, by the molar mass of one
  atom. That doubled their Tano values: H 129.64 → 64.82, N 13.68 → 6.84, O 12.82 → 6.41 J g⁻¹ K⁻¹. `Element` gains `atoms` (atoms per
  standard-state formula unit; default 1, and 2 for H, N, O), and `tano_metric` now divides by `atoms × molar_mass`. The values before the
  correction stay exposed as `tano_metric.TANO_METRIC_RETRACTED_CC71` (append-only). `seebeck_diagnostic` computes its own T_a and now
  uses the same normalisation; nothing changes for the metals it is meant for.
- Second route: the three gases counted from first principles reproduce the stored S° per mole of molecules: H₂ 130.682 against 130.68
  (ortho and para hydrogen; counting every J once instead gives 5.71 more, close to R ln 2), N₂ 191.564 against 191.61, O₂ 205.128 against
  205.14 (only odd rotational levels exist, because ¹⁶O has spin 0; R ln 3 for the electrons' triplet). So the division by 2M is the right one.
- Downstream, the default screen (`materials_screening()`, `python -m mtft screen`) moves with it: LaH₁₀ λ_eff 2.092 → 1.444 and T_c
  930.1 → 750.5 K, H₃S 743.1 → 598.7 K, and the O/Ba proxy for YBCO 125.6 → 109.8 K; MgB₂, Al and Pb are unchanged. The screen's other weakness, the form
  Θ e^(−1/λ) without the (1 + λ) renormalisation or μ*, is a model choice rather than part of this correction and stays as it is (DIAGNOSTIC).
- Found on 2026-10-28 in a review of Roger's Tano repository workbook, where the same doubling affects seven rows outside the package
  (H, N, O, F, Cl, Br, I). `tests/test_cc71_tano_diatomic.py` (11 tests).

## CC-72 (2026-11-15): the Burning Ship fold exponents are a caveat — Paper 4's off-axis table is not one cascade

- `mtft.constants.BETA_X = +0.29` and `BETA_Y = −0.21` (and `DELTA_X`, `DELTA_Y`, defined from them) were fitted to the five points of
  Paper 4's off-axis table, periods 48 to 768. The four gaps between those points are 2.42, 1.81, 2.08 and 1.34 × 10⁻⁴: they shrink 1.81
  times over three doublings, where a period-doubling cascade shrinks δ³ = 101.8 times, and the third gap is larger than the second. The
  five points are not consecutive doublings of one cascade, so the exponents are unconfirmed.
- No number changes and nothing is removed (append-only). `mtft.burning_ship.cc72_cascade_check()` states the comparison and lists the
  dependents (`AnisotropicScaling`, `satellite_mass_scale`, `MassRatios.m_tau_over_m_mu`); `mtft.burning_ship.real_axis_cascade()` is
  the control: on the real axis the Burning Ship's map is x² + c, and its superstable parameters do shrink by Feigenbaum's δ
  (δ₁₀ = 4.6691951, c∞ = −1.401155189092). `DELTA_X_MEASURED` and `DELTA_Y_MEASURED` (Paper 16) were not examined.
- What would settle it: an off-axis cascade search that confirms each doubling before any fit.
- Found on 2026-11-06 in the Modular Plane lab (record MP05, Room 11). The gaps are as transcribed there from Paper 4's table; the
  transcription is to be checked against the paper at audit. `tests/test_cc72_burning_ship_cascade.py` (4 tests).

## MA-01 (2026-11-15): `mtft.maass` — the Maass spectrum of Γ₀(143) by Hejhal's method, the cusps' scattering matrix, the Ford domain

A new subpackage, ported from the Modular Plane lab's records MP14 to MP20. Import it on demand (`import mtft.maass`); nothing is added
to the top-level namespace. Report: `docs/SM/MA01_MAASS_SPECTRUM.md`. It makes no statement about physics.

- `mtft.maass.kbessel`: K_{iR}(x) by a trapezoid rule (R ≤ 3), 2.1 × 10⁻¹⁴ against mpmath on the gate's sample; mpmath above R = 3.
- `mtft.maass.domain`: the reduction of a point into the Ford domain of Γ₀*(N) (Γ₀(N) with its Atkin–Lehner matrices, N squarefree),
  with the Atkin–Lehner characters and, for odd N, the sign of the η-twist χ = ν_η¹². Gated on Eisenstein series (no Maass form
  involved) and on words in S and T.
- `mtft.maass.hejhal`: the solver. `refine` finds the R at which a coefficient agrees on two horocycles; `sigma_min` is the pole-free
  detector. Gates: LMFDB's form 13.1 at level 13 (0.97081541), three forms at level 11, three at level 143, SL(2, ℤ)'s first two
  eigenvalues with mpmath's Bessel function (marked slow), and controls that a wrong sign, parity or R does not pass.
- `mtft.maass.spectrum`: the frozen table. 41 cusp forms below R = 2.0 (even 4, 4, 4, 7; odd 6, 6, 5, 5 by Atkin–Lehner sector) and 42
  forms of the η-twist below R = 1.8, each with 400 coefficients and SHA-1 provenance. The first eigenvalue is λ₁ = ¼ + 0.3738518157²
  = 0.3897652 (odd, signs (−,−)). Every eigenvalue was found again by the package's solver when the table was built, within 5 × 10⁻¹² of
  the lab's. `evaluate` gives a form's value anywhere in the upper half-plane (a grid of it is a picture); `rederive` repeats the solver for one
  form; `tones` and `sonify` hand the table to `mtft.music` as pitches, root·√(λ/λ₁).
- `mtft.maass.scattering`: Φ(s) = φ(s) M₁₁(s) ⊗ M₁₃(s) and the Eisenstein series of the four cusps (mpmath). Φ(s)Φ(1 − s) = I,
  unitarity on Re s = ½, residue 1/(56π) at s = 1, Φ(½) = −I: every Eisenstein series vanishes at the bottom of the continuum.
- `mtft.maass.ford`: the Ford domain of Γ₀*(N) (17 tops, 24 vertices and 16 distinct heights at N = 143; lowest vertex at height
  0.009851) with a panel quadrature: the area 14π to 2 × 10⁻¹¹, norms, overlaps, triple products, participation ratios (the lab's seven
  values to 10⁻³). Products whose character is not trivial are zero on the curve by symmetry and are returned as zero without integrating.
- Ledger (`mtft.matter.ledger`), wave MA-01 in piece IV: MA-1 upgrades the diagnostic figure of VE-5 (whose text stands) to λ₁ =
  0.3897652, with SPEC-01's finite elements as the independent route (0.3925, and the next eight eigenvalues within 1.5 %, the double one
  included); MA-2, the tower; MA-3, the edge of the continuum; MA-4, the trace formula's certificate.
- `mtft.maass.trace`: the Selberg trace formula by sector (odd; even with the continuum; the η-twist's total) on class-number weights
  frozen for lengths up to 17 (`_data/trace_weights.npz`, 276 kB, recomputed with PARI and equal to the lab's tables). By positivity no
  cusp form is missing below R = 1.959 in any untwisted sector and parity, nor below R = 1.607 in the twisted lists together, and no
  residual comes near the 2 an eigenvalue below ¼ would add (ledger MA-4). Controls: PSL(2, ℤ)'s whole formula cancels to 10⁻¹⁰; level
  11's three forms to 2 × 10⁻⁹; a removed form is located to 0.002; a spurious one breaks positivity. `l1_dirichlet` and `weight` are a
  second route for the class numbers with no PARI (150 weights to 10⁻⁹ in the tests).
- **What stays recorded and not gated.** The lab's sharper radii (1.98 untwisted, 1.8 twisted) used a narrower test function and lengths
  to 23.5, megabytes of weights that are not shipped. The class numbers are PARI's (conditional on GRH in general); the sample and the
  level-1 cancellation are the checks on them.
- Provenance: `studies/maass_2026nov/` (the two build scripts and the lab's scripts and number files). Tests: `test_maass_solver.py`
  (10), `test_maass_spectrum.py` (10), `test_maass_scattering.py` (5), `test_maass_ford.py` (5), `test_maass_trace.py` (8),
  `test_ma01_maass_rung.py` (3).


## LN-01 (2026-11-16): `mtft.lines` — the closed lines of Γ₀(143)\PSL(2, ℝ): thirty integers, every crossing, linking numbers

A new subpackage, ported from the Modular Plane lab's records MP12, MP13 and MP14. Import it on demand (`import mtft.lines`); nothing is
added to the top-level namespace. Report: `docs/SM/LN01_CLOSED_LINES.md`. It is mathematics of the curve's unit tangent bundle and makes
no statement about gravity.

- `mtft.lines.darts`: P¹(ℤ/N) and the graph on its darts (at 143: 56 vertices, 84 edges, 4 faces of lengths 1, 11, 13, 143; b₁ = 29).
  `Graph.cls` gives the class of an element of Γ₀(N) in ℤ^b₁ from its word in S and T; `Graph.generators` are the 29 generators of
  Γ₀(143). Levels with elliptic points are refused; `Darts` (cycle types of any matrix) works at every level.
- `mtft.lines.rademacher`: Dedekind sums by reciprocity, Rademacher's Φ and Ψ, the windings of η-quotients round the cusps (`windings`
  for the three units of level 143) and Ligozat's orders at the cusps, which the windings on the cusp loops must equal.
- `mtft.lines.closed`: `lifts(A, N)` builds the closed lines over a hyperbolic matrix, each with ℓ, its darts, γ = g A^ℓ g⁻¹, Ψ and its
  class; `lines_143()` is the lab's 28 over the cat map, Sphereland 10's motion and the Pell unit, equal to its table integer for
  integer. `transfer` gates the sums over a motion's lifts: classes to zero, fibre numbers to 168Ψ(A).
- `mtft.lines.periods`: the 26 periods of a class on the Hecke blocks (2 + 4 + 8 + 12) by two routes, the lab's frozen coordinates and
  the package's own Manin symbols in `mtft.hecke`, compared through zeros, ranks, relations and the lattice they span; `thirty(γ)`
  returns all thirty integers of an element. `lattice_143a1` gives Ω₁ = −0.984962 + 2.015772i, Ω₂ = −1.969923 from the q-series
  (CERTIFIED 1e-8), and `stored_periods_in_lattice` checks them against `mtft.periods`' PARI periods (index 1).
- `mtft.lines.involutions`: W₁₁, W₁₃, W₁₄₃ and the mirror as integer matrices on ℤ²⁹, their signs on the blocks, and
  `psi_shift`: an Atkin–Lehner involution does not keep the fibre number, Ψ(W_Q γ W_Q⁻¹) = Ψ(γ) + the winding of (η(Qτ)/η(τ))²⁴
  along γ. `conjecture_p13` computes the lab's missed conjecture (below).
- `mtft.lines.cutting`: the cutting sequence of a closed geodesic along the fundamental domain in rational arithmetic (4 arcs a period
  for the cat map, 3 for Sphereland 10's motion, 25 for the Pell unit); proper powers are refused.
- `mtft.lines.crossings`: `Crossings(N, motions)` counts every crossing of the closed geodesics over some motions, with signs, in two
  half-open conventions that must agree; `linking()` is Lk(K_γ, K_δ) = −I(γ, δ) for lines taken with their reverses (the lab's
  derivation). At 143 with the three motions: 24 geodesics, 3214 crossings between them, 2918 self-crossings, the cat unions linking
  0, −6, −8, −24; with three more motions, MP14's 38 geodesics. The mirror's action on the geodesics and the vanishing of the Pell rows
  of the signed matrix are gated.
- `mtft.lines.lorenz`: Ghys's Lorenz template. `lorenz_word`, `lorenz_linking`, and `transfer_sum`: the total of crossings between the
  geodesics over two motions at level N is ψ(N) times the template's sum. Gated on all fifteen pairs of six motions at levels 143, 11
  and 23, and on a pair no record holds.
- Ledger (`mtft.matter.ledger`), wave LN-01 in piece II: LN-1 to LN-8. **LN-6 records a miss**: the lab's pre-registered conjecture that
  every line with no 143a1 component is negated by an Atkin–Lehner involution is false (lines 22 and 24); it is also negative N23 of
  the puzzle map. LN-8 (the linking formula) is STRUCTURAL with a caveat: its sign assumes an orientation that was not checked.
- **Found during the build.** Five tests failed on their first runs, all through something their author had written without computing
  it; none through the port. Four were corrected and say so. The fifth was a wrong belief, that the involutions keep Ψ; the law above
  replaced it (LN-5).
- **What is not here.** A mechanism for the two zeros of LN-6; levels with elliptic points (so the lab's shadows on X₀(13)); any
  physical reading of the thirty integers.
- Provenance: `studies/lines_2026nov/` (the build script and the lab's scripts and number files for MP12 and MP13; MP14's are already
  under `studies/maass_2026nov/lab/mp14`). Data: `src/mtft/lines/_data/lines_143.json` (28 kB) with SHA-1 in `PROVENANCE.txt`.
  Tests: `test_lines_darts.py` (9), `test_lines_closed.py` (6), `test_lines_periods.py` (7), `test_lines_crossings.py` (10),
  `test_ln01_closed_lines.py` (4).

## `mtft.plane` (2026-11-16): the Modular Plane lab's kernel — one matrix read eight ways, and the lattice of quantities

A new subpackage with two modules, ported from the core of the Modular Plane lab (version 8.3). Import it on demand
(`import mtft.plane`). It has no ledger entries: nothing in it is a claim about X₀(143) or about physics.

- `mtft.plane.kernel`: `read(M)` returns the lab's eight readings of a real 2 × 2 matrix as a motion of the plane (a motion of τ; a
  form and a momentum; a number a + εb in ℝ[M], complex, dual or split; a flow; a field equation; one mode's vacuum; two levels; its
  numerical range), with its invariants and kind; `read(M, generator=True)` reads it as a generator, whose motion is `expm(L, t)`.
  Entries may be integers, Fractions or floats; exact entries give exact rational quantities. For an integer matrix of determinant 1
  it adds Rademacher's Ψ and the Lorenz word, and with `level=N` the cycle type on P¹(ℤ/N), all through `mtft.lines`.
  `numerical_range`, `jacobi_monodromy`, `fixed_points`, `decompose`, `region` are also callable alone.
- Gates from the lab's pre-registered table (record MP06): the Pell unit [[12, 1], [143, 12]] has rapidity 3.176313, length 6.352626,
  a numerical range centred at 12 with semi-axes 72 and 71 and numerical radius 84, ‖M‖/w = 1.7144, and Jacobi monodromy of trace
  574 = tr M²; the dual unit's range is the disk of radius ½, where Crouzeix's bound 2 is sharp; the cat map's is the segment [1/φ², φ²].
- **Three identities, stated and gated** (elementary; noticed while porting). Write M = αI + s₁S₁ + s₂S₂ + wJ with strain
  s = √(s₁² + s₂²). Then W(M) is the ellipse centred at α with semi-axis s along the real axis and |w| along the imaginary one, so the
  lab's seventh reading (two levels: gap 2s, rotation w) and its eighth (the range) carry the same two numbers; the spectral norm is
  |α + iw| + s; and for det M = 1, s² = sinh²(d/2), the quanta of the sixth reading (for the Pell unit, 72² = 5184). Checked on 500
  random matrices against the lab core's own computation (kept as `lab_numerical_range`) to 10⁻¹² and against numpy with no formula for
  the range (the top eigenvalue of the Hermitian part of e^{−iφ}M) to 10⁻⁹. The numerical radius has a closed form here; the lab
  searched for it.
- `mtft.plane.quantities`: the lattice ℤ⁴ of dimensions M^a L^b T^c I^d with the lab's catalogue (33 + 26 + 13 entries, identical to
  the lab's). `lattice()` computes, through `mtft.integral_lattice`, that ⟨c, ħ, ε₀⟩ has Smith invariants 1, 1, 2 (so natural units
  leave ℤ ⊕ ℤ/2: a mass dimension and a charge parity, `natural`) and that the Planck sublattice ⟨c, ħ, ε₀, G⟩ has invariants
  1, 1, 2, 2 and index 4. `planck_decompose` writes a site in those terms with integer coefficients; `site(v)` is the lab's prospector:
  what lives at a site, which pairs of catalogued quantities build it, and its class.
- What stays out: the lab's page, its matrix glossary and its instruments (the tidal catalogue, Newton's method, Morley's triangle).
  They are teaching material.
- **The port against the original.** The lab's own JavaScript core was run on 52 matrices and 49 sites and its output frozen in
  `studies/plane_2026nov/` with the core's SHA-256. Two tests compare the package with it: every number agrees to 3 · 10⁻¹⁶ relative
  and every label, integer and name exactly. The core itself is the lab's and is not shipped.
- Tests: `test_plane_kernel.py` (7), `test_plane_quantities.py` (6).

## GB-01 (2026-11-16): `mtft.maass.glueballs` — H24 against the lattice, with every registered forecast and two that are sealed

A new module in `mtft.maass`, ported from the lab's records MP15, MP18, MP19 and MP20. Report: `docs/SM/GB01_GLUEBALL_COMPARISON.md`.
Class PHENO: a comparison under a stated dictionary. **Its result is negative: H24 as written fails.**

- The lattice table (`LATTICE`): six values of M/√σ for SU(3) in the continuum limit, from Athenodorou and Teper (JHEP 11 (2020) 172),
  as the lab transcribed them on 12 November. They were not re-read against the paper in this release's build, on purpose: the same
  tables would show the states the sealed forecasts are about (`TRANSCRIPTION`).
- `reading_b()`: with the cusp forms only and C the Fricke sign, all five predictions miss (the lowest cusp form is odd; ratios 0.870,
  2.236, 1.247 against windows near 1.55, 1.72, 2.14; a C = −1 form at 2.152).
- `comparison()` and `needed_shifts()`: reading a with the η-twist as C = −1, state by state, in the form of `mtft.falsify`: the 0⁻⁺ at
  1.559 against 1.549(16) passes (the one registered pass, recorded as a coincidence without a mechanism); the 0⁺⁺* is 4.5σ off; the
  C = −1 states are 9.7σ and 7.5σ too light, a registered miss. Reading a fails on structure in any case: its 0⁺⁺ is the edge of a
  continuum.
- `self_interaction_fit()`, `one_coupling_fit()`, `minus_sector_fit()`, `twist_energy()`: MP20's fits recomputed from the frozen
  table: χ² = 10.2 on 2 for one coupling on C = +1; 170 on 4 for one coupling on all five; μ = 0.218(18) with χ² = 1.57 on 1 for a
  constant twist energy. The lab had kept no code for the 170: rebuilding it showed it used independent errors (178 with the shared
  error of the 0⁺⁺); both are in the function. `participation(recompute=True)` repeats the seven participation ratios with
  `mtft.maass.ford` (10⁻⁸ from the lab's).
- **Two open forecasts, sealed** (`SEALS`, GB-S1): a second odd C = −1 state at 2.87(7) and a second even one at 3.24(7), in units of
  the 0⁺⁺. Registered by the lab on 14 November before any lookup; the criterion was fixed and the text, numbers and rule hash-sealed
  on 16 November, still before any lookup. `evaluate_open` runs the sealed rule and refuses a seal that no longer verifies. The lookup
  has not been made.
- `register_forecast` and `compare_forecast`: a new dictionary is compared with the lattice only from a sealed record, forecast first.
- `register()`: the 23 registered forecasts and checks of the four records, each with the probability given beforehand where there
  was one, and its outcome.
- Ledger, wave GB-01 in piece IV: GB-1 to GB-5; negative N24 on the puzzle map. GB-5 (the twist energy) is a caveat and is not part of
  the negative, because its forecasts are open.
- Provenance: `studies/maass_2026nov/lab/` (MP15, MP18, MP19, MP20, and `needed_shifts.json`, added now). Tests:
  `test_maass_glueballs.py` (9), `test_gb01_glueball_comparison.py` (4).

## `mtft.gravity` (2026-11-16): knotted weak fields — the tube, the squared hopfion's lines of force, a frozen tetrad

A new subpackage with four modules, ported from the Modular Plane lab's records MP09 and MP21 to MP23. Import it on demand
(`import mtft.gravity`). Report: `docs/gravity/KNOTTED_FIELDS.md`. Status: model. Every field in it is a weak wave on flat space;
Newton's constant is not derived anywhere in MTFT, and the module's first page says so. It has no ledger entries yet: the ledger of
the gravity map is the next wave's work.

- **A search of the literature came first**, as the release plan required, and nothing here is called new. The (p, q) fields and their
  invariant are Kedia, Bialynicki-Birula, Peralta-Salas and Irvine's (2013); gravity's lines of force as the light's lines turned by a
  local duality angle are Swearngin et al.'s (2013); the integrable potential the later lines land on is Garnier's; and a helicity of the
  gravitational hopfion was computed in transverse-traceless gauge, with a different and time-dependent answer, by Huang, Long, Zhao
  and Zhou (arXiv:2608.30154) three months before the lab's record. What the search did not find is listed as "not found", which is
  not "new": the rotation number, the shear, the reduction of the later lines to the integrable system, the frozen tetrad.
- `fields`: the labels, F_pq, the squared hopfion, the null flow, and `integrate`, the package's own routine for lines of force (an
  extrapolation of the modified midpoint rule; no scipy), with crossings located to the routine's own accuracy.
- `knotted`: ν = √(2pq(p + q)), the circuit matrix, the pairs where it is the identity ((1, 1), (1, 8), (1, 49), (32, 49) among coprime
  p ≤ q ≤ 60; 22 pairs without "coprime"); the profile ι(h) on every torus; the shear at the core line,
  ν[3pq(p + q) − 2(p² + pq + q²)]/(12pq(p + q)), by three routes (closed form, the integral at 30 digits, the section integrated as a
  differential equation); the derivation exactly, for symbolic p and q. Out of sample: helicity 2π² pq·p! q!/(p + q)! and energy (p + q)
  times it, by quadrature on the label sphere and in space; their ratio is Kedia et al.'s.
- `hopfion`: the reduced system H = |P|² + |s² − σ²|²/t², the second invariant I₂, the winding as a ratio of two periods, MP22's table of
  W − 1 (to 10⁻⁹ against the lab's file), the two closed lines. In space, on the field itself: I₂ is kept to 10⁻¹³ of its range over 100
  circuits; the closed lines are linked −1 and +1; the winding measured along an orbit equals the integrals to 10⁻¹⁵.
- `tetrad`: B_μ = Re[(i a² (∂_μ b) dγ − b² (∂_μ a) dβ)/d²]. `exhibit` shows it frozen into the null flow with curvature equal to the squared
  hopfion and zero Ricci tensor: to 2·10⁻¹⁵ in double precision, to 10⁻⁴⁰ at 40 digits, and exactly over the Gaussian rationals.
  `certificate` makes the exact zeros a statement about every event: the differences are polynomials of degree at most 22, so ten
  random rational events leave a chance below 3·10⁻²¹ of a false pass (Schwartz–Zippel). The helicities (2/3)π² and π², on the label
  sphere and on two slices of space.
- Kept: every recorded miss of the four records (report, section 5). Open: T-G3 not attempted; H-G7 not tested.
- Found in the records during the build: MP09's "positive factor" is negative; MP23's file has the sphere's helicities with a minus
  sign, from the order of its coordinates.
- Provenance: `studies/gravity_2026nov/`. Tests: `test_gravity_knotted.py` (10), `test_gravity_hopfion.py` (7), `test_gravity_tetrad.py` (6).

## GB-02 (2026-11-16): E5 — the two sealed forecasts of the twist energy, looked up

The experiment GB-01 left sealed. Report: `docs/SM/GB02_E5_SEALED_FORECASTS.md`. Class PHENO. **H24's verdict does not change.**

- How a table is read was written down before the paper was opened (which table, what "first excited" means, single-spacing values are
  no data, two readings that must agree). The seal GB-S1 verifies as sealed; its outcome fields, which are outside its hash, are filled.
- **The odd forecast passed, at 1.4σ.** The continuum table of Athenodorou and Teper lists one excited C = −1 state of A₁, 10.39(37), the
  second A₁⁻⁻ state: 3.051(110) in units of the 0⁺⁺ against the forecast 2.866(72). **The even forecast has no lattice value**: A₁⁺⁻ has no
  excited state in the continuum table. By the seal's own words a pass is weak support.
- **Why it is worth little** (`e5()`): at the lattice's error any first excited state from 8.87 to 10.65 passes, and the ground state is at
  9.78(46); the pair is nearly degenerate and noisy (9.93(49), 10.03(47) in the authors' later paper); μ is one fitted number.
- **The lab's six transcribed values are right**, digit for digit (`transcription_check()`). What its record left out: the paper assigns the
  A₁⁺⁻ ground state to the 4⁺⁻ ground state, with "significant uncertainty". Under the dictionary's D2 that state is not a J = 0 state: the
  even twisted form of GB-01's comparison then has no lattice partner, and the twist energy was fitted to two ground states of which one
  probably should not have been used. Neither helps H24.
- **A limit of the reading, stated.** The paper was read as flattened text by automated extraction, twice; the PDF's layout was not seen.
  That 10.39(37) stands under P=−,C=− is identified (`identify_column()`): by redoing the paper's continuum extrapolation from its own
  finite-spacing tables (it reproduces 3.405(21), 9.32(28) and 10.39(37), and nothing under P=+,C=− gives the last without a poor fit),
  and by the same authors' later SU(3) table, whose cells are kept. The sealed rule's overall verdict is the same for either column.
- Kept: `LATTICE` is not edited; GB-5 and the GB-01 report stand as written, each with a dated pointer. The earlier tests that pinned
  "pending" and "open" now pin the outcome.
- Ledger, wave GB-02 in piece IV: GB-6 (the outcome; caveat) and GB-7 (the spin of the A₁⁺⁻ ground state). N24 is unchanged.
- **Read since.** Roger Tano read the PDF's page (a screenshot of Table 14), and the file's text with its layout kept shows the same:
  10.39(37) stands under P=−,C=−, as identified (`EXCITED_READ`, `TABLE14_A1_COLUMNS`). Every number held from the paper for this
  experiment was compared with the layout and agrees.
- Tests: `test_gb02_e5_sealed_forecasts.py` (6).
