# mtft v0.11.4 — dynamical units wave (2026-08-05)

The two-clock dimensionless ledger and the cycle-space map. Zero new
conjectures; one honest negative upgraded to a theorem. Wave code:
Claude Opus 5, run and certified in-session; independently audited
before push (Addendum BT, Kimi K3 — every headline number re-derived
on the auditor's own machinery: independent P^1(Z/143) chain,
independent exact integer linear algebra, independent path-based
Hecke, independent a_p tables to p <= 2500).

## du01 — the two-clock ledger (`studies/du01_two_clock_ledger.py`)

Before any y-to-units map can be proposed, the box's unit-free content
is frozen: the complete table of dimensionless spectral invariants any
admissible units map must preserve.

* Graph clock (CURVE), zero parameters. Exactly one self-loop on the
  56-vertex Farey dual => tr L = 166 = 168 - 2 (EXACT). Integer
  eigenvalues {0,1,2,4,5}, each simple, proved by fraction-free
  Bareiss rank over Z (3 and 6 absent). Kirchhoff spanning-tree count
  tau = 3518081582959364640 = 2^5 * 3 * 5 * 17 * 941 * 101921 * 4495339
  (EXACT integer determinant; E2 float route at 2e-14; no 7/11/13
  factors — observed, no interpretation filed). Spine identity
  b1 = E - V + 1 = 29 = 2g + (cusps - 1) = 26 + 3 (EXACT).
* Hecke clock (ITERATE), zero parameters. The 13 T2 lines, each
  doubled (max pair split 1.3e-14), trace -2; orbit traces
  cross-certified against the mtft.x0_143 oracle to n <= 50 at 9.2e-14
  (E2: Manin-symbol route vs LMFDB-validated corpus route). Ramanujan
  pass-band certificate: max |a_p|/(2 sqrt p) = 0.959220 at good
  p <= 47; bad p in {11,13} carry |a_p| = 1 (Atkin-Lehner, EXACT).
  Uniform-weight a2 variances EXACT from the integer Hecke polynomials
  via Newton identities: Var(f1) = 0, Var(f2) = 35/16, Var(f3) = 10/3
  (E2 vs float engine at 5.3e-15).
* The ledger. Graph-clock ratios to the Fiedler line (span
  1 .. 20.950319); Hecke-side invariants frozen as the joint a_p table
  plus per-orbit spreads. Anchor count: 2 — one external anchor per
  uncoupled clock; the H4 lifetime non-correspondence (v0.11.0) is
  carried, not retried.
* Physical ladder (PHENO). Pinned rungs f(1S-2S) = 2466061413187035 Hz
  (Parthey 2011), f(21 cm) = 1420405751.7667 Hz (Hellwig 1970),
  c R_inf = 3.2898419602500e15 Hz (CODATA 2022); unpinned rungs left
  BLANK per repository policy. Pure-number targets:
  f(1S-2S)/f(21cm) = 1736166.8735287427, f(21cm)/cRinf =
  4.3175501101e-07, f(1S-2S)/cRinf = 0.749598747594.

## du02 — the cycle-space map (`studies/du02_cycle_space_map.py`)

The spine identity b1 = 29 = 26 + 3 located the only geometric
coupling channel between the two clocks; this study computes it.

* LEFSCHETZ DUALITY from raw combinatorics. Each dual edge crosses the
  oriented geodesic of its Manin symbol once (+1 on x, -1 on x.sigma);
  the 29 x 29 pairing between an integer fundamental-cycle basis and
  the Manin quotient basis has det D = -1 (EXACT, Bareiss) — perfect,
  unimodular Poincare-Lefschetz duality for X0(143). All 140 sigma/tau
  relations are killed by all 29 basis cycles, worst residue 0.
* Cusp links realized. The unique cycle pairing as the boundary
  functional of each cusp is integral (forced by det D = +-1); the
  four links sum to zero; link lattice rank 3; the width census over
  168 cosets is {1:143, 11:13, 13:11, 143:1}; and du01's self-loop IS
  the link of the divisor-143 (width-1) cusp — one triangle wraps the
  width-1 cusp (EXACT identification).
* Transported Hecke. T2* = D^{-T} T2^T D^T is INTEGRAL on the cycle
  lattice (Hecke preserves integer cycles), preserves the link
  lattice, the Eisenstein 3-block has charpoly exactly (x - 3)^3, and
  charpoly(T2*) = charpoly(cuspidal T2) * (x - 3)^3 — the transported
  clock carries the SAME 13 lines (E2 route pair: duality transport vs
  the engine's independent restrict_to_cuspidal).
* Hodge decomposition of edge space. Delta_1 = d1^T d1 + d2 d2^T has
  dim ker = 26 (the harmonic homology); its nonzero spectrum is
  exactly {55 graph-clock levels} u {3 link-Gram levels} at 1.1e-14;
  trace(B1 B1^T) = 166 re-derived by a second route.
* THEOREM (free-level obstruction). On the shared 26-dim stage the
  free graph clock acts as ZERO (harmonic kernel, EXACT) while the
  Hecke clock runs all 13 lines doubled — so no free-evolution
  exchange rate chi_H/chi_g can be formed on homology. The H4 lifetime
  non-correspondence is thereby EXPLAINED, not merely recorded: it was
  structurally forced. An internal exchange rate requires an
  interaction lifting the harmonic degeneracy; the splitting pattern
  of the 26 harmonic modes against the 13 Hecke lines is the box's
  dispersion relation (the du03 program).

## Drawn-loop stage (React)

`viz/MTFT_DrawnLoop.jsx`, with the interface data exported by du02
stage F (`du02_interface_data.json`, regenerable; embedded verbatim in
the component). Draw
closed walks on the 56-triangle dual graph; on closure the stage
returns the 29 exact intersection numbers with the Manin arcs, the
cuspidal/Eisenstein decomposition (signature basis, DIAGNOSTIC
Euclidean normalization, cond 132.6), orbit weights (DIAGNOSTIC), and
Zeno survival curves with the EXACT rates 35/16 and 10/3. The four
cusp links ship as one-click pure-Eisenstein demos. Deferred honestly:
loop-x-loop Goldman intersection (sign-convention audit, session 3).

## Media

`viz/wavepacket_X0143.mp4` (1120x1008, 12 fps) and its auto-looping GIF
derivative: unitary wavepacket evolution on the dual graph — the du01
graph clock in motion. Embedded in the README under the hero block.

## Gates

`tests/test_dynamical_units.py` runs both study mains (about 4 s) with
ledger output redirected to tmp and re-asserts the headline EXACT
certificates from the outside: tr L = 166, integer spectrum, Kirchhoff
tau and its factorization, b1 = 29, variances 0 / 35/16 / 10/3, anchor
count 2, det D = -1, width census, self-loop = width-1 link, link
supports, T2* integrality, (x - 3)^3, charpoly factorization,
signature counts, link-Gram levels.
