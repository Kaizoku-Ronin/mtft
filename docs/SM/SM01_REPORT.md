# SM-01 report — steps 2–4 on X0(143), step 5 specified (2026-09-13)

Preregistration: SM01_PREREGISTRATION.md.  Tool: `mtft.surface.smflux` (v0.30.0).  Data: sm01_5stack_solutions.json.

## Step 2 — representation table (EXACT, integer arithmetic)
Model class: U(3)_c x U(2)_L x U(1)^k with one line-bundle flux per stack; bifundamental (x, y-bar) has net
left-handed multiplicity m_x − m_y with internal bundle S0 (x) L_x L_y^{-1}; hypercharge y_x in Z/6 (y_c = 1/6).
- Four stacks: NO model gives three chiral families (exhaustive |m| <= 8; the multiplicities telescope around
  the hypercharge-consistent cycle of representations, so the negative is a theorem of the class, for any curve).
- Five stacks: 30 solutions (|m| <= 6, y in [−1,1]); all have m_L = −3 and at least six vector-like doublet
  pairs (1,2,±1/2) — the Higgs candidates; right-handed neutrinos appear with net count in {±1, ±3}.
- Representative M1: degrees (c, L, a, b, d) = (0, −3, 3, 3, 0), y = (1/6, 0, −1/2, 1/2, −1/2):
    (c,L-bar) Q x3 | (c-bar,a) u^c x3 | (c-bar,b) d^c x3 | (L-bar,a) L x6 | (L-bar,b) L-bar x6 | (L-bar,d) L x3 |
    (a,d-bar) nu^c x3 | (b,d-bar) e^c x3     -> 3 families + 6 Higgs pairs + 3 nu^c.
- Realisation on X0(143): every odd degree from W13-fixed CM points (W11 breaks, as the parity theorem
  requires; W13 survives), e.g. L_c = O, L_d = O, L_L = O(−P1−P2−P3), L_a, L_b = O(P1+P2+P3) up to degree-0
  twists (Wilson lines).  Purity with S0 (from u = ±i/sqrt 13, two of each sign): the degree ±3 sectors (Q, u^c,
  d^c, L from d, nu^c, e^c) have h^1 = 0 exactly; the degree ±6 sectors likewise; index-0 sectors with trivial
  difference bundle carry h^0 = h^1 = 2 extra pairs unless separated by a generic Wilson line (then 0).

## Step 3 — anomalies (EXACT)
For M1: Y is anomaly-free (SU(3)^2 Y, SU(2)^2 Y, grav-Y all 0) and lies in the 3-dimensional anomaly-free
span of the five stack U(1)s; the two remaining combinations are anomalous and must be lifted by the parent's
Green–Schwarz sector (a declared parent input).  Stack ledger (SU3^2, SU2^2, grav): c (0, 9/2, 0),
L (−3, −12, −48), a (3/2, 3, 24), b (3/2, 3, 24), d (0, 3/2, 0).

## Step 4 — Higgs sectors (skeleton, exact where stated)
The six doublet pairs sit in the (L, a-bar) and (L, b-bar) sectors with opposite hypercharge; their mass
requires a Y = 0 VEV in the (a, b-bar)-type singlet sector, which every one of the 30 solutions possesses.
Electroweak breaking is condensation in the (L, x-bar) sectors (bundle recombination, `surface.condensation`)
with the post-Higgs spectrum given by cohomology of the recombined bundle; the identification of the surviving
doublet and the Higgs mass are step-5 quantities.

## Step 5 — specification (tool to build next)
Yukawas Q u^c H_u, Q d^c H_d, L e^c H_d, L nu^c H_u are triangle couplings of three bifundamental sectors
(automatically gauge-invariant); their values are cup products of sections of S0 (x) L_x L_y^{-1} — bundles of
the form O(6(0)+6(1/11) ± CM points), whose section spaces need Riemann–Roch bases (weight-2 forms divided by a
W13-even differential, with the vanishing at that differential's other zeros handled), then Petersson-type
normalisation with the compact metric factor.  This is new machinery of the same kind as `yukawa`/`petersson`.

## Caveats (recorded)
- Parent convention: multiplicities used chi(S0 (x) E) = deg E (a 6D Weyl fermion with internal spinor S0).  A
  supersymmetric parent counts chiral multiplets by h^1(L_x L_y^{-1}), shifting every index by −12; the same tool
  runs that search with the shifted formula — not done here.
- Green–Schwarz cancellation of the two anomalous U(1)s, the potential fixing which doublet pair survives, the
  KK scale, and the gauge couplings are parent inputs.  Nothing is fitted; nothing is a prediction yet.

## Correction appendix (2026-09-13, from Kimi's release audit)
- The net nu^c count across the 30 five-stack solutions ranges over {±1, ±3, ±5, ±7, ±9}, not {±1, ±3} as
  stated above; the JSON is authoritative.
- v0.29.1 changed two APIs (condensation_energy -> dict; petersson_gate keys) and left two existing tests
  stale (test_surface_yukawa.py, test_surface_petersson.py); Kimi repaired both in staging.  Process change:
  every candidate now gets the full fast suite and a clean-tree build (no __pycache__ in the sdist).
