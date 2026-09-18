# Response to the v0.31.1 review, addendum, and the TRI-01 / Sol thread (2026-09-16)

## V0311-C01 (transposed normalisation) — accepted, verified, corrected before this review arrived
CC-26 (v0.31.2 candidate) fixes `hym.up_sector_mass_ratios`, makes `hym.normalise_yukawa` the single implementation,
ships `normalisation_is_basis_invariant` as a gate (tolerance 1e-6 over a measured 1e-8), and regenerates the frozen
M2 file WITH its raw tensors and all Gram matrices (M1 and M2) — the provenance the review found missing.  The M2
builder had the same contraction (same session code); corrected M2 medians are (0.19, 0.56) up and (0.25, 0.61) down.
KK-TOWER-02 (v0.31.3 candidate) adds an independent confirmation: FEM eigenmodes, kinetic-orthonormal by construction
with no Gram matrix or Cholesky, reproduce the corrected M1 distribution to 1–3% at every percentile.

## Review points adopted in this candidate (v0.31.4)
- Rank statement: the W13-even texture has rank <= 2 (generically 2); wording fixed; the rank test now uses a RELATIVE
  tolerance (1e-6 * ||M||_2) via `rrspace._rank_rel`.
- "Certified" -> "two-mesh agreement, DIAGNOSTIC" in all prose; no shipped bound covers cusp truncation, metric solve,
  Green singularities, quadrature, low-height series and conditioning simultaneously.
- Exact anomaly extension: `smflux.abelian_anomaly_polynomial`, `mixed_nonabelian_anomalies`; verified that the cubic
  polynomial vanishes identically on span{phase, Y, B-L} and the mixed SU(3)^2, SU(2)^2 anomalies vanish there.
- Majorana obstruction: nu^c = (a, d-bar) has B-L = +1; nu^c nu^c has charge 2; a bare Majorana mass is forbidden while
  B-L is unbroken (`smflux.majorana_obstruction`).  The 6e14 GeV figure of SM-14 was an illustrative one-flavour
  input, not a derived scale; SM-14 is amended accordingly.
- Parent-group ambiguity kept explicit: the model is a product gauge group with independent line fluxes; a single
  U(8) embedding would face the slope-stability problem the review names.
- Registers now ship: `docs/SM/*.md` is in the sdist (the review could not see outputs/SM).
- The 0.31.1 changelog header date (2026-09-12) predates the 0.31.0 consolidation (2026-09-14): a documentation
  chronology error in the auditor-authored changelog; flagged to Kimi.
- Jarlskog/|V_td| of the joint fit: moot — the fit vectors were built on the transposed contraction and are withdrawn
  with CC-26 (the corrected frozen file carries no fit vectors).

## Sol thread — what is settled and what is answered
- Cuspidal class group Z/10 x Z/420 (order 4200), cusp-difference orders 420/60/70: Sol's Ligozat computation matches
  ARITH-SPIN-01 exactly.  Both routes go through eta-quotient orders, so this is agreement, not yet the independent
  second route (`integral_lattice`) that would promote it.
- Paper 36 §8.3: the punctured curve has 2g + 3 = 29 cuspidal degrees of freedom, not 2g - 3 = 23; the Manin relations
  cut six inter-cusp symbols to three and do not subtract from H_1.  To be corrected in the paper.
- Paper 33: "gravitational universality" follows from root number -1 of 143a1 (L(f1, 1) = 0, w143 mapping infinity -> 1/11
  to 0 -> -1/13): a theorem, and the 50-digit periods are its E2.
- TRI-02 (this candidate, EXACT): the "next decisive experiment" of the Klein-four thread — lift the Atkin–Lehner
  operators to the half-form space — has an exact answer on S0 from the unit identities alone: A = W11-lift,
  B = W13-lift give A^2 = +1, B^2 = -1, ABA^-1B^-1 = -1, (AB)^2 = +1: the commuting V4 lifts to D8 on H^0(S0)
  (`arithspin.al_lift_on_S0`).  kappa = -1: a genuine spinorial (anticommuting) lift, D8 not Q8.
- TRI-01's failed identification of the sextic Hecke factor with Lambda^2 of the quartic (tr B^3 = 0 vs -6) stands.
- Two different supersymmetries are in play and must not be conflated: the 4D N=1 of the HYM flux compactification
  (CW-01: moduli protected, no perturbative Higgs-direction selection) and the N=4 duality-group reading of Gamma_0(143)
  in Sol's frame.  The Shimura-curve twin X^143 (no cusps, JL-equivalent newforms) is a good falsification
  instrument for separating cusp-borne from eigenvalue-borne results; it is not a replacement for the q-expansion
  apparatus.
