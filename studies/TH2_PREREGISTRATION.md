# TH-2 Preregistration — Odd invariant theta gradients vs the D4 active 8

Filed: 2026-08-27, before any numerical evaluation of the 128 invariant
characteristics.  Package: mtft v0.23.0 (modules thetachar, thetafun,
liealg, homology).  Status of all inputs: theta census Cert (EXACT, E2);
D4 fingerprint and 8 + 1^5 decomposition CERT(tol, E2); frozen Riemann
matrix tau0 certified (19 PASS / 1 honest negative C9b, external queue).

No quantity defined below has been computed at filing time.  The engine
was smoke-tested exclusively on NON-invariant characteristics (ledgered:
one odd non-invariant null = 3.1e-16 through 1.9e7 lattice points).

## 1. Objects

G2 = <W11, W13, STAR> acts on theta characteristics in the periods
symplectic frame of tau0 (mtft.thetachar.x0143_periods_frame).  Its joint
fixed locus is an affine F_2^7: 128 characteristics, 96 even, 32 odd
(Cert EXACT, E2).  For each of the 32 ODD invariant characteristics
epsilon, theta[epsilon](0, tau0) = 0 identically and the gradient

    grad(epsilon) = grad_z theta[epsilon](z, tau0) |_{z=0}  in  C^13

is the canonical spin-Dirac datum of the corresponding invariant theta
divisor through the origin.  P8 denotes the certified projector onto the
active 8 of the D4 algebra g_fix (mtft.liealg.rep_summary P_active;
invariance leak < 1e-8 asserted at runtime), P5 = 1 - P8 the projector
onto the five singlets.

Primary observable, per odd invariant epsilon:

    s(epsilon) = || P5 grad(epsilon) || / || grad(epsilon) ||.

s is invariant under the global phase and scale of grad, hence under all
transport phases of the evaluation pipeline.

## 2. Frozen decision rules

Compute s(epsilon) for all 32 odd invariant characteristics with the
certified error intervals of Section 4.  Let s_hi(epsilon), s_lo(epsilon)
be the certified upper/lower bounds.

  H-A (gradients live in the active 8):
      PASS iff max over the 32 of s_hi(epsilon) < 0.01.
  H-B (gradients live in the singlet 5):
      PASS iff min over the 32 of s_lo(epsilon) > 0.99.
  H-C (mixed): otherwise.  The full distribution {s(epsilon)} is
      reported with intervals; no post-hoc threshold may be introduced.

Exactly one of H-A / H-B / H-C is recorded as the verdict.  A mixed
verdict is a finding, not a failure; under H-C the following
PREREGISTERED secondary summaries are reported (no others):
  (i)  the multiset {s(epsilon)} to 3 decimals;
  (ii) the rank of the 13 x 32 gradient matrix by spectral gap
       (mtft.liealg.spectral_gap_kernel on its singular values);
  (iii) the spectrum of pairwise |cos| between normalized gradients.

## 3. E2 protocol (two computational routes, disjoint steps)

Route A: direct evaluation at tau0 via mtft.thetafun.theta_grad
(LLL-reduced frame, template-split enumeration, calibrated transport).

Route B: evaluation in a preregistered twisted frame.  Let

    V = I_13 + E_{1,2} + E_{3,7} - E_{9,4} + E_{13,10} + E_{2,11}

(unit diagonal, unimodular; E_{i,j} the matrix unit, 1-indexed).  Set
tau_V = V^T tau0 V.  Each characteristic is transported through the
V-move by mtft.thetafun.reduce_char applied with a synthetic ready
structure (U = V, S = 0), then evaluated by a FRESH siegel_ready(tau_V)
reduction.  Route B therefore sums over a genuinely different lattice
presentation with different enumeration geometry, different calibrated
phases, and a different LLL basis.  Gradients pull back through V; the
observables compared are phase-free:

    |theta| values (even branch) and s(epsilon) (odd branch).

Agreement requirement: |s_A - s_B| <= (interval radii sum) for every
epsilon; any violation voids the affected epsilon and is reported as an
engine discrepancy, not adjudicated post hoc.

## 4. Error budget

Certified: the truncation tail bounds returned by the engine
(tail_bound_value / tail_bound_grad; conservative theta-series product
bounds; see module docstrings).  Reported alongside, NOT certified:
floating-point roundoff estimate r = 20 * eps_machine * absmass with
absmass = sum of |terms| (engine ``stats`` output).  Interval radii:

    delta_g = grad_bound + r_grad     (per gradient evaluation)
    s_hi = (||P5 g|| + delta_g) / max(||g|| - delta_g, 0)
    s_lo = (||P5 g|| - delta_g) / (||g|| + delta_g)   (clipped to [0,1])

Tolerance ladder (frozen): gradients at tol = 1e-4; if ||g|| < 100 *
delta_g, re-run that epsilon at tol = 1e-6; if still ||g|| < 100 *
delta_g, mark epsilon INDETERMINATE and report it as such (it counts
against H-A and H-B, i.e. both require all 32 determinate).

## 5. Even branch (same run, preregistered)

For all 96 EVEN invariant characteristics, compute theta[epsilon](0,
tau0) at tol = 1e-6 on both routes.  epsilon is declared a NUMERICAL
NULL iff |theta| <= 10 * (val_bound + r_val) on BOTH routes.  Report the
count N_van and the separation ratio

    min over non-null |theta|  /  max over null |theta|.

The count is asserted as a finding only if the separation exceeds 100;
otherwise the values are tabulated without a vanishing claim.  A nonzero
N_van (effective even characteristics) is geometric information about
semicanonical divisors of X0(143) and feeds R7.

## 6. Internal consistency checks (run before the primary computation)

(a) G2-equivariance / convention certificate.  The AL and STAR
operators fix tau0's frame; for gamma in G2 the classical transformation
law forces |theta[gamma . epsilon](0, tau0)| = |theta[epsilon](0, tau0)|.
Preregistered samples: the G2-orbits of the three NON-invariant seeds
(periods-frame coordinates, 26-bit vectors with 1s at positions listed,
0-indexed): s1 = {0}, s2 = {1, 20}, s3 = {4, 17}.  For each orbit,
max |theta| deviation across the orbit must be <= 10 * (val_bound +
r_val) at tol = 1e-8.  Failure falsifies the a/b frame convention and
BLOCKS the primary run (fix, re-file an amendment, then proceed).

(b) Projector certificate: rep_summary invariance_leak < 1e-8 and
active_dim = 8, common_fixed_dim = 5, re-asserted in the run.

## 7. Exploratory tier (labeled, no decision weight)

For each odd invariant epsilon and each op O in {W11, W13, STAR} (in
their complexified adapted-frame forms; STAR antilinear), report the
alignment |<U_O g, g>| / ||g||^2.  These probe how the induced action
organizes the 32 gradients; no thresholds, no verdicts.  The anomalously
small W143 non-normalization residual (0.196 / 0.327 across conventions)
remains a flagged follow-up; any structure here feeds it.

## 8. Compute budget and reporting

Measured baseline (this container): one genus-13 value at tol 1e-3 =
58 s / 1.9e7 points.  Projected: gradient at 1e-4 in minutes per
characteristic; full TH-2 (32 odd x 2 routes + 96 even x 2 routes +
orbit checks) is an overnight batch on the production machine.  All
raw values, gradients, intervals, npoints, absmass, radii, and timings
are written to a JSON ledger; the verdict section quotes only
preregistered quantities.  Deviations from this document, if any, are
recorded in an appended DEVIATIONS section before results are read.

Queued follow-ups (not part of TH-2): W143 structural probe; DSZ-type
pairing scan on the invariant set; R7 exact-catalog cross-reference.
