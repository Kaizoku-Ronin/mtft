# mtft v0.23.0 — theory into tools

Four new modules packaging the 2026-08-27 E2-certified machinery.

## mtft.homology
Canonical integral H_1(X0(143), Z) from the independent PARI/GP
mslattice route, shipped as package data with every structural fact
re-asserted exactly on load (integral commuting involutions W11, W13,
STAR; unimodular antisymmetric intersection form; symplectic /
anti-symplectic conditions). Exact integer linear algebra throughout:
Fraction Gauss-Jordan inverse, Bareiss determinant, integer symplectic
Gram-Schmidt producing U in GL(26, Z) with U^T P U = J on the nose.
`periods_frame_ops()` gives the operators exactly in the frame of the
frozen Riemann matrix tau0 (documented: in that frame the integral
lifts commute with STAR only mod 2; the GP frame carries fully
commuting integral lifts). Ledgered trap: the naive Manin-lattice
saturation is an index-4200 sublattice with degenerate mod-2 pairing.

## mtft.thetachar
Frame-agnostic theta-characteristic / spin-structure machinery: GF(2)
core, the affine action t -> (M^{-1})^T t + d with d_j = q0(M^{-1}e_j)
(self-checked), Arf parity, exact fixed-locus and parity counting,
Burnside and full-subgroup-lattice orbit inversion with divisibility
asserts. Certified census (Cert EXACT, E2) shipped as package data;
`x0143_periods_frame()` and `x0143_gp_frame()` both reproduce it, and
`invariant_characteristics()` returns the 128 (96 even / 32 odd) as
explicit vectors in tau0's frame.

## mtft.thetafun
Genus-g theta functions with characteristics, built for g = 13 at tau0:
own LLL reduction, exact Fincke-Pohst ellipsoid enumeration with the
innermost two coordinates vectorized, a template-split evaluator
(Q = Q_out + ||A x_in + B x_out||^2) for large g, rigorous
theta-series product tail bounds with optimized splitting, and exact
characteristic transport through GL and shift moves with a CALIBRATED
unimodular phase (certified constant per reduction). Values and
gradients; `stats=True` returns npoints / absmass / radius for error
budgets. Validated: g=1 vs mpmath (incl. gradient), g=2 through
nontrivial (U, S) vs direct summation, split-vs-direct bit-agreement
at g=8, and a g=13 smoke on a non-invariant odd characteristic whose
null cancels to 3.1e-16 across 1.9e7 lattice points (58 s, tol 1e-3).

## mtft.liealg
The hardened Lie-fingerprint toolkit, generic over u(n): greedy closure
with the ABSOLUTE two-tier gate and ambiguity exception, structure
constants with reconstruction certificate, center / derived / rank /
Killing, representation decomposition and commutants, normalizer and
kernel dimensions by SPECTRAL-GAP cut (never fixed tolerance), Cartan
and roots with cosine-bucket diagnostics, and the arithmetic symmetry
screen. `d4_report()` reproduces the full CERT(tol, E2) fingerprint of
the STAR-fixed triangle algebra in ~3 s: dim 28, growth 3-6-17-28,
Killing (28,0,0), 8 + 1^5, normalizer 54 across a 1.9e7 gap, 24 D4
roots, STAR identity, W11/W13/W143 non-normalizing.

## Studies
`studies/TH2_PREREGISTRATION.md`: frozen decision rules for the
odd-invariant-gradient probe s(epsilon) against the certified D4
projector, two-route E2 protocol, error budget, even-branch null
census, G2-equivariance gate, exploratory tier. Filed before any
invariant characteristic was evaluated.

## Tests
`tests/test_v0230.py` (pytest or direct): structural homology gates,
census in both frames, theta engine at g = 1 / 2 / 8, D4 report
asserts. MTFT_SLOW=1 adds the full census recomputation.
