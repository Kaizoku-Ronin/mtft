# W143 coupled-mode reduction and complementary-sector feedback

Frozen 8 September 2026 before reading the new canonical-mode and response
observables. Source is the unchanged supplied MTFT0.26.2. Inputs are the
previous direct and localization projectors and the Hodge-orthonormal frame.
This is finite operator mathematics, not a physical Hamiltonian identification.

## Objects and constructions

Use W143 freshly constructed from the package and transported through the
previous frame. Verify selfadjointness, involutivity, and complex linearity
before treating it as a complex13 matrix. The direct active projector has
complex rank8, complementary rank5. Write W in orthonormal sector bases as
[[A,B],[B*,D]]. No triangle, operator, sector, or signature search is allowed.

Primary canonical route: SVD of B, retaining singular values normalized by
||W||_F at or above1e-7, rejecting those at or below1e-9, with an ambiguous
verdict between the bands. Construct the corresponding coupled2x2 blocks,
and diagonalize the uncoupled restrictions. Do not assume their signs.
Record all singular values, a_j,b_j,d_j, off-block residuals and dimensions.

Independent route: use the earlier localization projector and the positive
eigenspace of R_+=(I+W)/2. Diagonalize its compression of P. Principal overlaps
rho_j determine a_j=2rho_j-1 and b_j=2sqrt(rho_j(1-rho_j)). Compare complete
canonical spectra and reconstruction residuals, allowing phase/unitary freedom
inside repeated or decoupled eigenspaces. Shared upstream inputs do not count
as independent physical data.

## Effective operator and regular response grid

For z outside spec(D), define Sigma(z)=B(zI-D)^-1 B* and
H_eff(z)=A+Sigma(z). Compare

  G_full(z)=U*(zI-W)^-1 U,
  G_eff(z)=(zI-A-Sigma(z))^-1,
  G_involution(z)=(zI+A)/(z^2-1).

The final identity is valid for z!=+-1 even when zI-D is singular; the direct
Schur expression additionally excludes spec(D). The baseline that drops
feedback is G_naive(z)=(zI-A)^-1. Record its error without optimizing A or z.

Regular grid: x=linspace(-1.5,1.5,301), z=x+i eta for eta in{0.05,0.2,0.7},
plus z in{-2,0,2,i}. All907 points are fixed. Compare matrix norms, resolvent
errors, Sigma ranks, and condition numbers. Numerical verification gate:
relative full/effective and full/involution errors <=1e-9 on the regular grid.
Record every point and summarize maxima. No claim of a global certified bound.

## Registered pole-conditioning stress test

For each eigenvalue d_j of D whose distance from both+-1 exceeds1e-7, evaluate
z=d_j+i epsilon with epsilon in{1e-2,1e-3,...,1e-12}. This rule selects all
interior compression poles, without inspecting outcomes. Record norms and
condition numbers of the Schur denominator, Sigma, and the continued/full
resolvents. These conditioning diagnostics do not alter the regular-grid gate.

At real z=d_j, evaluate only the full and involution formulas and their
canonical channel responses. Do not numerically invert zI-D there. Separate
true full-system poles at+-1 from poles of the feedback representation and
zeros of the projected response. Never interpret a compression pole as a
new physical eigenvalue.

## Optional-clock mathematical control, included in this run

As an explicitly chosen dimensionless unitary example, set U(t)=exp(-itW).
This t has no identified physical units or MTFT clock interpretation.
Use t=linspace(0,2pi,257) for the analytic response curve and independent
matrix-exponential checks at t in{0,0.17,pi/4,pi/2,pi,2pi}.
For initially active coupled modes and zero complementary input, compare
projected propagation with cos(t)I-i sin(t)A and complementary norm fractions
b_j^2 sin(t)^2. Test the matrix identity
I-U_PP(t)* U_PP(t)=sin(t)^2 B B*.
The active-only exp(-itA) is a baseline, not the claimed reduced evolution.

Derive the complementary-sector memory kernel B exp(-itD) B* symbolically;
no long numerical Volterra simulation is required. State initial-state
assumptions. These calculations test exact finite-system elimination and do
not establish a QFT, spacetime signature, physical transition rate, or mass.

## Reporting

All period-dependent coefficients remain float64 DIAGNOSTIC. Exact identities
and the involution rank bound are algebraic statements under their hypotheses.
An independent symbolic check and independent principal-overlap reconstruction
must precede the final report. Preserve any numerical stress failures and
their conditioning context. No package Lie-closure gate is used or repaired.
Save the protocol, scripts, all outputs, frozen inputs, source hashes, scientific
figure, and reproducible runner. Label unregistered follow-ups separately.
