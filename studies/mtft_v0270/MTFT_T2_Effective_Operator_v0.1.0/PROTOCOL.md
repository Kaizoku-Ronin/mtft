# T2: feedback without an involution identity

Frozen locally on 8 September 2026 before measuring the new T2 response,
compression poles, memory kernel, or time evolution. This follows the completed
W143 experiment. Earlier work already found T2 mixing the same active and
complementary sectors and filling the full space after one augmentation.
Those earlier observations are inputs, not newly blinded discoveries.

## Fixed source and geometry

Use the supplied, unchanged MTFT 0.26.2 source archive, SHA256
35cfc00c32877bd6c12b16464ddc8190bcca9e2f054474f4ab0a6c605fe9da38.
Construct fresh `mtft.hecke.cuspidal_hecke(2)`. Transport with the frozen R,Ri
from MTFT_Active_Module_v0.1.0. Use its direct complex active projector for the
primary route and localization real projector for the independent route.
No operator, triangle, sector, physical scale, or signature search is allowed.

Before spectral or dynamical interpretation, verify raw transported T2 is
selfadjoint and complex-linear, and P is an orthogonal complex projector of
rank 8 (complement 5). Relative Frobenius defects must be <=1e-9. Stop the
Hermitian experiment if these gates fail; do not repair a substantial defect.
After passing, an explicitly recorded anti-Hermitian roundoff cleanup is
allowed. Do not force an involution or rescale the arithmetic operator.

In orthonormal sector bases write H=[[A,B],[B*,D]], dimensions 8+5 complex.
All period-dependent quantities remain float64/complex128 DIAGNOSTIC.
The chosen generator H=T2 and the parameter t are dimensionless mathematical
controls, with no identification of physical time, energy, or Hamiltonian.

## Coupling structure and instantaneous dark directions

Record full spectrum, all five singular values of B, and complex rank.
Normalize singular values by ||H||_F: <=1e-9 small, >=1e-7 resolved, with an
ambiguous verdict in between. Never claim exact rank from these thresholds.
Construct SVD-paired active/complement directions; record the residual when
retaining only independent 2x2 pairs and any remaining active block. Measure
[A,BB*] and [D,B*B], and each candidate pair's invariance residual. This tests
whether W143's independent-pair structure survives; do not assume that it does.

For the active nullspace V0=ker(B*) obtained by the same SVD, record B*V0 and
B*A V0. The former tests instantaneous darkness; the latter tests delayed
coupling. Record rank decisions for B*A V0. Compute the complement reachable
span [B*, DB*, ..., D^4 B*]. Full column rank of B already makes all five
complementary directions necessary for an exact finite self-energy realization.
The general algebra and its hypotheses must be checked independently.

## Response and poles

Define Sigma(z)=B(zI-D)^-1 B*, H_eff=A+Sigma, and G_eff=(zI-A-Sigma)^-1.
Compare with G_full=U*(zI-H)^-1 U and a full-eigensystem spectral sum.
Baseline: G_naive=(zI-A)^-1. No W143 involution formula is presumed valid.

Fixed regular grid: x=linspace(-3,3,601), z=x+i eta with eta in {0.05,0.2,0.7},
plus z in {-4,4,2i,0.17+0.37i}: 1807 points total. If either real point is
actually spectral, report it and skip that inverse rather than substitute a
pseudoinverse. Registered response gate: relative Frobenius full/Schur and
full/spectral errors <=1e-9 on all admissible regular points. Report every
point, condition numbers, and baseline error without tuning the grid.

Diagonalize D. For each eigenvalue d_j, record the positive semidefinite
self-energy residue R_j=B v_j v_j* B*. Group degeneracies within
1e-9*max(1,||H||_2), using whole spectral-projector residues rather than
basis-dependent rank-one claims at a repeated eigenvalue. Check sum R_j=BB*.

Compression-pole stress rule: select every distinct d_j whose distance from
spec(H) exceeds 1e-7*max(1,||H||_2). Approach with z=d_j+i epsilon, epsilon
in {1e-2,...,1e-12}. Compare direct Schur, full solve, and full spectral sum;
record ill-conditioning and error without applying the regular-grid gate to
this stress test. At real d_j evaluate only the full solve and spectral sum.
Test the null relation G_PP(d_j) B v_j=0, normalized by ||G_PP|| ||B v_j||.
For repeated eigenvalues test the whole residue range. Do not invert zI-D at
these points. Compression poles that coincide with full eigenvalues must be
reported separately, without asserting a regular projected zero there.

## Memory and a dimensionless dynamical control

Derive K(t)=B exp(-itD) B* and compare it with sum exp(-it d_j) R_j.
Retain the forcing B exp(-itD)y(0) for nonzero initial complementary input.
The projected evolution is U_PP(t)=U* exp(-itH) U, not exp(-itA).

For t=linspace(0,2pi,257), use the full spectral decomposition to calculate
C(t)=F* exp(-itH) U. Record average complementary squared norm
||C||_F^2/8 and worst-input value ||C||_2^2. Also record average and worst
values restricted to the registered instantaneous dark subspace V0, normalized
by dim(V0) for the average. Report maxima only as sampled-grid maxima.
No fitted clock, optimized input search beyond this stated singular-value
definition, or asymptotic irreversibility claim is allowed.

Independent scipy matrix-exponential controls at t in
{0,0.17,pi/4,pi/2,pi,2pi}: compare full spectral evolution, projected response,
unitarity, I-U_PP*U_PP=C*C, and the memory-kernel spectral expansion.
Record the active-only exponential error. Numerical control tolerance 1e-9.

For short times t=0.1/2^k, k=0,...,5, check C(t)/t approaches -i B*.
For dark inputs check C(t)V0/t^2 approaches -B*A V0/2. Report normalized
errors and ||C V0||_F^2/(dim(V0)*t^4) against its derived coefficient.
Use the analytic expansion, not a fitted log-log exponent. Record the
corresponding ordinary active average divided by t^2. The generic dark
squared-norm expansion has an O(t^5) remainder unless extra symmetry is proved.

## Independent verification and reporting

An independent numerical route uses the full real 26-dimensional Hodge frame
and localization projector, rather than the primary complex block inverse.
It obtains full response from the real eigensystem and exact characteristic
polynomial from the arithmetic matrix, and checks the same regular grid.
Compare ambient projected responses, spectral values, and coupling singular
values. Shared upstream geometry is explicit, not independent physical data.

A separate symbolic control on a generic non-involutive rational Hermitian
example must check Schur reduction, pole cancellation/null direction, memory
and composition identities, and delayed leakage of dark inputs.

Preserve the frozen protocol, source and input hashes, all raw results,
scripts, independent checks, proof, scientific figure, and replay runner.
Report failed/ambiguous outcomes and stress loss. No new physical QFT,
particle assignment, or arithmetic certificate follows from a small residual.
