# Finite spacetime covariance addendum

7 September 2026. Recorded before running `local_qed_covariance.py`.

Reuse the stored 13-complex-dimensional Hecke internal operators with
`M = M0 + 0.2 balanced_V`. Supply a flat periodic Euclidean lattice of size
3^4, four Dirac spin components, and the Wilson regulator with a = r = 1.
These spacetime and spin structures are assumptions of the experiment.

For seeds 2026090716 through 2026090719, generate real Gaussian link
potentials A, independent site phases chi, and complex Gaussian spinors.
Set e = sqrt(4 pi MTFT.GAUGE.alpha), U_mu(x) = exp(i e a A_mu(x) Q),
and Omega(x) = exp(i chi(x) Q). Test the supplied Wilson–Dirac operator
under U'_mu(x) = Omega(x) U_mu(x) Omega(x+mu)^dagger and
psi'(x) = Omega(x) psi(x).

Compare three fixed trial charge assignments on (ell, old, q4, q6):
(-1,-1,-1,-1), (0,-1,-1,2), and (0,-1,0,2). The first is implemented as
Q = -I exactly; the latter two use the frozen sector projectors. Keep M
fixed first. Gauge covariance should hold in the first two cases. In the
third, the defect must agree with [M,Omega(x)] psi(x). As a control,
transform M into the site-dependent spurion Omega M Omega^dagger; this
must restore covariance. It does not construct a charged Higgs field.

Record absolute norms and relative norms normalized by the larger of
||D'_U psi'||_2 and ||Omega D_U psi||_2, using the Euclidean norm over
all sites, spins, and internal components. Defect-agreement errors use
the same denominator; additionally record their relative error against
the defect when it is nonzero. Numerical acceptance is 1e-10 for
Clifford, unitarity, gamma5-adjoint, commuting-case covariance,
defect-identity, and spurion-covariance checks. Require the distinct-charge
case to have a relative covariance defect greater than 1e-4 for every
seed. This is a finite operator test, not a continuum limit, an electron
identification, or a Standard Model derivation.
