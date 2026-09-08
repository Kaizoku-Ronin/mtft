# MTFT electron–photon investigation

**7 September 2026 · MTFT 0.26.0 with two local audit corrections**

We now have an executable electromagnetic compatibility test bed for the MTFT internal operators. Its strongest MTFT-specific result is a charge selection rule: **the old–quartic mixing found in the previous Hecke experiments can enter a neutral mass term only if the mixed components have equal electric charges.** We verified that rule algebraically, through a tree-level gauge identity, and inside a finite spacetime fermion operator with deliberate failing controls.

This makes the internal geometry useful for constraining candidate interactions. The experiment does not identify an electron within that geometry or derive the Standard Model. A Maxwell field, Dirac spin, spacetime, charge assignments, and a dimensional mass anchor remain explicit inputs. The work therefore gives us a concrete way to test a proposed physical identification before accepting it.

## The action and its inputs

Use the conditional action template, in natural units,

\[
\mathcal L=-\frac14F_{\mu\nu}F^{\mu\nu}
+\bar\Psi\bigl[i\gamma^\mu(\partial_\mu+ieQA_\mu)-M\bigr]\Psi
+\mathcal L_{\mathrm{Pauli}}.
\]

Here Q is a Hermitian internal charge operator; the electron convention is Q = −1. M is a Hermitian, spacetime-independent mass operator for the neutral-mixing tests. The Pauli term is an optional gauge-invariant dipole interaction. Maxwell and Dirac dynamics are standard supplied structure; their spin and propagation properties do not follow merely from the number of internal dimensions. [Tong, Standard Model symmetries](https://davidtong.org/pdfs/teaching/standard-model/standardmodel1.pdf)

| Ingredient | Provenance in this experiment |
|---|---|
| Rational Hecke sector projectors | Reused frozen X0(143) data; exact projector identities rechecked |
| Positive Hodge metric and complex structure | Frozen MTFT data, with numerical compatibility checks |
| Trial internal mass/mixing operator | Chosen affine T2 baseline plus controlled perturbations |
| Coupling magnitude | Existing MTFT three-term alpha relation, interpreted conditionally as the low-energy electromagnetic coupling |
| Spacetime, local U(1), Dirac representation | Supplied physical structure |
| Sector charges | Trial assignments subject to consistency constraints |
| Electron rest energy for eV conversion | External CODATA input, 510998.95069 eV |
| Additional Pauli coefficient | Free diagnostic coefficient; no MTFT matching calculation supplied |

We tested two instances of the template: thirteen internal complex components for the structural covariance experiments, and a prescribed single electron for electromagnetic benchmarks. **The step that selects and normalizes that physical electron from the thirteen-dimensional internal space remains open.** These are compatible component tests, not a complete interacting continuum solution.

## What Hecke mixing permits

The real sector dimensions 2, 4, 8, 12 become complex dimensions 1, 2, 4, 6 after using the Hodge complex structure. These are internal dimensions, not particle or spin counts. All adjoints and norms are taken in the same Hodge metric, transformed to an orthonormal frame.

Restrict attention to four scalar sector charges,

\[
Q=q_{\rm ell}P_{\rm ell}+q_{\rm old}P_{\rm old}
+q_4P_4+q_6P_6.
\]

The projectors are exact rational idempotents, sum to the identity, and commute with T2. For any mass matrix,

\[
P_i[Q,M]P_j=(q_i-q_j)P_iMP_j.
\]

Thus each nonzero neutral mixing block equates the two corresponding charges. In the restricted four-charge family, the allowed real coefficient space has the following dimensions:

| Neutral mixing pattern | Charge constraint | Remaining dimension |
|---|---|---:|
| No mixing between sectors | None | 4 |
| Old–quartic mixing | q_old = q_4 | 3 |
| A connected mixing graph spanning all sectors | All four charges equal | 1 |

These dimensions describe a Lie-algebra parameter space. A globally compact U(1) also requires a compatible charge lattice and normalization; the spacetime examples use integer charges. The calculation does not classify the full operator commutant, derive charge quantization, or assign observed particles to sectors.

For numerical realization we chose

\[
M_0=I+\frac{T_2}{4\|T_2\|_{\rm op}},\qquad M=M_0+\epsilon V,
\]

with Hermitian V of complex Frobenius norm one. The old–quartic perturbation preserves the two Atkin–Lehner involutions. Its only off-diagonal sector blocks connect old and quartic components. Sixteen fixed seeds, four mixing strengths, and four charge differences give **256 scan points**. They obey

\[
\|[Q,M]\|_F=\epsilon\,|q_{\rm old}-q_4|
\]

to within 3.32 × 10⁻¹⁵. The minimum scanned trial mass eigenvalue is 0.587109 in the chosen dimensionless units. Positivity here belongs to this selected finite operator, not to a derived particle spectrum or a continuum mass-gap theorem.

The stripped tree vertex is Γ^μ = γ^μ ⊗ Q. With S⁻¹(p) = slash(p) ⊗ I − I ⊗ M, its Ward-identity defect is I ⊗ [M,Q]. For mixing strength 0.2, equal old/quartic charges give a residual below 7 × 10⁻¹⁵; changing their charge difference to one gives the expected residual **0.4**, including the spinor Frobenius factor. This failing control shows that the check can detect the forbidden interaction.

## The same constraint inside spacetime

We supplied a periodic Euclidean 3⁴ lattice, four Dirac spin components, thirteen internal components, and a Wilson regulator. The field has **4,212 complex entries**. Links use the same chosen coupling, e = √(4π alpha_MTFT) = 0.30282212059140423. Four seeds and three charge assignments give twelve configurations.

Under site transformations Ω(x), links transform as U′ = Ω U Ω† at the adjacent sites. The fixed-mass fermion covariance defect is exactly [M,Ω]Ψ. This is a finite operator test; it does not establish continuum Lorentz reconstruction or solve the quantum path integral.

| Case | Relative fermion covariance defect |
|---|---:|
| Common charge, or equal charges on mixed sectors | At most 8.21 × 10⁻¹⁶ |
| Different old/quartic charges | 0.008394–0.009443 |
| Difference between measured defect and [M,Ω]Ψ | At most 8.15 × 10⁻¹⁶ |
| Transform the mass source as M′ = ΩMΩ† | At most 9.08 × 10⁻¹⁶ |

The denominator is the larger norm of the two transformed fermion-operator outputs, not the small defect itself. Full absolute norms are retained in the data. Clifford, unitarity, and gamma5-adjoint checks also pass.

Making M transform is a **spurion control**: it shows what transformation law a charged source would need. It supplies neither a dynamical scalar nor a Higgs potential. It also does not justify mixing particles of different charge while leaving electromagnetism unbroken.

## Photon and atomic controls

For the supplied spatial Maxwell discretization, each sampled nonzero momentum has two transverse modes and one gauge direction. We examined three physical momenta on grids with 16, 32, 64, and 128 points per direction, using Fourier symbols rather than allocating every spatial site. The worst relative algebra residual is 9.21 × 10⁻¹⁶. Dispersion errors decrease by approximately four when spacing halves, approaching the assumed c = 1. Zero momentum is handled separately as a homogeneous zero-frequency sector.

The Coulomb calculation numerically diagonalizes the radial Hamiltonian for six states through principal quantum number n = 3, with orbital angular momenta 0, 1, 2. It assumes a point nucleus of infinite mass. At the finest step, 0.0125 Bohr radii:

| State | Numerical energy, Hartree | Continuum energy, Hartree |
|---|---:|---:|
| 1s | −0.499980470276 | −0.5 |
| 2s | −0.124998779321 | −0.125 |
| 2p | −0.125000406909 | −0.125 |
| 3s | −0.055555314431 | −1/18 |
| 3p | −0.055555743103 | −1/18 |
| 3d | −0.055555571631 | −1/18 |

All six states show approximately second-order convergence. The apparent splitting of equal-n states in this table is a discretization effect, not predicted fine structure. The ground-state relative continuum error is 3.906 × 10⁻⁵. Increasing the boundary radius from 100 to 150 Bohr radii changes any reported energy by at most 2.16 × 10⁻¹² Hartree. An independent closed-form solution of the infinite discrete radial grid agrees with the ground-state solver within 3.63 × 10⁻¹³ Hartree.

Using the MTFT three-term alpha and the external electron mass, the analytic leading Coulomb energy is **−13.6056930910 eV**. This is the infinite-nuclear-mass nonrelativistic value. Recoil, relativity, Lamb shift, nuclear structure, and hyperfine effects are excluded; those are required for precision spectroscopy. [Tong, Atoms](https://davidtong.org/pdfs/teaching/topics-in-quantum-mechanics/topicsinqm3.pdf)

## Magnetic response is an additional matching problem

Define the extra Pauli coefficient κ_P by its contribution to the zero-momentum magnetic form factor. At the retained order,

\[
F_1(0)=1,\qquad F_2(0)=\frac{\alpha}{2\pi}+\kappa_P,
\qquad |g|=2[1+F_2(0)].
\]

This combines the standard one-loop minimal-QED contribution with a tree-level Pauli contribution. Renormalized mass and charge are held fixed. The Pauli vertex is transverse because q_μ σ^(μν) q_ν = 0, which we verified exactly in the supplied Dirac representation and numerically on 32 momenta. Consequently, the Ward check cannot fix κ_P. The separation of charge and magnetic form factors is standard in electromagnetic effective theory. [Paz, An Introduction to NRQED](https://arxiv.org/pdf/1503.07216)

| Added κ_P | Magnetic g magnitude at retained order | Mass and charge |
|---:|---:|---|
| −0.01 | 1.98232281946 | Fixed |
| 0 | 2.00232281946 | Fixed |
| +0.01 | 2.02232281946 | Fixed |

These are deliberately contrasting models, not experimentally allowed fits. The leading spin-independent Coulomb energy also stays fixed, while magnetic response changes. Higher-order atomic effects can change. Minimal QED fixes its own radiative terms; a more general effective theory needs its additional dipole coefficient matched to an underlying model or measurement. Terms from Pauli insertions in loops, including order alpha × κ_P, are omitted here. No full electron g−2 prediction is claimed.

For MTFT, this identifies a productive target: derive a definite electromagnetic vertex or dipole coefficient from a specified interaction construction. Mass and charge agreement alone cannot complete that task.

## Existing alpha formulae remain distinct hypotheses

We propagated every existing option without selecting or fitting a new formula. The reference is the recommended CODATA 2022 inverse alpha, 137.035999177 ± 0.000000021. [NIST constants table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt)

| Coupling input | Inverse alpha | Offset / reference uncertainty | Leading binding-energy relative shift |
|---|---:|---:|---:|
| MTFT three-term, current default | 137.035999337732 | +7.65 | −2.346 × 10⁻⁹ |
| MTFT four-term | 137.035999165439 | −0.55 | +1.687 × 10⁻¹⁰ |
| MTFT Monster | 137.035514275416 | −23090.55 | +7.077 × 10⁻⁶ |

The uncertainty column is an offset against the experimental reference only. It is not a significance from a complete MTFT likelihood: no theory-error model has been derived, and the four-term expression includes the package's specified numerical confinement-depth input. Hydrogen and magnetic quantities obtained by reusing alpha are correlated consequences, not independent confirmations.

The finest radial grid error is much larger than the three-term coupling's parts-per-billion energy shift. That shift is an analytic rescaling result; the finite-grid spectrum cannot resolve it. The package's rounded electron catalog entry also differs from the reference mass by about 2.05 ppm. We therefore keep catalog precision, numerical precision, and physical accuracy separate.

## Corrections and reproducibility

Two local source corrections accompany the study:

- **SM-AUDIT-01:** public Higgs quartic changed from 0.5179175466708401 to 0.12947938666771003, enforcing m_H² = 2λv². The historical value remains available under a retired accessor.
- **SM-AUDIT-02:** the random SU(N) initializer now divides its single corrected row by the full determinant. The historical defective initializer remains available for reproducing earlier runs. Existing hot-start trajectories made with it need rerunning before being interpreted as SU(N) simulations.

Ten targeted pytest cases pass. Forty-eight group draws across N = 2, 3, 5 satisfy the group-membership gates. Both implemented lattice action terms remain invariant under a local SU(3) transformation on a separate 2⁴ configuration, with relative differences below 2.1 × 10⁻¹⁶. This tests finite implementation consistency; it establishes no continuum mass gap. The full MTFT test suite was not run.

The bundle contains the pre-run plans, inputs, exact rational certificate, all sampled data, corrected source archive, unified patch, scripts, figure, source hashes, and an independent review record. The upstream archive and previous experimental reports are preserved. Exact statements concern finite algebra; all floating-point residuals are diagnostics rather than rigorous interval certificates.

The next construction to target is a specified geometric interaction that selects an electron subspace, a charge normalization, and an electromagnetic vertex together. We can now subject that construction to the charge, propagation, atomic, and magnetic checks developed here.
