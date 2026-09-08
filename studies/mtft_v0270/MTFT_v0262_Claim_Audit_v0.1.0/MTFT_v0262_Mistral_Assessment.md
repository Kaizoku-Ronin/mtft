# MTFT 0.26.2 and the Mistral review

**7 September 2026 · audit of the attached release and the supplied review text**

Mistral points toward useful operator, theta, and statistical questions. Its central claim that D4 triality can identify the Clifford signature from an eight-dimensional sector is not justified. Several other suggestions need corrections before they become sound experiments. We checked those claims against the attached 0.26.2 source and ran bounded numerical and exact controls.

The full atlas reviewed by Mistral was not attached to this message. This report assesses the quoted claims and available code; it does not certify the atlas's theorem statements or its claimed completeness.

## Release changes and a correction to our earlier wording

Release 0.26.2 incorporates the SU(N) initializer fix as **CC-21**. A row multiplier changes a determinant by its full value, so dividing one row by det(Q), rather than its Nth root, is the correct operation. The release preserves the defective draw under `random_su_n_defective_v0261`.

For the Higgs quartic, **CC-22** clarifies two conventions and preserves the original API:

| Accessor | Value | Mass relation |
|---|---:|---|
| `HIGGS.lambda_quartic` | 0.5179175466708401 | m_H² = lambda v²/2 |
| `HIGGS.lambda_quartic_pdg` | 0.12947938666771003 | m_H² = 2 lambda v² |

The latter agrees with the default Hosotani self-coupling. I retract the earlier characterization of the larger number as intrinsically wrong: the demonstrated issue was a convention mismatch between APIs. Our earlier local patch adopted the second convention for the first accessor; the release instead preserves the formula and names the conversion. The numerical electromagnetic experiments are unaffected. `CORRECTION_TO_PRIOR_HIGGS_AUDIT.md` records this appended clarification.

Three targeted release tests pass, covering CC-21, CC-22, and four-way version consistency. An additional 24 sampler draws at N=2,3,5 pass group-membership checks. The full package test suite was not run.

Release 0.26.1, included here, also added the fast exact Ising density-of-states route C. Therefore the current 56-spin N=143 partition function is available through exact integer coefficients; the genus-seven limit concerns the completed cross-check of the much more expensive Pfaffian route. It is not a genus-seven ceiling on route C.

## 1. D4 is useful evidence, but triality does not choose spacetime signature

We reproduced `liealg.d4_report()`:

| Diagnostic | Reproduced value |
|---|---:|
| Lie dimension / rank | 28 / 4 |
| Closure growth | 3, 6, 17, 28, 28 |
| Killing-form inertia, negative/zero/positive | 28 / 0 / 0 |
| Complex module decomposition | active 8 + five singlets |
| Roots | 24, equal length |
| Root-angle cosine deviation | 6.61 × 10⁻¹⁶ |
| Largest rejected absolute closure residual | 9.26 × 10⁻⁸ |
| Relative bracket reconstruction residual | 1.07 × 10⁻⁵ |

These values are a replay of the package's numerical fingerprint, with its registered two-tier gate. They are not an exact algebraic reconstruction or an interval proof. In particular, the tiny root-angle error must not be substituted for the larger closure residual.

The construction operates in u(13), using anti-Hermitian matrices. Its negative Killing form supports the **compact real form so(8)** within that construction. Compactness here is compatible with the chosen adjoint structure; it is not a discovery of Lorentzian spacetime.

Compact Spin(8) has three distinct eight-dimensional representations: the vector representation and two half-spin representations. Triality permutes them within the same compact setting. An eight-dimensional support count cannot distinguish those labels, and selecting one label does not select a spacetime signature. [Baez, Spinors and Trialities](https://math.ucr.edu/home/baez/octonions/node7.html)

We supplied an exact control: a rational 4×4 matrix of order three preserves the positive Euclidean D4 root metric and all 24 roots while cycling the three eight-weight sets. This demonstrates the issue without relying on the numerical MTFT fingerprint. It is a standard representation-theory control, not an arithmetic triality operator constructed from MTFT.

Mistral's displayed Clifford isomorphism also fails a dimension check. The even algebra on four generators has real dimension 2³=8; Mat₂(H) has real dimension 16. Under the convention e_i²=−1, the full four-generator algebra is Mat₂(H), while the even algebra is H⊕H. A signature with five generators is a different dimensional question. [Baez, Clifford algebras](https://math.ucr.edu/home/baez/octonions/node6.html)

The description “Hecke-invariant subalgebra” is also unsupported. The object is the **STAR-fixed triangle algebra**. Its own symmetry screen reports Atkin–Lehner non-normalization residuals of approximately 0.938666, 0.938595, and 0.327152 for W11, W13, and W143. Working in Hecke-adapted coordinates does not establish invariance under those operators.

There is a useful conditional obstruction. If an exact nonzero fixed subspace F is annihilated by the Lie algebra, every algebra element squares to zero on F. The unital associative algebra it generates acts scalarly on F, so two invertible elements cannot anticommute there. Thus a construction confined to that algebra cannot supply a nondegenerate Clifford pair on the full module. The implication is exact linear algebra; the fixed-five decomposition currently has numerical evidence. Restricting to the active block or adjoining other operators changes the problem.

## 2. TH2 is a support test, and no result is packaged

The attached release contains `studies/TH2_PREREGISTRATION.md`. It says the primary quantities had not been computed at filing. No packaged TH2 gradient ledger or H-A/H-B/H-C result was found. We independently reproduced its finite input census in both provided frames: affine dimension 7, 128 invariant characteristics, 96 even and **32 odd**. We did not evaluate theta functions.

If H-A passes, it would show that the evaluated gradients meet the registered active-subspace support condition. It would not establish equivariance of the gradient construction, distinguish vector from half-spin labels, construct gamma operators, or identify their real quadratic form.

Before running it, the gradient covectors and their dual metric must be transported from the period frame into the projector's adapted frame. Phase/scale invariance alone does not make a ratio independent of arbitrary coordinate changes. The preregistration also labels its roundoff estimate **not certified**; any interval-style conclusion inherits that limitation even if its truncation-tail component is bounded.

## 3. Vanishing one-forms closes a channel, not all spectral geometry

The bimodule instrument changes the represented algebra through an operator W, but computes ordinary commutators [D,pi(a)]. It is not using a Connes–Moscovici twisted commutator. When all represented commutators vanish, this construction has no represented internal one-form fluctuations. That does not prove that the Atkin–Lehner operators are gauge transformations or that their physical connections have zero curvature. Inner fluctuations require the represented one-forms themselves. [Connes–Marcolli, spectral triples and inner fluctuations](https://www.math.fsu.edu/~marcolli/bookhyper.pdf)

The statement Aut(K4)=Aut(K6)={1} constrains automorphisms of those specified rational fields. It does not rule out every larger algebra, multiplicity-space transformation, scalar extension, or alternative bimodule.

We tested that distinction constructively. In the synthetic algebra A=R³ of diagonal matrices, take W to be a three-cycle and M=W+Wᵀ. In the doubled representation, exact integer checks satisfy the tested order-zero, first-order, self-adjointness, reality, and grading conditions, while each elementary generator has squared one-form commutator norm **4**. Identity and transposition controls give zero.

This is a counterexample to an overbroad “all such routes are trivial” conclusion. It is explicitly **not an X0(143) construction**, nor a verification of every spectral-triple axiom. Finding an arithmetic analogue remains an open construction problem.

## 4. Finite Ising data do not determine critical exponents

For N=143, route C reproduces 85 exact density-of-states coefficients with sum 2⁵⁶, and an even-subgraph polynomial with coefficient sum 2²⁹. Its elimination width is 11. The reference evaluation A(1/√3)=14.751840885839925 is a floating diagnostic of that exact polynomial.

The Laplacian statement about simple eigenvalues 0,1,2,4,5 is correct; we checked the characteristic polynomial exactly. The smallest positive eigenvalue is approximately 0.2726371273. The proposed inference to critical exponents is not justified: a finite Ising partition sum is analytic at real finite inverse temperature. A thermodynamic critical exponent requires a specified limit. [Tong, Phase Transitions](https://davidtong.org/pdfs/teaching/statistical-physics/statphys5.pdf)

For the finite heat trace K(t)=Tr exp(−tL), we measured

\[
d_{\rm eff}(t)=\frac{2t\sum_k\lambda_k e^{-t\lambda_k}}
{\sum_k e^{-t\lambda_k}}.
\]

| t | Including the zero mode | Omitting the zero mode |
|---:|---:|---:|
| 0.1 | 0.53223 | 0.54513 |
| 1 | 1.90330 | 2.11996 |
| 10 | 0.85021 | 6.44872 |
| 100 | 8.25 × 10⁻¹¹ | 54.78541 |

Including zero, the diagnostic tends to zero at both ends on this finite connected graph. Removing zero changes the large-time behavior to approximately 2 lambda_gap t. Neither operation extracts a universal dimension from a finite plateau. The finite positive-spectrum zeta sum is entire; it has no continuum dimension pole.

The most decisive topology control is exact: **N=6 and N=11 have identical dual adjacency matrices and identical Ising density of states [2,2,2,6,4,0,0], but their surface genera are 0 and 1.** We also checked their density of states independently by brute force. Hence their scalar zero-field thermodynamics agrees at every temperature, despite the different genera. Genus-dependent spin-structure decompositions are additional information, not something the scalar sum must recover.

The productive version of the proposal is finite diffusion diagnostics now, followed by a preregistered growing/refining graph family with specified metric and coupling scaling. Spectral dimension alone is not a universal formula for all Ising exponents.

## 5. The Mellin concern is right; the quoted equation is not

For the definitions used in the package, absolute convergence gives, for Re(s)>1,

\[
Z_D(s)=\frac{(2\pi)^s}{\Gamma(s)}
\int_0^\infty y^{s-1}Z_L(y)\,dy
=-\zeta(s)\zeta'(s+1).
\]

The corresponding inverse transform is

\[
Z_L(y)=\frac{1}{2\pi i}\int_{c-i\infty}^{c+i\infty}
\Gamma(s)Z_D(s)(2\pi y)^{-s}\,ds,\qquad c>1.
\]

Mistral combines the forward transform with the inverse transform's contour and prefactor, omitting the normalization required by exp(−2πyn). These are Mellin integral-transform relations. They do not alone establish a spacetime Wick rotation or analytic continuation of quantum-field correlation functions. [NIST DLMF, Mellin transform and inversion](https://dlmf.nist.gov/2.5)

Three finite-weight checks at 50-digit working precision, including one complex s, reproduce the correctly normalized transform with relative error below 3.32 × 10⁻⁵¹. Those quadrature checks are numerical diagnostics; the infinite-series identity follows from absolute convergence. The package's `wick_rotate(y,beta)` computes and compares both ensembles at user-supplied parameters. It does not perform the Mellin integral or derive a unique y↔beta point map.

## 6. Prediction accounting needs a dependency model

The current `prediction_table()` returns **23** rows. `honest_report()` labels its metadata grouping as **17** effective independent predictions. Excluding the explicitly definition-level generation-count row leaves **16 metadata families**. The Higgs quartic row 17 is exactly the squared/scaled Higgs mass relation in row 8; merging it into the Higgs family leaves **15**.

These are bookkeeping counts, not proven statistical independence. The rational torque comparison and the numerical-depth comparison are mathematical targets rather than separate measured particle observables. Several other rows share measured masses, and some predictions are conditional on empirical inputs: for example, the Koide tau calculation supplies electron/muon masses, and the cosmological row uses `PDG.alpha_inv` in its prediction.

Thus Mistral's 13–15 estimate is not a demonstrated count of independent empirical tests. Neither that estimate nor the package's label supports a joint success probability. A defensible analysis needs the dependency graph, experimental covariance and scheme conventions, an error model, and a documented formula-selection/search process. Theory-tolerance bands are not automatically probability distributions. No joint p-value or look-elsewhere correction was invented in this audit.

## 7. A systole can normalize, but cannot supply missing algebra

The trace calculation is reproducible exactly: trace 3 is excluded modulo 143, while a determinant-one trace-4 element exists, giving 2 arccosh(2)=2.633915793849633 in curvature-radius units. The metric here is the cusped quotient Gamma0(143)\\H, whose area is 56π. It should not silently be identified with the smooth compact genus-13 uniformizing metric, whose area at curvature −1 is 48π.

A declared systole normalization is basis-independent within its specified geometry. It does not choose a unique map from operators to geodesic length, change the inertia of an anticommutator form under positive rescaling, or remove a non-scalar anticommutator residual.

## Recommended next step

Pursue an **arithmetic-selected Clifford action on a specified active module**, using D4 and TH2 as constraints. First fix the real structure, projector transport, and candidate construction. Then compute the full anticommutators, their scalar Gram matrix, its rank and inertia, and every non-scalar residual. Merely finding standard gamma matrices somewhere inside a copy of so(8) would not show that arithmetic selects them.

`CLIFFORD_GATE_NEXT.md` gives the bounded proposal and elementary module obstructions, including the exact impossibility of an invertible anticommuting pair of complex 13×13 matrices. That statement does not prohibit real-linear operators on the 26-dimensional real homology or operators on an enlarged space.

All release-source files were left unchanged. The companion bundle preserves scripts, source hashes, exact controls, numerical results, and the correction to our previous Higgs wording. No full 4¹³ Pfaffian sum or theta-gradient experiment was run.
