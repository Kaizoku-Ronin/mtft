# Matter-03 kinetic normalization and mass gate

The Matter-02 residue certificate establishes a holomorphic rank lower bound. It does not supply the physical scale or the mode metrics. This audit specifies the field-coordinate conversion needed for a mass calculation and checks it independently with complex, nonorthogonal examples.

## Explicit tensor convention

Use three columns of holomorphic four-dimensional fields, with a lower-index superpotential tensor:

\[
W=\kappa\sum_{a,b,c}Y_{abc}L_a E_b H_c,
\qquad
K_{\rm quad}=L^\dagger G_L L+E^\dagger G_E E+H^\dagger G_H H.
\]

Here the positive Hermitian matrices are field-coordinate kinetic matrices. The coefficient \(\kappa\) includes any overall parent superpotential normalization that was omitted when the residue pairing was written only up to proportionality. It also needs to be derived or supplied. It is set to one in the numerical audit as an explicit INPUT.

Write each raw field as \(L=C_L L_{\rm can}\), and similarly for E and H. The conditions are

\[
C_X^\dagger G_X C_X=I.
\]

If \(G_X=R_XR_X^\dagger\) is a Cholesky factorization, one valid choice is \(C_X=(R_X^\dagger)^{-1}\). Substitution in the holomorphic superpotential gives

\[
Y^{\rm can}_{ijk}
=\sum_{a,b,c}(C_L)_{ai}(C_E)_{bj}(C_H)_{ck}Y_{abc}.
\]

There is no complex conjugation on a holomorphic field slot. In particular, for a canonically normalized Higgs direction \(n\), \(n^\dagger n=1\), and neutral expectation value \(\langle H_{\rm can}^0\rangle=(v_H/\sqrt2)n\),

\[
M_{\rm can}=\frac{\kappa v_H}{\sqrt2}\,
C_L^T\left[\sum_cY_c(C_Hn)_c\right]C_E.
\]

The nonnegative singular values of this matrix are the tree-level Dirac masses in this declared effective-field-theory convention. They are not the eigenvalues of the generally complex, nonsymmetric Yukawa matrix. For one raw Higgs coordinate with kinetic coefficient \(G_H>0\), the Higgs factor reduces to \(G_H^{-1/2}\). For several Higgs coordinates a scalar rescaling is insufficient.

The scale \(v_H\) is the expectation value of the particular down-Higgs sector. It is not automatically the total electroweak scale if other Higgs sectors also acquire expectation values. The audit chooses \(v_H=1\) in arbitrary units and does not infer an electroweak value.

## Covariance under a raw basis change

Let \(L=B_L L'\), \(E=B_E E'\), \(H=B_H H'\), with every B invertible. Then

\[
G'_X=B_X^\dagger G_X B_X,
\qquad
Y'_{ijk}=\sum_{a,b,c}(B_L)_{ai}(B_E)_{bj}(B_H)_{ck}Y_{abc}.
\]

For canonical frames computed separately in the two raw bases, define

\[
U_X=C_X^{-1}B_XC'_X.
\]

Every U is unitary. The same physical Higgs expectation value is represented by \(n'=U_H^\dagger n\), and the mass matrices obey

\[
M'_{\rm can}=U_L^T M_{\rm can}U_E.
\]

Thus singular values and rank are invariant. Keeping the same numerical n after changing the Higgs basis generally chooses a different vacuum direction; that is not a basis-invariance test.

The historical `surface.hym.normalise_yukawa` routine works with a different section-basis convention and a multiplication tensor whose Higgs index is a target/dual index. Its use of an inverse Cholesky matrix on section indices and a Cholesky or dual conversion on the target index is the archived CC-26 convention. The procedure above is specifically for lower-index field-coordinate W. This audit does not change or reinterpret the historical routine.

## Why rank does not fix a mass hierarchy

For the Matter-02 selected minor,

\[
Y_0=\begin{pmatrix}1&0&0\\0&0&1\\0&1&0\end{pmatrix},
\]

identity kinetic matrices give three equal singular values. The equally positive supplied matrices

\[
G_L=\operatorname{diag}(1,4,25),\quad
G_E=\operatorname{diag}(1,9,16),\quad G_H=1
\]

give singular values \(1,1/8,1/15\) for the normalized coupling. Both examples have rank three. These are examples of dependence on missing data, not predicted MTFT lepton ratios. A 3 × 3 minor selected from the 18 × 44 tensor is also not an identified three-family low-energy sector.

If light subspaces were selected independently, columns \(S_L\) and \(S_E\) would embed them in the complete mode spaces. Their metrics would be \(S_L^\dagger G_LS_L\) and \(S_E^\dagger G_ES_E\); the raw projected coupling would be \(S_L^T Y(H)S_E\). Those embeddings need an actual vacuum and a consistent chirality dictionary. Selecting the residue-witness columns solely because they make a nonzero determinant does not solve the light-family problem.

## Metric inventory

| Data | Present status | Relation to Matter-02 |
|---|---|---|
| Full cusp basis and ordinary-cusp jets | Frozen exact input; exact finite certificate | Defines the chosen algebraic geometry and residue witness |
| 18 × 44 × 10 holomorphic Yukawa tensor | Only rank-three slice computed | Full tensor remains needed for a complete mass operator |
| 18 × 18 lepton kinetic matrix | Not constructed | Requires harmonic representatives on the blown-up surface |
| 44 × 44 antilepton kinetic matrix | Not constructed | Requires its own harmonic representatives and line metric |
| 10 × 10 down-Higgs kinetic matrix | Not constructed | A curve differential Gram alone does not suffice |
| Overall parent superpotential coefficient κ | Pairing previously stated up to proportionality | Needed for absolute coupling magnitudes |
| Smooth compact curve HYM/Green-function machinery | Existing `surface.hym`; DIAGNOSTIC | Useful computational machinery for older curve models |
| Petersson differential Gram and period comparison | Existing `surface.petersson`; DIAGNOSTIC | Curve one-form norms, not B3 charged (0,1)/(2,0) mode metrics |
| Frozen M2 tensors | Existing older model; DIAGNOSTIC | Different spaces and metric; cannot be imported as B3 results |
| Higgs vacuum and light-family embeddings | OPEN | Determines which mass operator should be normalized |

A physical construction needs a smooth Kähler metric on the blow-up, a compatible Hermitian line-bundle metric solving the selected HYM problem, global harmonic mode representatives, the actual parent kinetic prefactors, and normalization of the Higgs expectation value. Algebraic cohomology dimensions and local principal parts do not provide these metric integrals. Schematically the charged (0,1) modes use inner products \(\int_S\langle a_i,a_j\rangle_{g,h}\,d\mathrm{vol}_g\), while the Higgs modes require the corresponding (2,0)-form inner product. The trace normalization and couplings from the parent action also enter.

## Positive D-terms and the vacuum scale

In a canonical global supersymmetric gauge theory, eliminating auxiliary D fields gives a sum of nonnegative squares. For one weak doublet with hypercharge magnitude 1/2 and no FI term, the quartic piece is

\[
V_D=\frac{g_2^2}{8}\sum_A(H^\dagger\sigma^A H)^2
+\frac{g_Y^2}{8}(H^\dagger H)^2
=\frac{g_2^2+g_Y^2}{8}(H^\dagger H)^2.
\]

The last equality follows from the exact Pauli identity checked symbolically. This term supplies a positive quartic, not a negative quadratic mass. For neutral up/down components u and d, the electroweak restriction is

\[
V_D=\frac{g_2^2+g_Y^2}{8}(|u|^2-|d|^2)^2.
\]

The ray \(|u|=|d|\) is D-flat at all scales. A positive quartic or a flat ray by itself does not select a finite nonzero Higgs expectation value. These electroweak examples omit the additional Abelian generators and parent constraints of B3; they are illustrative INPUT effective theories, not a solved B3 vacuum.

For a supplied ordinary single-doublet potential \(V=m^2H^\dagger H+\lambda(H^\dagger H)^2\), \(\lambda>0\), a broken minimum requires \(m^2<0\), with \(v^2=-m^2/\lambda\). Nothing in the Matter-02 residue rank certificate derives this negative mass coefficient, a soft sector, a μ term, or quantum effective potential. In the existing unbroken split construction the chosen cohomology modes are massless fluctuations; a later off-diagonal Higgs expectation value can change the zero-mode complex and invalidate reuse of the split dimensions.

## Verification

`check_kinetics.py` passes 117 checks. Eleven are exact symbolic identities/positivity statements; the remaining checks are finite numerical examples and input rejection. The exact witness determinant after supplied rational-complex normalization is −4/81. Twenty-four independent complex changes of every field basis preserve the mass singular values when the same physical Higgs direction is transported. The complex example detects the incorrect use of a conjugate transpose on the first holomorphic slot.

All numerical residuals in `kinetic_checks.json` are diagnostics for these supplied finite examples. They are not bounds on the missing B3 mode metrics. The test does not run a vacuum search or predict masses.

## Primary reference

Stephen P. Martin, *A Supersymmetry Primer*, arXiv:hep-ph/9709356v7 (27 January 2016), sections 3.2–3.5, 4.10 and 8.1: https://arxiv.org/pdf/hep-ph/9709356 . Standard superpotential, Kähler metric, auxiliary-field, and Higgs-potential framework. The field-basis covariance proof and explicit examples above are derived for this audit.
