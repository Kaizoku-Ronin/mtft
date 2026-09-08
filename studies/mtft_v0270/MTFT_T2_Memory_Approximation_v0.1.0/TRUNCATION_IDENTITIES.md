# Exact identities for a structure-preserving memory approximation

This note concerns finite-dimensional Hermitian matrices and direct-sum amplitude coordinates. Its statements do not depend on MTFT measurements. The companion script uses independent synthetic rational examples. The measured T2 study is a separate numerical application.

## 1. The retained auxiliary model

Let

\[
H=\begin{pmatrix}A&B\\B^*&D\end{pmatrix}=H^*,
\qquad A\in\mathbb C^{n\times n},\quad D\in\mathbb C^{m\times m}.
\]

Choose an orthonormal auxiliary basis \(Q_r\in\mathbb C^{m\times r}\), and complete it by \(Q_s\) so that \([Q_r\ Q_s]\) is unitary. Define

\[
B_r=BQ_r,\quad D_r=Q_r^*DQ_r,
\qquad H_r=\begin{pmatrix}A&B_r\\B_r^*&D_r\end{pmatrix}.
\]

This Galerkin compression is Hermitian for every chosen subspace. No fit of a general non-Hermitian response matrix is required. Its feedback is

\[
\Sigma_r(z)=B_r(zI-D_r)^{-1}B_r^*,\qquad
G_r(z)=\left[zI-A-\Sigma_r(z)\right]^{-1}.
\]

The feedback formula requires \(z\notin\operatorname{spec}(D_r)\); the response additionally requires \(z\notin\operatorname{spec}(H_r)\). Both conditions hold whenever \(\operatorname{Im}z>0\).

Write the spectral decomposition \(D_r=\sum_d d\Pi_d\). Then

\[
\Sigma_r(z)=\sum_d\frac{B_r\Pi_dB_r^*}{z-d}.
\]

Every residue is positive semidefinite; a residue may vanish. For \(z=x+i\eta\), \(\eta>0\),

\[
-\operatorname{Im}\Sigma_r(z)
=\eta B_r[(xI-D_r)^2+\eta^2I]^{-1}B_r^*\succeq0.
\]

Consequently the Schur denominator has imaginary part at least \(\eta I\). The full reduced resolvent obeys

\[
-\operatorname{Im}(zI-H_r)^{-1}
=\eta(zI-H_r)^{-*}(zI-H_r)^{-1}\succ0.
\]

These identities preserve analyticity and the appropriate resolvent sign in the upper half-plane. They do not imply that a small auxiliary model is accurate. Its real eigenvalues can differ substantially from the full model's eigenvalues, though Galerkin compression keeps them within the full spectral interval.

For real dimensionless \(t\), \(e^{-itH_r}\) is unitary. If the reduced initial vector is \((x_0,0)\), its active and retained auxiliary amplitudes satisfy

\[
\|x_r(t)\|^2+\|y_r(t)\|^2=\|x_0\|^2.
\]

The active block alone need not be unitary or obey a composition law. Norm conservation of the reduced full model does not mean its active response agrees with the original model, nor does it imply that active norm loss is always overestimated or always underestimated.

## 2. A full-space embedding and an auditable error scale

In active, retained, discarded coordinates, put

\[
B_s=BQ_s,\quad M=Q_r^*DQ_s,\quad D_s=Q_s^*DQ_s.
\]

The unitarily transformed full matrix and a decoupled embedding of the approximation are

\[
\widetilde H=
\begin{pmatrix}A&B_r&B_s\\B_r^*&D_r&M\\B_s^*&M^*&D_s\end{pmatrix},
\qquad H_0=H_r\oplus D_s.
\]

Thus

\[
E=\widetilde H-H_0=
\begin{pmatrix}0&C\\C^*&0\end{pmatrix},
\qquad C=\begin{pmatrix}B_s\\M\end{pmatrix}.
\]

The top zero block has dimension \(n+r\). Since

\[
E^2=\operatorname{diag}(CC^*,C^*C),
\]

the spectral operator norms obey the exact identity

\[
\boxed{\|E\|_2=\|C\|_2},\qquad
\|C\|_2^2=\lambda_{\max}(B_s^*B_s+M^*M).
\]

This expression separates two effects: direct coupling from active coordinates to discarded coordinates, and mixing between retained and discarded auxiliary coordinates. Rotating either auxiliary basis internally leaves the norm unchanged.

Let \(P_a\) extract the original active coordinates, and write

\[
G(z)=P_a(zI-\widetilde H)^{-1}P_a^*,\qquad
G_r(z)=P_a(zI-H_0)^{-1}P_a^*.
\]

The second equality uses the fact that the discarded block of \(H_0\) is decoupled. For \(\operatorname{Im}z=\eta>0\), both full resolvents have operator norm at most \(1/\eta\). Applying the resolvent identity gives

\[
\boxed{\|G(z)-G_r(z)\|_2\le\frac{\|C\|_2}{\eta^2}}.
\]

The corresponding Duhamel identity, with unitary factors on both sides of \(E\), gives

\[
\boxed{\|P_a e^{-it\widetilde H}P_a^*
-P_a e^{-itH_0}P_a^*\|_2
\le\min(2,|t|\,\|C\|_2)}.
\]

The bound by 2 follows because both projected propagators are contractions. These are **absolute spectral operator-norm bounds**. They are not relative Frobenius errors and do not by themselves certify a 1% relative response budget. A norm conversion and a justified denominator bound would be needed for that conclusion.

There is a useful optional refinement. Set \(R=(zI-\widetilde H)^{-1}\), \(R_0=(zI-H_0)^{-1}\). The exact expansion

\[
R-R_0=R_0ER_0+R_0ERER_0
\]

has zero retained-retained first term: \(R_0\) is block diagonal and \(E\) is block off-diagonal. Therefore

\[
\|G-G_r\|_2\le
\min\left(\frac{\|C\|_2}{\eta^2},
\frac{\|C\|_2^2}{\eta^3}\right).
\]

This mathematical refinement changes no preregistered selection rule or evaluation gate.

## 3. Why a small discarded singular value is insufficient

A singular-vector truncation of \(B\) controls \(\|BQ_s\|\). It need not control \(M=Q_r^*DQ_s\). Thus a weakly or even zero directly coupled auxiliary direction can still receive an amplitude through a retained auxiliary direction and later feed back into the active sector.

The projected propagator expansions give an exact local comparison:

\[
P_a e^{-itH}P_a^*
=I-itA-\frac{t^2}{2}(A^2+BB^*)+O(t^3),
\]

\[
P_a e^{-itH_r}P_a^*
=I-itA-\frac{t^2}{2}(A^2+B_rB_r^*)+O(t^3).
\]

Their difference begins as

\[
\boxed{-\frac{t^2}{2}B(I-Q_rQ_r^*)B^*+O(t^3)}.
\]

If \(B_s=0\), the first three active matrix moments agree. The fourth moments instead satisfy

\[
(\widetilde H^4)_{aa}-(H_r^4)_{aa}=B_rMM^*B_r^*,
\]

so the leading possible difference is

\[
\frac{t^4}{24}B_rMM^*B_r^*+O(t^5).
\]

The coefficient can vanish in special cases, so fourth order is a possibility, not a universal strictly nonzero term.

For a concrete exact counterexample, take one active and two auxiliary coordinates:

\[
H=\begin{pmatrix}0&1&0\\1&0&1\\0&1&0\end{pmatrix},
\qquad H_r=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]

The discarded column of \(B=(1,0)\) is exactly zero. Nevertheless,

\[
U_{aa}(t)=\frac{1+\cos(\sqrt2t)}2,
\qquad U_{r,aa}(t)=\cos t,
\qquad U_{aa}-U_{r,aa}=\frac{t^4}{24}+O(t^6).
\]

The discarded amplitude is \((\cos(\sqrt2t)-1)/2\), whose squared norm starts at \(t^4/4\). The active resolvents differ by

\[
G(z)-G_r(z)=\frac{1}{z(z^2-2)(z^2-1)}.
\]

The original memory kernel is \(K(t)=\cos t\), while the retained memory is \(K_r(t)=1\). Matching the initial kernel therefore does not ensure agreement at later times. All these equalities are checked symbolically.

## 4. Exact dimension versus approximation

For the exact feedback,

\[
\Sigma(z)=z^{-1}BB^*+O(z^{-2}).
\]

An \(r\)-coordinate realization with constant couplings has a leading coefficient of rank at most \(r\). Therefore, **if** \(\operatorname{rank}B=5\), fewer than five coordinates cannot reproduce the exact feedback as a rational function of \(z\). This is a statement about exact realization, independent of any finite-grid approximation experiment.

It does not rule out a useful smaller approximation on a specified spectral domain and finite time horizon. Such a claim requires a declared error norm, parameter range, budget, selection procedure, and evaluation points that were not used for selection. Passing a finite grid is evidence on that grid; it is not a uniform continuum guarantee. The bounds above are available as uniform but potentially conservative controls.

## 5. Validation and interpretation

`truncation_identities.py` checks a rational 2-active/2-auxiliary example with a nontrivial orthonormal auxiliary rotation and nonzero retained/discarded mixing. It checks the block embedding, norm identity, resolvent identity, Schur formula, positive residue, upper-half-plane sign, and short-time moments. A separate three-coordinate example checks exact finite unitary evolution and the indirect-coupling counterexample. Results and the script hash are recorded in `truncation_identity_results.json`.

These synthetic controls verify identities; they do not validate measured T2 coefficients or establish its approximation error. The general inequalities are proved above, rather than inferred from numerical examples.

Throughout, sectors are summands of an amplitude space. There is no tensor-product environment, density-matrix partial trace, or asserted physical Hamiltonian in this construction. Real parameters called time here are dimensionless parameters of a chosen Hermitian matrix evolution.
