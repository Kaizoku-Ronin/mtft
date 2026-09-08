# Exact identities for the W143 feedback experiment

These are finite-dimensional algebraic statements, conditional on the stated
operator identities. They do not certify the numerically reconstructed MTFT
projectors. The companion `feedback_identities.py` verifies 27 identities exactly
for a generic synthetic two-dimensional coupled channel using Sympy.

## Setup and exact reduction

Let the orthogonal splitting be \(\mathcal H=P\mathcal H\oplus Q\mathcal H\),
and write a self-adjoint involution as

\[
W=\begin{pmatrix}A&B\\B^\dagger&D\end{pmatrix},\qquad W^\dagger=W,\quad W^2=I.
\]

Here \(B:Q\mathcal H\to P\mathcal H\). Multiplication gives

\[
A^2+BB^\dagger=I_P,\quad D^2+B^\dagger B=I_Q,\quad AB+BD=0.
\]

For \(z\notin\operatorname{spec}(D)\), define

\[
\Sigma(z)=B(zI_Q-D)^{-1}B^\dagger,\qquad H_{\rm eff}(z)=A+\Sigma(z).
\]

If also \(z\notin\operatorname{spec}(W)\), block elimination yields

\[
G_{PP}(z):=P(zI-W)^{-1}P\big|_{P\mathcal H}
=[zI_P-H_{\rm eff}(z)]^{-1}.
\]

The complement-invertibility condition is part of the standard Feshbach–Schur
construction; see Theorem 1.2 and equations (1.15)–(1.16) of
[Dusson, Sigal and Stamm](https://arxiv.org/abs/2105.02058).

The involution supplies a separate closed form, valid for every \(z\ne\pm1\):

\[
(zI-W)^{-1}=\frac{zI+W}{z^2-1},\qquad
\boxed{G_{PP}(z)=\frac{zI_P+A}{z^2-1}.}
\]

This last formula remains valid at interior eigenvalues of \(D\), where literal
Schur evaluation is undefined. No inverse of a singular complement is taken.

On the common domain of the expressions,

\[
(zI_P+A)\Sigma(z)=I_P-A^2,
\]

and, when \(zI_P+A\) is invertible,

\[
\Sigma(z)=(I_P-A^2)(zI_P+A)^{-1},\qquad
H_{\rm eff}(z)=(zA+I_P)(zI_P+A)^{-1}.
\]

These rational identities justify continuation, rather than evaluating the
original inverse at a forbidden argument.

## Canonical coupled channels and their poles

Suppose \(Dq=dq\), \(\|q\|=1\), and \(|d|<1\). Set
\(b=\sqrt{1-d^2}>0\) and \(p=Bq/b\). The block identities imply
\(\|p\|=1\), \(Ap=-dp\), and \(B^\dagger p=bq\). Thus each such mode gives

\[
W_j=\begin{pmatrix}a_j&b_j\\b_j&-a_j\end{pmatrix},\qquad
a_j=-d_j,\quad a_j^2+b_j^2=1.
\]

There are exactly \(\operatorname{rank}B\) coupled channels, counting
multiplicity. All remaining modes of \(A\) and \(D\) have eigenvalues in
\(\{\pm1\}\) and no inter-sector coupling.

For one coupled channel,

\[
H_{{\rm eff},j}(z)=a_j+\frac{b_j^2}{z+a_j}
=\frac{1+a_jz}{z+a_j},\qquad
G_{PP,j}(z)=\frac{z+a_j}{z^2-1}.
\]

At \(z=-a_j=d_j\), the effective operator has a simple pole with positive
residue \(b_j^2\), while the corresponding projected full resolvent has a
zero. The full two-dimensional block has eigenvalues \(+1,-1\), regardless
of \(a_j\). The compression poles therefore do not introduce extra full-system
eigenvalues. Coincident \(d_j\) combine residues into a higher-rank residue at
one pole. Complement eigenvalues \(\pm1\) have zero coupling residue; literal
complement inversion is still undefined there, and \(\pm1\) may remain actual
poles of the full resolvent.

If \(\dim P=8\), \(\dim Q=5\), \(\operatorname{rank}B=2\), and the positive
eigenspace of \(W\) has dimension two, the two coupled blocks exhaust all
positive directions. Consequently the six uncoupled active directions and
three uncoupled complementary directions all carry \(W=-I\). The minimal
coupled model is **four complex dimensions: two active and two complementary**.
It consists of two two-dimensional channels, not a two-dimensional full model.

For completeness, if \(W=-I+2R_+\), let \(\lambda_j\) be the two eigenvalues of
\(R_+PR_+\) restricted to \(\operatorname{ran}R_+\). They give
\(a_j=2\lambda_j-1\) and \(b_j=2\sqrt{\lambda_j(1-\lambda_j)}\), when both
angles are nontrivial. These are properties of the chosen pair of subspaces.

## Exact feedback and norm exchange

Introduce the dimensionless mathematical evolution \(i\dot\psi=W\psi\),
with \(\psi=(x,y)\). This is a toy generator; physical time or energy has not
been supplied. Eliminating \(y\) gives

\[
y(t)=e^{-itD}y(0)-i\int_0^t e^{-i(t-s)D}B^\dagger x(s)\,ds,
\]

\[
\boxed{
i\dot x(t)=Ax(t)+Be^{-itD}y(0)
-i\int_0^t B e^{-i(t-s)D}B^\dagger x(s)\,ds.}
\]

The initial-complement term only vanishes when \(y(0)=0\). The exact memory
kernel is

\[
K(t)=Be^{-itD}B^\dagger=e^{itA}(I_P-A^2),
\]

using \(AB=-BD\). Since \(e^{-itW}=\cos t\,I-i\sin t\,W\), initial data
\((x_0,0)\) evolve as

\[
x(t)=(\cos t\,I_P-i\sin t\,A)x_0,\qquad
y(t)=-i\sin t\,B^\dagger x_0,
\]

and hence

\[
\boxed{\|x(t)\|^2=\|x_0\|^2-\sin^2t\,\|B^\dagger x_0\|^2},\qquad
\|y(t)\|^2=\sin^2t\,\|B^\dagger x_0\|^2.
\]

The norm lost from the active sector is exactly gained by its complement.
For a normalized single-channel input the loss is \(b_j^2\sin^2t\).
The reduced propagator \(U_{PP}(t)\) has the exact composition defect

\[
U_{PP}(t)U_{PP}(s)-U_{PP}(t+s)
=\sin t\sin s\,(I_P-A^2).
\]

Evolution with \(e^{-itA}\) alone preserves the active norm and misses this
feedback. The revival at \(t=\pi\) follows automatically from the involution;
it is not a discovered physical resonance.

## Numerical interpretation

The frozen protocol supplies the actual spectral grid. On its upper-half-plane
points, self-adjointness guarantees all required inverses exist. Use linear
solves, record residuals and conditioning, and compare the direct full solve,
Schur solve, and independent involution formula. For \(\operatorname{Im}z=\eta>0\),

\[
\operatorname{Im}\Sigma(z)
=-\eta B[(\operatorname{Re}z\,I-D)^2+\eta^2 I]^{-1}B^\dagger\preceq0,
\]

so the Schur matrix has imaginary part at least \(\eta I\). Approaches to
\(d_j+i\eta\) as \(\eta\downarrow0\) are separately registered conditioning
diagnostics. Their large self-energy and small projected channel response are
predicted cancellation behavior. Evaluation at \(z=d_j\) should use the full
solve or its closed form, never a pseudoinverse silently substituted into the
Schur formula.

All numerical claims about the actual MTFT matrices inherit the accuracy of
the Hodge metric, projectors, and frame transport. These exact identities do
not convert floating-point residuals into arithmetic certificates.
