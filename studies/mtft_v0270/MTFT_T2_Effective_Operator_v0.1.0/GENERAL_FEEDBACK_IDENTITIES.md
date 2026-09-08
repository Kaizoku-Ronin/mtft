# General Hermitian feedback identities

This note treats a finite Hermitian operator, without assuming an involution.
It supplies exact algebra for the T2 experiment; it does not certify the
numerically reconstructed MTFT subspaces. The companion script uses synthetic
rational matrices and reads no MTFT outputs.

## Schur reduction and its domain

For an orthogonal direct sum \(\mathcal H=P\mathcal H\oplus Q\mathcal H\), write

\[
H=\begin{pmatrix}A&B\\B^\dagger&D\end{pmatrix}=H^\dagger,
\qquad B:Q\mathcal H\longrightarrow P\mathcal H.
\]

For \(z\notin\operatorname{spec}(D)\), define

\[
\Sigma(z)=B(zI_Q-D)^{-1}B^\dagger,\qquad
H_{\rm eff}(z)=A+\Sigma(z),\qquad S(z)=zI_P-H_{\rm eff}(z).
\]

If also \(z\notin\operatorname{spec}(H)\), ordinary block elimination gives

\[
\boxed{G_{PP}(z):=P(zI-H)^{-1}P\big|_{P\mathcal H}=S(z)^{-1}.}
\]

Complement invertibility is an explicit hypothesis of the standard
Feshbach–Schur construction; see
[Dusson, Sigal and Stamm, Theorem 1.2 and equations (1.15)–(1.16)](https://arxiv.org/abs/2105.02058).
One must not silently replace an undefined inverse by a pseudoinverse at a
compression eigenvalue.

For \(\operatorname{Im}z=\eta>0\), all these inverses exist, and

\[
\operatorname{Im}\Sigma(z)
=-\eta B[(\operatorname{Re}z\,I-D)^2+\eta^2I]^{-1}B^\dagger\preceq0.
\]

Thus \(\operatorname{Im}S(z)\succeq\eta I_P\). This ensures invertibility but
does not ensure small floating-point forward error near a compression pole.

There is no general analogue of \(BB^\dagger=I-A^2\). The compression A alone
does not determine the feedback: keeping A and B fixed while changing D
already changes \(\Sigma\). At large spectral parameter,

\[
\Sigma(z)=\sum_{j\ge0}\frac{BD^jB^\dagger}{z^{j+1}},\qquad |z|>\|D\|.
\]

## Spectral residues and projected zeros

Let \(D=\sum_d d\Pi_d\), where the sum is over its distinct eigenvalues. Then

\[
\boxed{\Sigma(z)=\sum_d\frac{C_d}{z-d},\qquad
C_d=B\Pi_dB^\dagger\succeq0.}
\]

For a simple normalized eigenvector v, the residue is
\(C_d=(Bv)(Bv)^\dagger\). With degeneracy, one uses the entire eigenspace
projector; the result is independent of the eigenbasis. Its rank is
\(\operatorname{rank}(B\Pi_d)\). Zero residue means a removable singularity
of this self-energy term, although the original complement inverse is still
undefined at that argument.

Suppose \(d\in\operatorname{spec}(D)\) but
\(d\notin\operatorname{spec}(H)\). The full inverse is then regular. The
upper-right block of \((zI-H)^{-1}(zI-H)=I\) gives

\[
G_{PP}(z)B=G_{PQ}(z)(zI-D).
\]

Therefore \(G_{PP}(d)Bv=0\) for every \(v\in\ker(D-dI)\). In fact,

\[
\boxed{\ker G_{PP}(d)=B\ker(D-dI).}
\]

For the reverse inclusion, if \(G_{PP}(d)u=0\), write
\((dI-H)^{-1}(u,0)=(0,v)\); multiplication by \(dI-H\) gives
\(u=-Bv\) and \((dI-D)v=0\). Moreover B is injective on this eigenspace:
otherwise \((0,v)\) would be a full eigenvector at d, contradicting the
hypothesis. Hence the nullity of \(G_{PP}(d)\) equals the multiplicity of d
in D, under this disjoint-spectrum hypothesis.

This statement concerns a **kernel direction of the projected response**, not
the vanishing of the entire projected matrix. When d also belongs to the full
spectrum, these regular-resolvent arguments do not apply; the full spectral
residues must be examined separately. Poles of the effective operator are not
additional full-system eigenvalues. The full projected resolvent has possible
poles only at eigenvalues of H.

## Exact amplitude memory and initial-state forcing

Choose the dimensionless mathematical evolution \(i\dot\psi=H\psi\),
\(\psi=(x,y)\). Variation of constants gives

\[
y(t)=e^{-itD}y(0)-i\int_0^t e^{-i(t-s)D}B^\dagger x(s)\,ds,
\]

\[
\boxed{i\dot x(t)=Ax(t)+Be^{-itD}y(0)
-i\int_0^t K(t-s)x(s)\,ds,\qquad
K(t)=Be^{-itD}B^\dagger=\sum_d C_d e^{-itd}.}
\]

The forcing term disappears only when the initial complement vanishes.
Neither a physical clock nor an energy scale is supplied by these equations.
The projection is onto a direct-sum space of amplitudes, not a density-matrix
partial trace over a tensor-product environment.

Let \(U(t)=e^{-itH}\). Block multiplication of its exact composition law yields

\[
\boxed{U_{PP}(t)U_{PP}(s)-U_{PP}(t+s)=-U_{PQ}(t)U_{QP}(s).}
\]

Near \((t,s)=(0,0)\), the right-hand side is

\[
ts\,BB^\dagger+O\bigl(|ts|(|t|+|s|)\bigr).
\]

Thus if \(B\ne0\), the projected propagator cannot equal the exponential of
any time-independent matrix on the same active amplitudes at every time. It
fails the required composition law. Retaining auxiliary amplitudes or the
memory equation preserves the exact response.

For initially active data, unitarity gives the exact norm balance

\[
U_{PP}(t)^\dagger U_{PP}(t)+U_{QP}(t)^\dagger U_{QP}(t)=I_P.
\]

## Initially dark amplitudes can leak later

If \(B^\dagger x_0=0\), the direct first-order transfer vanishes. Expanding the
full exponential nevertheless gives

\[
y(t)=-\frac{t^2}{2}B^\dagger Ax_0+O(t^3),
\]

\[
\boxed{\|y(t)\|^2=
\frac{t^4}{4}\|B^\dagger Ax_0\|^2+O(t^5).}
\]

The quartic coefficient need not be strictly positive for every dark vector: it
vanishes if that vector also satisfies \(B^\dagger Ax_0=0\). Permanent absence
of leakage requires \(B^\dagger A^j x_0=0\) for all j, equivalently for
\(0\le j<\dim P\) by Cayley–Hamilton.

The general squared-norm remainder is **O(t^5), not O(t^6)**. If
\(v=B^\dagger Ax_0\) and
\(w=B^\dagger A^2x_0+DB^\dagger Ax_0\), its fifth-order coefficient is
\(\operatorname{Im}(v^\dagger w)/6\). Real H and real initial data make the
norm curve even and therefore remove this term. The symbolic controls include
a complex Hermitian example with
\(\|y(t)\|^2=t^4/4+t^5/6+O(t^6)\), showing that the weaker general remainder
is necessary.

## Minimum exact auxiliary space

Define the complement Krylov space

\[
\mathcal C=\operatorname{span}\{D^jB^\dagger x:
x\in P\mathcal H,\ 0\le j<\dim Q\}.
\]

Cayley–Hamilton makes this space invariant under D. Hermiticity also makes
\(\mathcal C^\perp\) invariant, and B vanishes on that orthogonal complement.
Consequently restricting D and B to \(\mathcal C\) preserves the exact
self-energy and memory kernel.

This auxiliary dimension is minimal for an exact finite linear realization
of the same rational self-energy. To see the lower bound, form the finite
block Hankel matrix with blocks \(BD^{i+j}B^\dagger\),
\(0\le i,j<\dim Q\). It is the Gram matrix of the block Krylov columns
\([B^\dagger,DB^\dagger,\ldots]\), so its rank is \(\dim\mathcal C\).
Any realization with k auxiliary coordinates factors this same Hankel matrix
through a k-dimensional space and must have \(k\ge\dim\mathcal C\).

In particular, **if B has full column rank five, all five complementary
amplitudes are necessary**: \(\operatorname{ran}B^\dagger=Q\mathcal H\)
already at Krylov order zero. The simpler high-frequency identity
\(\lim_{z\to\infty}z\Sigma(z)=BB^\dagger\), of rank five, gives the same
lower bound. This does not rule out approximation, or an energy-dependent
effective operator with no explicitly stored auxiliary states; it rules out
an exact realization with fewer than five linear auxiliary amplitudes on the
same active space.

## Independent symbolic controls

The primary synthetic example is

\[
H=\frac1{49}
\begin{pmatrix}-41&72&66\\72&109&6\\66&6&30\end{pmatrix},
\qquad \operatorname{spec}(H)=\{-2,1,3\}.
\]

Its exact rational orthogonal eigenbasis allows direct verification of the
Schur identity, unitary composition defect, both forms of the initial-state
memory equation, and norm exchange. The compression eigenvalue \(30/49\)
lies outside the full spectrum. The normalized dark input
\((1,-11)/\sqrt{122}\) satisfies
\(B^\dagger Ax_0=-180/(7\sqrt{122})\ne0\), so its transfer starts at quartic
order in squared norm. A second exact Hermitian control verifies the possible
fifth-order term. These controls validate the algebra; they are not arithmetic
data or a construction of physical interactions.
