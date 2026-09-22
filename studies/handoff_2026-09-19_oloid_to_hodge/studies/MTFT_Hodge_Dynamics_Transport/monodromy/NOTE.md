# Explicit monodromy and transport growth

This calculation studies the elliptic period connection

\[
\frac{dY}{dc}=A(c)Y,\qquad
A(c)=\begin{pmatrix}
-\frac1{4c}&\frac1{4c(c+1)}\\
-\frac14&\frac1{4c}
\end{pmatrix}
\]

on the complex parameter domain \(\mathbb C\setminus\{0,-1\}\).
It does not prescribe a physical time or a logistic parameter evolution.

## Paths and convention

The basepoint is \(c=1\). All small circles have radius \(0.2\) and
counterclockwise orientation.

- The loop around zero follows the positive real segment from 1 to 0.2,
  traverses its circle, and retraces that segment.
- The loop around minus one uses the stem
  \(1\to1+0.75i\to-0.8+0.75i\to-0.8\), traverses its circle, and
  retraces the stem in reverse.

Let \(S\) and \(U\) denote the resulting matrices around zero and minus
one, respectively. The ODE convention is left multiplication: following
the zero loop and then the minus-one loop produces \(US\).

## Exact character data

Eliminating the second component gives the scalar Picard–Fuchs equation

\[
u''+\left(\frac1c+\frac1{c+1}\right)u'
+\frac{4c-1}{16c^2(c+1)}u=0.
\]

The local exponents are:

| Place | Exponents / asymptotics | Monodromy trace |
|---|---|---:|
| \(c=0\) | \(-1/4,1/4\) | 0 |
| \(c=-1\) | \(0,0\), with logarithm | 2 |
| infinity | \(u\sim c^{-1/2}\), repeated with possible logarithm | −2 |

The trace of \(A\) is zero, so every fundamental matrix initialized at
the identity has determinant one. For these standard generators, the
global loop relation and the infinity exponents give

\[
\operatorname{tr}S=0,\quad S^2=-I,\qquad
\operatorname{tr}U=2,\quad U=I+N,\quad N^2=0,
\qquad \operatorname{tr}(US)=-2.
\]

The nilpotent \(N\) is nonzero: the local connection residue at minus one
is \(\begin{psmallmatrix}0&-1/4\\0&0\end{psmallmatrix}\), and its equal
eigenvalues cause no positive-integer obstruction to the usual analytic
reduction to residue form. Thus local monodromy has a genuine logarithmic
term. Also \(US\neq-I\), since otherwise \(U=-S^{-1}=S\), contradicting
their distinct traces.

Consequently

\[
\operatorname{tr}(U^kS)
=\operatorname{tr}((I+kN)S)=-2k.
\]

For \(k=2\), the exact characteristic polynomial is
\(\lambda^2+4\lambda+1\). Its eigenvalues are

\[
\lambda_\pm=-2\pm\sqrt3,
\qquad\rho(U^2S)=2+\sqrt3.
\]

This proves the following distinct growth behaviors in any fixed basis:

| Repeated transport | Exact reason | Asymptotic norm growth |
|---|---|---|
| \(S^n\) | \(S^4=I\) | bounded |
| \(U^n\) | \(I+nN\) | linear |
| \((US)^n\) | \((-1)^n[I-n(US+I)]\) | linear |
| \((U^2S)^n\) | eigenvalue magnitude \(2+\sqrt3\) | exponential |

The exponential growth per repetition of the last word is

\[
\lim_{n\to\infty}\frac1n\log\|(U^2S)^n\|
=\log(2+\sqrt3)=1.3169578969248166\ldots.
\]

It is essential to specify the word: simple alternation \(US\) is
parabolic, not hyperbolic. The hyperbolic result does not require a
chaotic driver. It is the growth of a linear cocycle under a chosen,
periodically repeated parameter path. Dividing by elapsed time would
require an additional traversal-time convention.

## Numerical checks

`compute_monodromy.py` integrates the unmodified connection with SciPy
DOP853 in complex128 at relative tolerances \(10^{-8},10^{-10},10^{-13}\),
using absolute tolerance one tenth the relative tolerance. It does not
normalize matrices to determinant one. The finest run finds:

| Quantity | Numerical value / defect |
|---|---:|
| \(\operatorname{tr}S\) | \(-1.15\times10^{-14}-1.24\times10^{-14}i\) |
| \(\operatorname{tr}U-2\) | \(-2.07\times10^{-14}-3.33\times10^{-16}i\) |
| \(\operatorname{tr}(US)+2\) | \(-6.40\times10^{-14}-1.29\times10^{-14}i\) |
| \(\operatorname{tr}(U^2S)+4\) | \(-7.46\times10^{-14}-1.25\times10^{-14}i\) |
| \(\rho(U^2S)\) | 3.732050807568958 |
| maximum determinant defect among reported matrices | \(5.44\times10^{-15}\) |
| maximum entry-matrix norm change from previous tolerance | \(3.59\times10^{-11}\) |

These are numerical consistency checks, not interval certificates.
Near the defective unipotent matrices, a generic eigenvalue routine
splits the repeated eigenvalue by roughly \(10^{-7}\), despite matrix
errors near \(10^{-14}\). Such tiny apparent eigenvalue growth must not
be mistaken for genuine hyperbolicity. Exact trace identities and
Cayley–Hamilton determine the parabolic classification.

In two dimensions, determinant preservation is equivalent to
\(M^TJM=J\) for \(J=\begin{psmallmatrix}0&1\\-1&0\end{psmallmatrix}\).
The transpose here is not conjugate transpose. This complex symplectic
property does not mean that transport preserves the ordinary Euclidean
Hermitian norm in this algebraic frame.

## Nonconstant gauge check

Set \(\widetilde Y=G(c)Y\), with

\[
G(c)=\operatorname{diag}(e^{0.31ic},e^{-0.31ic}).
\]

Then

\[
\widetilde A=G'G^{-1}+GAG^{-1},\qquad
\widetilde M=G(1)MG(1)^{-1}
\]

for every based closed loop. Independently integrating the transformed
connection around minus one gives a conjugacy defect of
\(3.55\times10^{-14}\). Trace, eigenvalues, and repeated-loop exponential
rates are invariant under this change of frame. Arbitrarily rescaling
states along an unbounded nonclosed trajectory requires additional norm
control; that broader Lyapunov-invariance claim is not used here.

The displayed matrices are in a complex de Rham period frame. This
experiment does not identify a specific integral homology basis, nor
assert that their entries themselves are integer matrices. The exact
character calculation above does not need such a basis.

## Reproduction

```sh
python3 compute_monodromy.py
```

Outputs are `monodromy_results.json` (full matrices, diagnostics, paths,
and tolerances) and `monodromy_growth.csv` (finite repeated-word norm
values). Dependencies are NumPy and SciPy. The CSV powers illustrate
the exact growth classifications; very long floating-point powers of
parabolic matrices are unsuitable as independent asymptotic proofs.
