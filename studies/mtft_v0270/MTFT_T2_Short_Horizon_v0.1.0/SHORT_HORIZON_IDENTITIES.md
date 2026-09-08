# A continuous bound for finite Hermitian active-space propagation

Let H be Hermitian on C^(p+q), let K be Hermitian on C^(p+r), and let E and F
be the isometric inclusions of the first p coordinates. Define

\[
U(t)=E^*e^{-itH}E,\qquad V(t)=F^*e^{-itK}F.
\]

These are direct-sum amplitude compressions. They are not density-matrix
partial traces. The results below require neither a tensor-product structure
nor an identification of H with a physical Hamiltonian.

## 1. A denominator that cannot vanish

Write the first p columns of exp(-itH) as the vertical stack (U,L).
Unitarity gives U*U+L*L=I_p. Since rank(L)<=q, the matrix U*U has at least
p-q eigenvalues equal to1 when p>q. Consequently

\[
\|U(t)\|_F\ge\sqrt{p-q}\quad(t\in\mathbb R).
\]

Here p=8,q=5, giving sqrt(3). This lower bound applies to the exact explicitly
stored block-coordinate matrix, whose coordinate inclusion is exactly
isometric. It does not assume the original floating Hodge frame is exactly
unitary.

## 2. Taylor moments and a unitary remainder

Define M_k=E*H^kE-F*K^kF. Taylor's integral formula is

\[
e^{-itH}=\sum_{k=0}^m\frac{(-itH)^k}{k!}
 +\frac{(-i)^{m+1}}{m!}\int_0^t(t-s)^m e^{-isH}H^{m+1}\,ds.
\]

For t>=0, unitary invariance and contraction under left projection give

\[
\|E^*R_m(H,t)E\|_F
\le\frac{t^{m+1}}{(m+1)!}\|H^{m+1}E\|_F.
\]

The same bound holds with |t| for negative t. Thus for all |t|<=T,

\[
\frac{\|U(t)-V(t)\|_F}{\|U(t)\|_F}
\le \frac{1}{\sqrt{p-q}}\left[
\sum_{k=0}^{m}\frac{T^k}{k!}\|M_k\|_F
+\frac{T^{m+1}}{(m+1)!}
 (\|H^{m+1}E\|_F+\|K^{m+1}F\|_F)\right].
\]

Every term is nonnegative, so evaluating the right side at T bounds the
whole interval. A second valid, usually looser, remainder replaces the two
column norms by sqrt(p)(rho_H^(m+1)+rho_K^(m+1)), where each rho is any
operator-norm upper bound. For Hermitian matrices the maximum absolute row
sum is one such bound. The experiment records both.

This is a relative Frobenius matrix-error guarantee. It is not a guarantee
of the same relative error for every initial state, observable, or probability.

## 3. Exact arithmetic for a numerical model

Every finite float64 real component is a dyadic rational. Each frozen matrix
is encoded as (R+iI)/2^e with integer R,I. Hermiticity is checked by exact
integer equality. Matrix powers use Gaussian integer multiplication, with
their powers of2 denominators tracked exactly. Frobenius norm squares are
therefore rational numbers.

For any nonnegative rational x=a/b, the code uses

\[
l=\lfloor\sqrt{\lfloor a10^{60}/b\rfloor}\rfloor/10^{30}.
\]

Then l^2<=x; either l^2=x or (l+10^-30)^2>x. These give exact lower and upper
enclosures. Norms in the numerator are rounded upward; sqrt(3) in the
denominator is rounded downward. All final comparisons with1/100 are exact
rational comparisons. A software replay independently recomputes the
certificate from the stored integer matrices. It requires only Python's
standard library.

The enclosures certify the finite stored matrices. They do not enclose the
error of MTFT's numerical period matrices, Hodge frame, or original geometry.
The independent reviewer also recomputes the decisive bound using rational
matrix arithmetic through a separate implementation.

## 4. Why the result cannot be repeatedly restarted for free

In general U(t+s) differs from U(t)U(s): the complement stores amplitude
which can return later. The short-horizon bound is for one propagation from
the active space with zero initial complement amplitude. Repeating a reduced
active-only step would reset or discard memory and needs a new error bound.
No extension of this result to arbitrary total time follows by stitching
short intervals together.

## 5. Interpretation of the study

For the stored T2 model H has size13, K has size12, and both have the same
8x8 active block. The reduced model retains four of five complement
directions. Among the preselected horizons, T=1/2 passes the1% bound; T=3/4
does not obtain a certificate from this bound even though its sampled error
is below1%. A failed upper bound does not prove the actual error exceeds the
budget. At T=1 the sampled endpoint error itself exceeds1%.

These identities are proved here by elementary finite-dimensional linear
algebra. The numerical inputs and the motivation for choosing this subspace
come from the frozen T2 Memory Approximation parent study.
