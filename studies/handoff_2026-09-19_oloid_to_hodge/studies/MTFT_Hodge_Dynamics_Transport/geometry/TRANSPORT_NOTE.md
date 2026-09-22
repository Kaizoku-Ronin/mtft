# Exact transport and splitting through genus three

This calculation works over `Q(c)` and concerns the explicitly chosen quadratic
family `f_c(x)=x^2+c`. It does not derive that choice from MTFT, construct a
Feigenbaum renormalization operator, or assert a physical gauge theory.

## 1. Three curves and two quotient maps

Write `P_m=f_c^{\circ m}` and take smooth projective compactifications of

\[
C_2:\quad w^2=P_2(u)=(u^2+c)^2+c,
\]
\[
C_3:\quad y^2=P_3(x)=P_2(x^2+c),
\]
\[
D_2:\quad v^2=Q(u,c)=(u-c)P_2(u).
\]

Their genera are respectively 1, 3, and 2. Throughout this note,

\[
c(c+1)H(c)\ne0,\qquad H(c)=c^3+2c^2+c+1.
\]

The quotient maps from `C_3` are

\[
\pi_+(x,y)=(x^2+c,y),\qquad
\pi_-(x,y)=(x^2+c,xy).
\]

For `sigma(x,y)=(-x,y)` and the hyperelliptic involution
`iota(x,y)=(x,-y)`, the first quotient is by `sigma`, and the second by
`sigma iota`. These are different involutions. Pullback gives the rational
de Rham/Hodge splitting

\[
H^1(C_3)\simeq\pi_+^*H^1(C_2)\oplus\pi_-^*H^1(D_2),
\]

and consequently `J(C_3)` is isogenous to `J(C_2) x J(D_2)`. This is an
isogeny-level statement; it is not an assertion of an integral product of
polarized Jacobians.

## 2. A residue-free basis matters

On the odd-degree genus-two curve `D_2`, the four classes

\[
\eta_j=u^j\,du/v,\qquad 0\le j\le3,
\]

form a de Rham basis. The first two are holomorphic; the last two are
differentials of the second kind. At the unique point at infinity, using
`u=t^{-2}`, their Laurent powers are even, so there is no residue term.

For the even-degree curve `C_3`, using all monomials naively would introduce
third-kind differentials. In particular, `x^3 dx/y` has nonzero residues at
the two points at infinity, and `x^5 dx/y` alone also generally has residues.
The following six pulled-back forms are instead a valid residue-free basis:

\[
\begin{aligned}
b_0&=2x\,dx/y,\\
b_1&=2x(x^2+c)^2\,dx/y,\\
b_{j+2}&=2(x^2+c)^j\,dx/y\quad (0\le j\le3).
\end{aligned}
\]

For example `b_1=2(x^5+2cx^3+c^2x) dx/y`; its residue contributions cancel.
The inherited two-dimensional sector is spanned by `b_0,b_1`. The new
four-dimensional sector is spanned by `b_2,...,b_5`. Their holomorphic parts
have dimensions one and two, respectively.

## 3. Explicit genus-two Gauss–Manin matrix

For a locally transported closed cycle, let `I_j=\oint\eta_j` and take
`I=(I_0,I_1,I_2,I_3)^T`. Then `dI/dc=A_D(c) I`.
The full matrix is specified by

\[
A_D(c)=\frac{B(c)}{8cH(c)},
\]

with the following four rows of `B`:

\[
B_{0,*}=\left(
-6c^3-13c^2-5c-3,\quad
\frac{c(5-2c-5c^2)}{c+1},\quad
\frac{2c^3+c^2-2c+1}{c+1},\quad
\frac{-3c^2+3c+9}{c+1}
\right),
\]
\[
B_{1,*}=\left(
c(-c^2+c+3),\quad
-4c^3+3c^2+c-1,\quad
-c(3c+1),\quad3(3c+1)
\right),
\]
\[
B_{2,*}=\left(
-c(2c^3+c^2-2c+1),\quad
c(11c^2+5c+3),\quad
-2c^3+c^2+c+1,\quad3c(3c+1)
\right),
\]
\[
B_{3,*}=\left(
c^2(3c^2+4c+1),\quad
c(8c^3-c^2-3),\quad
c(5c^2+3c+3),\quad
12c^3+9c^2+3c+3
\right).
\]

The trace vanishes identically. The only finite poles occur among the bad
parameters `0`, `-1`, and the three roots of `H`. The polynomial discriminant
of `Q` is

\[
\operatorname{Disc}_u Q=256c^5(c+1)H(c)^2.
\]

In particular, this matrix has no additional finite apparent poles in the
chosen basis. Trace zero here is a verified property of this basis; it does
not imply unitary transport or conservation of the Hodge norm.

Two exact specializations are

\[
A_D(-2)=
\begin{pmatrix}
3/16&-11/8&7/16&9/16\\
3/8&41/16&-5/8&-15/16\\
-7/8&-37/8&19/16&15/8\\
5/4&71/8&-17/8&-63/16
\end{pmatrix},
\]

\[
A_D(-3)=
\begin{pmatrix}
19/88&-17/88&19/264&9/176\\
9/88&131/264&-1/11&-1/11\\
-19/44&-87/88&61/264&3/11\\
6/11&57/22&-39/88&-83/88
\end{pmatrix}.
\]

## 4. Exact certificate, including varying pullbacks

Every row is certified by a polynomial `R_j(u,c)` of degree at most four,
whose rational coefficient functions are included in the JSON output:

\[
-\tfrac12u^jQ_c
=Q\sum_{k=0}^3 A_{D,jk}u^k+
R'_jQ-\tfrac12R_jQ_u.
\]

Dividing by `v^3` and multiplying by `du` proves

\[
\partial_c\eta_j=\sum_k A_{D,jk}\eta_k+d_u(R_j/v).
\]

Thus the certificate is independent of numerical quadrature or a branch
choice. The program checks all nine coefficients of each identity in the
exact rational-function field `Q(c)`.

For the elliptic basis `(du/w,u^2 du/w)`, the inherited matrix is independently
certified as

\[
A_C(c)=\begin{pmatrix}
-1/(4c)&1/[4c(c+1)]\\-1/4&1/(4c)
\end{pmatrix}.
\]

The genus-three matrix in the six-form basis above is precisely

\[
A_{C_3}=\operatorname{diag}(A_C,A_D).
\]

Because the quotient maps themselves depend on `c`, it would be unsafe to
verify this by substitution alone. The program differentiates each pulled
form and certifies all six identities directly on `C_3`. Its exact primitive
contains the appropriate additional contraction term from `du/dc=1`.

Concretely, if the elliptic reduction primitive is `R/w`, the genus-three
primitive is `[R(x^2+c)+(x^2+c)^j]/y` for `j=0,2`. In the `D_2` sector it is
`[R_j(x^2+c)+(x^2+c)^j]/(xy)`. The latter is polynomial-over-`y`, because
`R_j(c)=-c^j`, so the numerator is divisible by `x^2`.

## 5. The projectors are flat; the Hodge splitting is compatible

The involution `sigma` is an automorphism of the entire smooth family over
the parameter base. Naturality of the Gauss–Manin connection therefore gives

\[
\nabla\sigma^*=\sigma^*\nabla.
\]

Consequently the two projectors

\[
P_\pm=\tfrac12(1\pm\sigma^*)
\]

are parallel endomorphisms, and each summand is a flat subbundle. In the
displayed basis, `sigma^*=diag(1,1,-1,-1,-1,-1)`, so this is also visible
directly from the zero off-diagonal blocks of `A_{C_3}`. Holomorphic pullback
preserves the Hodge filtration. The inherited and new sectors therefore
split both the local system and its Hodge structure.

This does not say that the holomorphic subbundle itself is flat: differentiating
its forms produces second-kind classes, as the matrices explicitly show.
Nor does it define an iteration on the parameter `c`. Transport along a
chosen path is now computable, but selecting a path, base flow, or
renormalization remains additional input.

## Reproduction and scope

Run `python3 geometry/compute_transport.py` from the bundle root. It uses only
Python's standard library and writes `geometry/transport_results.json`.
It certifies four genus-two identities, two elliptic identities, six direct
genus-three identities, and the trace-zero identity. Specializations are
exact fractions. These are algebraic verification checks, not interval
certifications of analytic monodromy or Lyapunov exponents.

General mathematical background: functorial Gauss–Manin transport and Hodge
bundles are discussed in Forni–Matheus–Zorich,
[*Lyapunov spectrum of invariant subbundles of the Hodge bundle*](https://webusers.imj-prg.fr/~anton.zorich/Papers/MatheusForniZorich%20Lyapunov%20Spectrum%20ETDS%202012.pdf).
The displayed matrix and polynomial certificates are calculations performed
for this explicitly specified family; no novelty claim is made.
