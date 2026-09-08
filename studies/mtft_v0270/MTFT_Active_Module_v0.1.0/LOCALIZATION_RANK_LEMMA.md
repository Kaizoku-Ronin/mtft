# Localized harmonic support bounds the active module

8 September 2026. An elementary linear-algebra result extracted during the
active-module audit. No novelty or priority claim is made.

## Proposition

Let W be a real full-row-rank 2g by m harmonic matrix, Gg=W W^T, and let G
be a positive definite metric with a compatible complex structure J:
J^2=-I and J^T G J=G. Let diagonal masks D_i have support in an edge set S,
and let the columns W[:,S] have rank r. Choose a basis H for their span.

Define the localized operators and their two successive projections by

\[
V_i=G_g^{-1} W D_i W^T,\quad
A_i=\tfrac12(V_i+V_i^{\dagger_G}),\quad
K_i=\tfrac12(A_i+J A_i J),\qquad
V_i^{\dagger_G}=G^{-1}V_i^T G.
\]

Then every K_i has its image in

\[
E=\operatorname{span}_{\mathbb R}
\{G_g^{-1}H,\ G^{-1}H,\ JG_g^{-1}H,\ JG^{-1}H\},
\qquad \dim_{\mathbb R} E\le 4r.
\]

The common kernel of the K_i contains the J-stable space E^{perp_G}.
Its complex dimension is therefore at least max(0,g-2r). The same kernel
bound holds for all commutators [K_i,K_j] and the Lie algebra they generate.
Thus the complex dimension of the corresponding active sector is at most 2r
(and, trivially, at most g).

## Proof

The support assumption gives im(V_i) contained in im(Gg^-1 H). Transposing
V_i shows im(V_i^{dagger_G}) contained in im(G^-1 H). The image of A_i is
therefore in the span of these first two sets of columns. Multiplication on
the left by J adds only the last two sets; multiplication on the right does
not increase the image. This proves the displayed range bound.

Compatibility implies J^{dagger_G}=-J. Therefore K_i is G-selfadjoint and
K_i J=-J K_i. The space E is J-stable because J^2=-I. Selfadjointness and
im(K_i) contained in E imply K_i vanishes on E^{perp_G}: for x in that
complement and every y, <K_i x,y>_G=<x,K_i y>_G=0.

The complement is J-stable since J is a G-isometry. Its real dimension is
at least 2g-4r, hence its complex dimension is at least g-2r whenever this
is positive. Commutators of maps annihilating the same subspace also
annihilate it, and this property survives linear combinations and further
commutators. This proves the claim.

## Application to the supplied MTFT construction

The unchanged triangle choices are (1,11,12). They touch the zero-indexed
edge set [1,2,11,12,13,59,81,82]. On the packaged exact rational harmonic
matrix these columns have rank r=4, certified by a nonzero 4x4 minor and
exact reconstruction coefficients in `kernel_check_results.json`.
The independent columns are [2,11,12,13]; one nonzero minor has determinant
-10610616427/14658673262330686.

For g=13 the proposition gives at least five complex inactive directions,
before any numerical Lie closure is performed. The 26x16 matrix spanning E
has numerical rank 16, and its orthogonal projector agrees with the measured
active projector to 4.90e-13 in the orthonormal real frame. The corresponding
common kernel has numerical complex dimension five.

**Evidence classes:** the support-rank certificate is EXACT on the supplied
rational data. The proposition is exact under its displayed assumptions.
Saturation of the period-dependent span, projector equality, and the result
of substituting the numerical Hodge inputs remain DIAGNOSTIC. Neither the
rank bound nor its saturation identifies the Lie algebra as D4 or assigns a
Clifford/spacetime signature.
