# Independent review of the charged aa rank41

The proposed rank41 and the resulting degree-one hypercohomology count237 per electric-charge sign are valid conditional holomorphic counts for the frozen paired representative. The argument must use the global filtration below; simply adding local jet ranks to guessed curve and elliptic contributions would not establish an exact rank because extension cross terms could otherwise contribute.

Let p:S→C be the blow-up followed by the product projection. For a line bundle F on C and k≥0, consider the sheaf on Y=C×T

\[
 F\boxtimes O_T\otimes I_{(P,Q)}^k.
\]

Its direct image to C is F(-kP). The fat-point exact sequence shows

\[
 0\longrightarrow T_k\longrightarrow R^1p_*\longrightarrow F\longrightarrow0,
 \qquad
 T_k=\bigoplus_{j=1}^{k-1} F|_{(k-j)P}\,w^j.
\]

The length of T_k is k(k−1)/2. Locally the j0 fat jets are exactly the image of F→F|_(kP); the remaining jets have positive w powers. Since F is locally free and T_k has finite support, Ext1_C(F,T_k)=H1_C(F* tensor T_k)=0. The resulting sheaf splitting is noncanonical. No functorial splitting of a multiplication map should be assumed.

Leray over the curve gives an exact filtration

\[
 0\to H^1(C,p_*A)\to H^1(S,A)\to H^0(C,R^1p_*A)\to0.
\]

Maps induced by the fixed Higgs respect this filtration. Two elementary linear-algebra facts prevent extension cancellations or extra rank: if the quotient map is injective, the full rank is the quotient rank plus the subspace rank; if the subspace map is surjective onto the target subspace, the same rank sum holds. These do not require split diagrams.

## Down maps

For a=0,1,2, put k_s=3a and k_t=3a+5. The local Higgs coefficient is z5 times a unit. Multiplication T_(ks)→T_(kt) is injective. Its source lengths are0,3,15. The induced free quotient map on H0 is nonzero from the one-dimensional H0(O_C(aP)), so the full quotient map is injective of rank length(T_ks)+1. An off-diagonal extension term does not change injectivity or this rank.

The curve subspace map is

\[
 H^1(O_C(-2aP))\to H^1(K_C(-(2a+3)P)).
\]

The target dimension is h0(O_C((2a+3)P))=1 because the certified ordinary-cusp gaps include1 through7. The dual map multiplies the nonzero canonical differential alpha_d and is therefore nonzero. This curve map has rank1. The injective-quotient fact now yields full ranks2,5,17. In particular an unexamined quotient kernel cannot feed additional image into the curve target: the quotient map has no kernel.

## Up maps

For a3 or4, put k_s=3a, k_t=3a−7. The local coefficient is a unit, so the torsion map is onto T_(kt), of lengths1 and10. The source free quotient has one-dimensional H0 and a nonzero image, giving quotient ranks2 and11.

The curve map is

\[
 H^1(O_C(-2aP))\to H^1(K_C(-(2a-5)P)).
\]

The target dimension is again1, and dual multiplication by alpha_u is nonzero, so the curve map is onto. Extension cross terms arising from quotient kernels already lie in this filled one-dimensional curve target and cannot enlarge the total rank. Thus the full ranks are3 and12. For a2 the target is K_S, whose pushforward and higher direct image are K_C and K_C. The curve map is rank1 by nonzero dual multiplication, and the free elliptic quotient contributes rank1, giving rank2.

Summing the six disjoint root blocks gives2+2+5+3+17+12=41. Their disjointness follows from the paired2×2 off-diagonal Higgs block: each spectator's two charged roots map to distinct target roots, with no cross-spectator entries. The two-term-complex long exact sequence then gives

\[
 \dim\mathbb H^1_-=(216-41)+(63-1)=237,
 \qquad
 \dim\mathbb H^1_+=(217-41)+(62-1)=237.
\]

The subtraction1 is the nonzero charged weak gauge generator map at H0. The previously certified cusp input, sheaf pushforwards, exact neutral paired representative, classical parent and deformation-complex dictionary remain assumptions. The count is neither a metric mass calculation nor observed family selection. It is a concrete excess-mode obstruction of this restricted branch.
