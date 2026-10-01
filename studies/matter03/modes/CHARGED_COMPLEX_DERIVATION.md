# Matter-03: charged Higgs-complex modes

The specified paired neutral Higgs representative does not select three light charged families. The two-term holomorphic deformation complex leaves **237 degree-one hypercohomology classes of each electric charge sign**. This is an exact conditional geometric count; it is not a physical mass spectrum. Its identification with physical zero modes additionally requires a Hermitian BPS vacuum and the chosen bulk parent dictionary.

The isolated Matter-02 18 × 44 Yukawa slice does have rank at least three. It is only one block of the complete charged operator. Its rectangular kernel difference 26 is not an electric-charge chirality index.

## Fixed representative and charge convention

Use the Matter-02 background

\[
S=\mathrm{Bl}_{(P,O)}(X_0(143)\times E),\quad
\mathcal L=\pi^*\operatorname{pr}_C^*\mathcal O_C(P)(-3D),\quad P=\infty.
\]

The U(8) component exponents and electromagnetic weights are

| Component | Exponent n | Q |
|---|---:|---:|
| c₁,c₂,c₃ | 0 | 1/6 |
| L+ | 1 | +1/2 |
| L− | 1 | −1/2 |
| A₁ | 2 | −1/2 |
| A₂ | −3 | −1/2 |
| B | −1 | +1/2 |

A root i ← j has q=nᵢ−nⱼ and Q=Qᵢ−Qⱼ. The two neutral Higgs entries are

\[
\Phi=s\,\phi_d E_{L+,B}+t\,\phi_u E_{B,L+},\qquad s,t\ne0.
\]

Here φ_d is the Matter-02 q=2 representative with canonical differential of order three at P; relative to its base frame its coefficient begins z⁵. Choose φ_u in H⁰(K_Sℒ⁻²) with canonical differential of order two; its relative base coefficient begins 1. Both sections exist by the exact ordinary-cusp certificate. Their nonzero scales are inputs to this count. No a-field expectation value is added.

## Full split charged spaces

In the chosen bulk convention the a and Φ chiral spaces are H¹(ℒᑫ) and H⁰(K_Sℒᑫ). For q=0 the latter two have dimensions 14 and 13; its H⁰(ℒ⁰)=1 gauge generator is separate. The complete charge-minus-one roots are

| Root | q | a dimension | Φ dimension |
|---|---:|---:|---:|
| L− ← L+ | 0 | 14 | 13 |
| A₁ ← L+ | 1 | 18 | 12 |
| A₂ ← L+ | −4 | 72 | 9 |
| L− ← B | 2 | 32 | 10 |
| A₁ ← B | 3 | 55 | 8 |
| A₂ ← B | −2 | 25 | 11 |
| Total | | **216** | **63** |

The conjugate charge-plus-one orientations give a dimension 217 and Φ dimension 62. Thus the complete split chiral totals are **279 per sign**, not 18 and 44. In particular N(−3)=44+10=54. These statements use the BHV chirality dictionary explicitly; the archived Euler-index orientation has the opposite sign.

## Deformation complex

For a holomorphic two-form Higgs background the relevant two-term complex is

\[
\operatorname{End}V\xrightarrow{\operatorname{ad}\Phi}
K_S\otimes\operatorname{End}V.
\]

The low-degree hypercohomology sequence gives

\[
0\longrightarrow\operatorname{coker}d_0\longrightarrow
\mathbb H^1\longrightarrow\ker d_1\longrightarrow0,
\qquad d_p:H^p(\operatorname{End}V)\to H^p(K_S\operatorname{End}V).
\]

Global Serre duality identifies the H¹ target with the dual a-space of opposite electric charge. Thus d₁ is the complete holomorphic aa mass pairing for the fixed Higgs background. A matrix identity in the root algebra restricts this pairing to six independent rectangular blocks:

| Spectator | Higgs entry | a exponents | Matrix size | Exact rank |
|---|---|---|---:|---:|
| L− | down | 0,−2 | 14 × 25 | 2 |
| L− | up | 2,0 | 32 × 14 | 2 |
| A₁ | down | 1,−3 | 18 × 44 | 5 |
| A₁ | up | 3,−1 | 55 × 15 | 3 |
| A₂ | down | −4,2 | 72 × 32 | 17 |
| A₂ | up | −2,4 | 25 × 87 | 12 |
| Total | | | 216 × 217 | **41** |

All ranks apply to the stated canonical orders and nonzero scales. They are holomorphic ranks, preserved by invertible kinetic normalization. They do not fix singular values or physical masses. Even without computing these ranks, the block sizes bound the total rank by 118; that weaker bound alone leaves at least 160 charged hypercohomology classes per sign.

## Why the six ranks are computable

The rank calculation reduces by global Serre duality to multiplication maps

\[
T_d(a):H^1(\mathcal L^a)\longrightarrow H^1(K_S\mathcal L^{a+2}),
\quad a=0,1,2,
\]

and

\[
T_u(a):H^1(\mathcal L^a)\longrightarrow H^1(K_S\mathcal L^{a-2}),
\quad a=2,3,4.
\]

Project the blown-up surface to its curve factor, p:S→C. All the bundles in these maps have vanishing higher direct image for the blow-down; this includes K_S, whose exceptional factor is 𝒪(D) with R¹π_*=0. We can therefore compute p-direct images on the product using the fat-point ideals. For a line bundle F_C on C and k≥0, the fat-point exact sequence gives

\[
p_*(F_C\boxtimes\mathcal O_E\otimes I_p^k)=F_C(-kP),
\]

\[
0\longrightarrow T_k\longrightarrow
R^1p_*(F_C\boxtimes\mathcal O_E\otimes I_p^k)
\longrightarrow F_C\longrightarrow0,
\qquad
T_k\simeq\bigoplus_{j=1}^{k-1}\mathcal O_{(k-j)P}\,w^j.
\]

Here the torsion length is k(k−1)/2; T₀=T₁=0. The source of the torsion is exactly the w-dependent part of the local fat-jet cokernel. Global elliptic sections are independent of w, so the w⁰ quotient is instead incorporated into F_C(−kP). The displayed extension splits as an 𝒪_C-module extension because Ext¹(F_C,T_k)=H¹(F_C⁻¹T_k)=0. This splitting is noncanonical and is not assumed to commute with multiplication.

The curve base has dimension one, so the Leray d₂ target H²_C(p_*F) vanishes. Consequently

\[
0\longrightarrow H_C^1(p_*F)\longrightarrow H_S^1(F)
\longrightarrow H_C^0(R^1p_*F)\longrightarrow0.
\]

We use the maps on this exact sequence, explicitly controlling possible off-diagonal extension contributions.

### Down Higgs: ranks 2,5,17

For a=0,1,2, the source has F_C=𝒪(aP), k_s=3a; the target has G_C=K_C((a+2)P), k_t=3a+5. The relative section coefficient is z⁵ times a unit. On torsion it injects each source module of length k_s−j into the target module of length k_s+5−j by multiplication by z⁵. Hence its torsion ranks are 0,3,15.

The free quotient H⁰_C(𝒪(aP)) is one-dimensional because P is ordinary and a≤2. Its multiplication into H⁰_C(K_C((a+2)P)) is nonzero. Therefore the complete map on H⁰(R¹p_*) is injective, of rank length(T_s)+1. This remains true for a noncanonical splitting: the torsion image and the nonzero free projection are independent.

The curve subspace map is

\[
H^1_C(\mathcal O_C(-2aP))\longrightarrow
H^1_C(K_C(-(2a+3)P)).
\]

Its target has dimension h⁰(𝒪_C((2a+3)P))=1: ordinary-cusp gaps up to seven suffice. Its Serre-dual map is multiplication by the nonzero canonical section α_d∈H⁰(K_C(−3P)); the domain of the dual is one-dimensional. Thus the curve map has rank one.

Since the map on the quotient H⁰(R¹p_*) is injective, every kernel class of the full surface map lies in the curve subspace. Possible extension cross terms cannot add rank. The full rank is precisely rank(curve)+rank(quotient), giving

\[
\operatorname{rank}T_d(0)=2,\quad
\operatorname{rank}T_d(1)=5,\quad
\operatorname{rank}T_d(2)=17.
\]

The local w-dependent rank-three submatrix in Matter-02 is contained in the a=1 rank-five map. The other two directions are global curve and elliptic contributions; they are not inferred from an arbitrary splitting.

### Up Higgs: ranks 2,3,12

At a=3,4 the source has k_s=3a and the target has G_C=K_C((a−2)P), k_t=3a−7. The relative coefficient of φ_u is a unit at P. Thus its torsion map is surjective onto T_t, with target lengths one and ten. The source free quotient has one global section, whose nonzero multiplication adds an independent rank-one free projection. The map on H⁰(R¹p_*) therefore has rank two or eleven.

The curve subspace map is

\[
H^1_C(\mathcal O_C(-2aP))\longrightarrow
H^1_C(K_C((5-2a)P)).
\]

The target dimension is h⁰(𝒪_C((2a−5)P))=1 for a=3,4. Its dual multiplication by α_u is nonzero, so the curve map is surjective onto this target subspace. Any extension cross term from a kernel class of the quotient is already in the filled curve target and cannot add a further direction. Hence the full ranks are one plus the quotient ranks: three and twelve.

At a=2 the target is K_S. Its curve and free-elliptic direct images are both K_C and there is no torsion. The curve target H¹(K_C) has dimension one, with nonzero dual multiplication by α_u. The free quotient image also has rank one. The resulting rank is two.

These injective/surjective exact-sequence arguments handle precisely the concern that filtered ranks cannot in general simply be added. No compatible splitting of the source and target extensions is assumed.

The arithmetic script constructs both canonical representatives from the archived q-expansion matrix, projects the actual rational series onto equivalent finite jet cokernels, and computes the local ranks by rational row reduction. The rank-one global contributions follow from the geometric proof and curve Serre duality. The arithmetic gates are not a formal proof-assistant certificate of sheaf cohomology.

## Gauge quotient and the remaining modes

Each electric-charge sign has one H⁰ gauge generator, the W root within the original weak U(2). Its commutator with the paired Higgs is nonzero:

\[
[\Phi,E_{L-,L+}]=-s\phi_d E_{L-,B},\qquad
[\Phi,E_{L+,L-}]=t\phi_u E_{B,L-}.
\]

Hence d₀ has rank one per sign. The Φ cokernels have dimensions 63−1=62 and 62−1=61. The a kernels have dimensions 216−41=175 and 217−41=176. The degree-one hypercohomology dimensions are therefore

\[
\dim\mathbb H^1_{Q=-1}=175+62=237,
\qquad
\dim\mathbb H^1_{Q=+1}=176+61=237.
\]

The bulk superpotential is linear in Φ. With the chosen a expectation value zero, it supplies aa pairings and the gauge quotient but no independent ΦΦ superpotential mass. Charged Φ classes therefore cannot be removed by simply normalizing the existing Yukawa tensor. Additional interactions, a changed background, or different parent data are needed to remove the excess modes. A scalar potential or fermion spectrum computed in the full Hermitian solution may organize multiplets further; this holomorphic count is the gate it must reconcile.

## Limits

This calculation is conditional on the specified holomorphic representative, parent twist and two-term complex. It neither fixes an electroweak scale nor derives a positive kinetic metric. No observed charged-family selection is established. Existence of a Hermitian BPS representative for the specified stable paired branch is established separately under its parent assumptions; the actual metric, kinetic integrals, and physically admissible gauge kernel remain open. The 237 classes are a concrete obstruction for this specified Higgs profile, not a theorem excluding other MTFT matter backgrounds. The ≥160 bound applies to any nonzero pair of sections in the same aligned weak direction.

Primary parent reference: Beasley, Heckman, Vafa, *GUTs and Exceptional Branes in F-theory—I*, arXiv:0802.3391, §§3.3.2–3.3.3 and Appendix C. New counts and map ranks are derived here from the frozen Matter-02 archive and the stated geometry.

## Review record

The projection-to-C proof was independently reviewed by the vacuum investigator and root investigator on 1 October 2026. They checked the injective down quotient and surjective up curve map, including extension cross terms. The exact rank41 and resulting237 count are approved only as conditional holomorphic statements. No Hermitian norm or observed-particle claim was reviewed or inferred.
