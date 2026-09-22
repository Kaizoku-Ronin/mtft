# Triangle involution, Feigenbaum scaling, and a Hodge tower

Research note · 19 September 2026 · MTFT v0.32.0 source snapshot

## Result

The triangle reflection can be carried into an explicit tower of algebraic curves associated with quadratic iteration. The tower has exact maps, an exact decomposition of its holomorphic differentials, and singular fibers controlled by the critical orbit. Its first elliptic member also gives an explicit degree-two isogeny and a rational parametrization of the classical modular curve X₀(2).

This answers a concrete version of the question “can Hodge structures build on themselves as the algebra scales?” They can in this chosen tower: each double cover retains the previous Hodge structure and adds a Prym component. The construction does **not** derive the quadratic iteration, its parameter, a physical coupling, or X₀(143) from the triangle invariant.

We also computed the main period-doubling cascade through primitive period 1024. Its finite ratios approach the two usual Feigenbaum constants. These computations are high-precision numerical diagnostics, not interval certificates and not a new proof of universality.

## 1. What is already present in MTFT

The source files `src/mtft/combinatorial.py` and `src/mtft/quadratic_forms.py` share the exact identity

\[
T(n)=\frac{n(n+1)}2,\qquad \sigma(n)=-1-n,\qquad
T(\sigma(n))=T(n),\qquad (2n+1)^2=8T(n)+1.
\]

In the coordinate u=2n+1, the reflection is u↦−u. Thus

\[
\mathbb Q[n]^\sigma=\mathbb Q[T],\qquad
\mathbb Q[n]=\mathbb Q[T]\oplus(2n+1)\mathbb Q[T].
\]

Applied to three integer coordinates, the same change of variables gives

\[
N=T(a)+T(b)+T(d)
\iff 8N+3=(2a+1)^2+(2b+1)^2+(2d+1)^2.
\]

The separate self-similarity of the three-square forbidden indicator, F(4n)=F(n), follows from the obstruction n=4ᵃ(8b+7). It is a different statement from invariance under reflection. Neither statement defines time evolution.

The repo’s own combinatorial introduction describes the analogy with Hodge eta-parity as structural, explicitly withholding a derivation link. This investigation constructs an additional family in which the reflection/Hodge relationship can be made explicit.

## 2. Exact quadratic dynamics after an additional choice

Extend T to a real polynomial and choose to iterate

\[
F_r(n)=2rT(n)=rn(n+1).
\]

With x=−n, this is exactly the logistic map x↦rx(1−x). For r≠0, the further affine coordinate z=r(n+1/2) gives

\[
z\longmapsto f_c(z)=z^2+c,\qquad c=\frac{r(2-r)}4.
\]

The relevant real inverse branch is r=1+√(1−4c). The triangle reflection becomes z↦−z. Notice the precise identity: f_c(−z)=f_c(z). This is a folding symmetry of the input; it is not equivariance f_c(−z)=−f_c(z).

The parameter r and the instruction to iterate F_r have been added. Iterating the unscaled triangle polynomial itself corresponds to r=1/2. On n∈[−1,0], it tends toward 0; on the nonnegative integers, 0 and 1 are fixed and n≥2 grows. The literal triangle polynomial does not select the period-doubling regime.

Feigenbaum universality is an established renormalization result for suitable families and transverse parameter changes. Lyubich identifies the parameter scaling constant with an unstable renormalization eigenvalue and proves universality under stated hypotheses.[1] An arbitrary reflection-invariant map need not satisfy these hypotheses.

## 3. A natural curve tower relative to the chosen map

Write

\[
P_m(z,c)=f_c^{\circ m}(z),\qquad p_m(c)=P_m(0,c).
\]

The proposed family is the smooth compactification, at good parameters, of

\[
\boxed{C_{m,c}:\ y^2=P_m(z,c).}
\]

This is a definite construction, not a uniqueness claim: other curves can be attached to a dynamical system. Here the branch points are precisely the depth-m preimages of the critical point 0, which makes the choice particularly transparent.

For generic c, P_m has degree 2ᵐ and distinct roots. The double cover of the projective z-line therefore has

\[
g_m=2^{m-1}-1,\qquad h^{1,0}=h^{0,1}=g_m.
\]

For m≥2 there is an actual degree-two map

\[
\pi_m:C_{m,c}\to C_{m-1,c},\qquad (z,y)\mapsto(z^2+c,y).
\]

Its deck involution is σ(z,y)=(−z,y). It is different from the hyperelliptic involution ι(z,y)=(z,−y).

At common good parameters the map has four branch points: two at z=0 and two at infinity. Riemann–Hurwitz gives g_m=2g_{m−1}+1. A basis of holomorphic differentials is

\[
\omega_j=z^j\frac{dz}{y},\qquad 0\le j\le g_m-1,
\qquad \sigma^*\omega_j=(-1)^{j+1}\omega_j.
\]

The invariant forms are the pullbacks of the previous curve’s forms:

\[
\pi_m^*\!\left(w^k\frac{dw}{y}\right)
=2z(z^2+c)^k\frac{dz}{y}.
\]

Their distinct odd leading powers prove independence and identify the whole invariant subspace. Consequently,

\[
\dim H^{1,0}(C_m)^+=g_{m-1},\qquad
\dim H^{1,0}(C_m)^-=g_{m-1}+1.
\]

| Iteration depth m | Genus gₘ | Inherited holomorphic forms | New holomorphic forms |
|---:|---:|---:|---:|
| 2 | 1 | 0 | 1 |
| 3 | 3 | 1 | 2 |
| 4 | 7 | 3 | 4 |
| 5 | 15 | 7 | 8 |
| 6 | 31 | 15 | 16 |

Over rational cohomology, this is a decomposition into Hodge substructures. Equivalently, the Jacobian is isogenous to J(Cₘ₋₁) times a Prym variety of dimension gₘ₋₁+1. The norm–pullback composition is multiplication by 2. This is standard double-cover/Prym theory applied to this explicit tower.[2] It is an isogeny statement, not an integral product identity or an assertion that the induced Prym polarization is principal; these covers have four branch points.

The new part also has an explicit curve model. Quotienting by the product involution τ=σι, rather than σ, gives

\[
D_{m-1,c}:v^2=(u-c)P_{m-1}(u,c),\qquad
(u,v)=(z^2+c,zy).
\]

Its generic genus is gₘ₋₁+1. Its differentials pull back as

\[
u^k\frac{du}{v}\longmapsto2(z^2+c)^k\frac{dz}{y},
\]

which span the σ-odd subspace. Thus the full decomposition can be written

\[
\boxed{J(C_{m,c})\sim J(C_{m-1,c})\times J(D_{m-1,c}).}
\]

Here ∼ denotes isogeny, a surjective homomorphism with finite kernel. These maps are defined over the parameter field; Hodge structures are considered after passing to complex fibers. The new component is therefore explicitly computable, not just a dimension count.

## 4. Critical orbit returns control singular fibers

The discriminant identity is

\[
\boxed{\operatorname{Disc}_z P_m
=(-1)^{2^{m-1}(2^m-1)}\,2^{m2^m}
\prod_{j=1}^m p_j(c)^{2^{m-j}}.}
\]

Proof: P′ₘ=2ᵐ∏ₖ₌₀ᵐ⁻¹Pₖ, where P₀=z. Use multiplicativity of the resultant and the fact that Pₘ at any root of Pₖ equals pₘ₋ₖ. The monic polynomial discriminant is (−1)ᴰ⁽ᴰ⁻¹⁾/² Res(Pₘ,P′ₘ), D=2ᵐ. This proves the displayed formula for all m; the program independently checks 28 specializations by Sylvester determinants and fraction-free elimination.

For example,

\[
\operatorname{Disc}P_1=-4c,\qquad
\operatorname{Disc}P_2=256c^3(c+1),\qquad
\operatorname{Disc}P_3=2^{24}p_1^4p_2^2p_3.
\]

If the critical point has primitive period L≥2, then p_L(c)=0 while p_j(c)≠0 for 1≤j<L. The first singular curve is C_L. Its unique double root is z=0, with expansion

\[
P_L(z,c)=2^{L-1}\!\left(\prod_{j=1}^{L-1}p_j(c)\right)z^2+O(z^4).
\]

Thus it has one ordinary node. Moreover p′_L(c)≡1 modulo 2: the recursion pₗ₊₁=pₗ²+c proves this immediately. Since p_L is monic and squarefree modulo 2, it is squarefree in characteristic zero. The parameter crosses the nodal degeneration transversely. Standard Picard–Lefschetz theory then gives unipotent monodromy and a logarithmic dual period. Not every normalized period diverges; the other local period can remain finite.

Because Cₗ₋₁ is still smooth at such a parameter, the new degeneration belongs to the added Hodge component. The quotient model Dₗ₋₁ displays the collision between u=c and a root of Pₗ₋₁. Isogenies can change integral cycle normalizations, so this observation does not identify their integral monodromy matrices without further choices.

**Where the infinite cascade lives.** A superstable center of primitive period 2ⁿ first makes C₂ⁿ singular. The corresponding generic genera begin 0, 1, 7, 127, 32767 for periods 1, 2, 4, 8, 16. Iteration depth m and period-doubling generation n must not be confused.

At the accumulation parameter c∞, the critical orbit is not periodic. Therefore every finite P_m is squarefree there and every finite C_m is smooth. The accumulation is visible across the unbounded tower, not as a singularity of one fixed member. A generically smooth finite algebraic family over the parameter line has only finitely many bad fibers; it cannot have an infinite set of isolated bad parameters accumulating inside that line.

## 5. The first elliptic member gives X₀(2)

For m=2,

\[
C_{2,c}:y^2=z^4+2cz^2+c^2+c,\qquad c\ne0,-1.
\]

The monic quartic has rational points at infinity, so choosing one makes it an elliptic curve over the parameter field. Its j-invariant is

\[
j(C_{2,c})=64\frac{(4c+3)^3}{c+1}.
\]

One derivation uses the binary quartic invariants I=4c(4c+3), J=16c²(8c+9), the polynomial discriminant 256c³(c+1), and j=256I³/Disc.

The product involution τ=σι sends (z,y) to (−z,−y). It has no fixed points at good parameters. Its quotient is a degree-two elliptic isogeny, represented by

\[
(u,v)=(z^2,zy),\qquad
E'_c:v^2=u(u^2+2cu+c^2+c),
\]

with

\[
j(E'_c)=64\frac{(3-c)^3}{(c+1)^2}.
\]

This is distinct from the quotient by σ alone, whose target C₁ has genus zero. Choosing compatible origins turns the quotient map by τ into a group homomorphism.

Set t=−64(c+1). The pair becomes

\[
\boxed{j=\frac{(t+16)^3}{t},\qquad
j'=\frac{(t+256)^3}{t^2}.}
\]

These are the classical rational functions parametrizing X₀(2), the moduli curve of elliptic curves with a cyclic degree-two isogeny. The explicit quotient establishes the isogeny, rather than relying on numerical proximity of j-values. The code also checks the classical modular polynomial Φ₂(j,j′)=0 exactly on six rational parameters. The Fricke transformation interchanging the isogeny and its dual is

\[
t\mapsto4096/t,\qquad c\mapsto-\frac{c}{c+1}.
\]

It has order two. It is not the period-doubling renormalization operator. X₀(2) is also not X₀(143); no correspondence to the latter is supplied by this construction.

## 6. A period that can actually be computed

Approach c=−1 from below, writing c=−1−ε, ε>0, s=√(1+ε), a²=s(s+1), b²=s(s−1). Then

\[
P_2=(z^2-a^2)(z^2-b^2),\qquad
I(c)=\int_b^a\frac{dz}{\sqrt{(a^2-z^2)(z^2-b^2)}}
=\frac{K(k)}a,\quad k^2=1-b^2/a^2.
\]

I is the magnitude of a branch-to-branch integral, half a closed period up to phase. The standard elliptic-integral asymptotic gives

\[
\boxed{I(-1-\epsilon)
=\frac{\log(1/\epsilon)}{2\sqrt2}+\frac{\log8}{\sqrt2}+o(1).}
\]

We evaluate K using its arithmetic–geometric mean formula and repeat at 60 and 90 Decimal digits.[3,4] For ε from 10⁻² to 10⁻³⁰, the largest precision-rerun discrepancy is below 10⁻⁵⁸. The final secant slope against log(1/ε) is approximately 0.353553390593273762203, approaching the exact coefficient 1/(2√2). This verifies the normalization of the local logarithmic period numerically. It does not identify that coefficient with a Feigenbaum constant.

The modulus used by this elliptic-integral expression need not be the direct branch cross-ratio used to compute j(C₂); quadratic elliptic transformations intervene. At c=0 the quartic is singular although j tends to 1728, so finite j alone is not a smoothness test for a chosen model.

## 7. Numerical period-doubling check

The scaling program solves p₂ⁿ(cₙ)=0 on the principal real cascade through n=10. It distinguishes these superstable centers from bifurcation boundaries. It avoids known lower-period roots and checks the halfway return is numerically nonzero. A sign-change scan is not a certified enumeration of all roots.

\[
\delta_n=\frac{c_{n-1}-c_{n-2}}{c_n-c_{n-1}},\qquad
d_n=f_{c_n}^{\circ2^{n-1}}(0),\qquad
\alpha_n=\left|\frac{d_{n-1}}{d_n}\right|.
\]

At primitive period 1024:

| Quantity | Numerical result |
|---|---:|
| c₁₀ | −1.40115478254661784122 |
| Corresponding logistic r₁₀ | 3.56994535548646858089 |
| δ₁₀, using the c parameter | 4.66919515603001717402 |
| α₁₀ magnitude | 2.50290757151227907728 |

The limiting comparison values are approximately 4.66920160910299 and 2.50290787509589. The final finite-order discrepancies are about 6.45×10⁻⁶ and 3.04×10⁻⁷, respectively. They are much larger than the 60/90-digit arithmetic rerun discrepancies, about 2.21×10⁻⁴² and 5.10×10⁻⁴³. Neither small residuals nor precision stability prove the infinite limit.

Finite c-parameter ratios differ from finite r-parameter ratios under the nonlinear change c=r(2−r)/4. Their asymptotic limits agree because this change is locally nonsingular at the accumulation point.

## 8. What this establishes, and what it leaves open

| Statement | Status |
|---|---|
| Triangle reflection and odd-square change of coordinates | Exact; already in MTFT |
| Fᵣ=2rT is conjugate to the logistic family | Exact after choosing this family |
| Cₘ has double-cover maps and inherited/Prym Hodge pieces | Exact; explicit construction plus standard theory |
| Critical orbit returns determine the discriminant | Proved by resultant factorization here |
| First elliptic member carries a degree-two isogeny | Exact; explicit quotient and j-functions |
| Nodal period has the stated logarithmic asymptotic | Standard theory plus explicit elliptic reduction; numerically checked |
| Finite scaling ratios approach known Feigenbaum values | High-precision numerical evidence here; universality proved in existing literature |
| MTFT’s operators select this r or its accumulation value | Not established by the inspected source |
| Feigenbaum δ is a local Picard–Lefschetz eigenvalue | Not supplied: local nodal monodromy is unipotent |
| The construction derives electromagnetism, gravity, or X₀(143) | Not established |

The next substantive research target is a compatible action of dynamical renormalization on this tower’s periods or correspondences. The tower maps by themselves increase iteration depth; they are not yet an operator whose unstable eigenvalue is δ. Any proposed action must state its domain, normalization, parameter map, and compatibility with pullback/norm. A physical interpretation would additionally have to select the dynamics from the existing MTFT structure.

## 9. Source audit and legacy issues

The read-only audit found genuine existing dynamics: `surface/dynamics.py` builds linear Hamiltonian flows from supplied energy matrices, and `periods/hamiltonian.py` constructs specified graph-derived quadratic Hamiltonians. No triangle-to-logistic selector was found in the inspected modules. `info_geometry.logistic_iterate` accepts r as an argument; its accumulation parameter and Feigenbaum constants are literals. The conditional parameter selection in `surface/oldsector.py` concerns a different operator family and cannot be substituted for r without a new derivation.

The literal v0.32.0 `info_geometry.py` still labels a nonzero expression as the Ricci scalar of a one-dimensional Fisher manifold. Intrinsic Ricci curvature in one dimension is zero. That legacy expression was not used here. No CC-18 pointer was found in this sdist; a correction may exist in other project artifacts. No repo corrections or version changes were made.

The supplied Wikipedia PDF also drops a factor of two from its cited asymptotic relation: Delbourgo, Hart, and Kenny report 3δ∼2α² for the period-multiplying context, not 3δ∼α².[5] The repo’s separate product conjecture evaluates numerically to δ²logα≈20.0018043867, rather than exactly 20. Neither approximation was used to generate the cascade or the curve tower.

## Reproduction and provenance

Run from this bundle’s root:

```bash
python3 scaling/compute_scaling.py
python3 geometry/verify_geometry.py
python3 plot_results.py
```

The first two programs use only the Python standard library; plot regeneration additionally needs Matplotlib. Results, conventions, and classification are stored in adjacent CSV/JSON files. The source snapshot is identified by the manifest; the calculations do not require installing MTFT. The algebraic proofs above supply the global claims; finite tests are supporting checks rather than substitutes for proof.

Reviewed archive: `mtft-0.32.0.tar.gz`.

SHA-256: `46ffc563a0cd81a7a82d975497bf64d82f8471eb424409e3d154532b5d8d6069`.

## References

1. M. Lyubich, *Feigenbaum–Coullet–Tresser universality and Milnor’s Hairiness Conjecture*, Annals of Mathematics 149 (1999), especially the universality statement and renormalization eigenvalue interpretation. [Author paper](https://arxiv.org/abs/math/9903201).
2. D. Mumford, *Prym Varieties I* (1974), double covers, norm/pullback, even/odd decomposition and polarization qualifications. [Author-hosted paper](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1974a--PrymVar-I-NC.pdf).
3. NIST Digital Library of Mathematical Functions, §19.8, quadratic transformations and the arithmetic–geometric mean. [DLMF](https://dlmf.nist.gov/19.8).
4. NIST Digital Library of Mathematical Functions, §19.12, elliptic-integral asymptotics. [DLMF](https://dlmf.nist.gov/19.12).
5. R. Delbourgo, W. Hart and B. G. Kenny, *Dependence of universal constants upon multiplication period in nonlinear maps*, Physical Review A 31 (1985), 514–516. [Original paper](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.31.514).

The tower specialization and elementary formulas were derived in this investigation using standard mathematics. No claim of historical novelty or new physical law is made.
