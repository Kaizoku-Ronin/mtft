# Scratch derivation: quadratic iteration, Hodge tower, and critical degenerations

This is an intermediate mathematical audit, not an MTFT implementation or a claim of novel theorems. All statements below concern the explicitly chosen family f_c(x)=x²+c. The triangle invariant alone does not select this family or its parameter.

## 1. Definitions and exact discriminant

Let P_m(x,c)=f_c^{∘m}(x), with P_0=x, and p_m(c)=P_m(0,c), so p_0=0 and p_(m+1)=p_m²+c. Put d=2^m.

The chain rule gives

    ∂_x P_m = 2^m ∏_(k=0)^(m−1) P_k.

If P_k(z)=0, then P_m(z)=p_(m−k). Resultants therefore give

    Disc_x(P_m) = (−1)^[2^(m−1)(2^m−1)]
                  · 2^[m·2^m] · ∏_(j=1)^m p_j(c)^[2^(m−j)].

Proof: Disc(P)=(-1)^[d(d-1)/2] Res(P,P') for monic P. Multiplicativity of the resultant and Res(P_m,P_k)=p_(m−k)^[2^k] give the formula. The interchange sign of these two resultants is + because deg P_m is even.

Explicitly:

* m=1: Disc(P_1)=−4c.
* m=2: P_2=x⁴+2cx²+c²+c, Disc(P_2)=2⁸c³(c+1).
* m=3: Disc(P_3)=2²⁴c⁷(c+1)²(c³+2c²+c+1).

Every p_m is squarefree over Q: p_m'(c)≡1 modulo 2, and p_m is monic integral, so a nonconstant monic common factor with p_m' would reduce to a nonconstant divisor of 1 modulo 2.

Its roots are exactly critical parameters of periods dividing m. The exact-period factors are the Gleason polynomials G_m, with p_m=∏_(e|m)G_e. Root multiplicity is never needed to remove lower periods: remove the lower-period factors explicitly.

## 2. A genuine tower of curves and Hodge structures

For c outside the discriminant zero set, let C_(m,c) be the smooth compactification of

    y²=P_m(x,c).

The double cover of the x-line has 2^m distinct finite branch points and none at infinity. Riemann–Hurwitz gives

    g_m=2^(m−1)−1,   m≥1.

Thus m=1,2,3,4,5 gives genera 0,1,3,7,15. Holomorphic forms are

    ω_j=x^j dx/y,  0≤j≤g_m−1.

There is an exact morphism

    π_m:C_(m+1,c)→C_(m,c),   (x,y)↦(x²+c,y),

because P_m(x²+c,c)=P_(m+1)(x,c). It has degree two, with deck involution

    σ(x,y)=(−x,y).

This is different from the hyperelliptic involution ι(x,y)=(x,−y). The latter acts as −1 on every holomorphic one-form and every H¹ class. The former has both signs and creates the old/new splitting.

At generic common-good c, π_m has four ramification points: the two points with x=0 and the two points at infinity. Hence g_(m+1)=2g_m+1.

On C_(m+1), σ*ω_j=(−1)^(j+1)ω_j. Its + eigenspace has complex holomorphic dimension g_m; its − eigenspace has dimension g_m+1=2^(m−1). Pullback of base forms is

    π_m*(u^j du/y)=2x(x²+c)^j dx/y,

and spans the + eigenspace. Transfer satisfies π_(m,*)π_m*=2. Thus, after rational coefficients,

    H¹(C_(m+1),Q) = π_m*H¹(C_m,Q) ⊕ H¹(C_(m+1),Q)^−

as Hodge structures at each complex good parameter. The anti-invariant piece has real/rational dimension 2(g_m+1).

On Jacobians, the norm map has connected-kernel component P_m=Prym(C_(m+1)/C_m). Then

    J(C_(m+1)) ~ J(C_m) × P_m,    dim P_m=g_m+1.

Here ~ means isogeny, not integral product or equality of principally polarized varieties. The maps, connected abelian kernel, and isogeny are defined over Q(c), and descend over char-zero specialization fields wherever both curves are smooth. Betti Hodge structures are considered after a complex specialization. The four-branch-point Prym does not automatically inherit the principal polarization familiar from unramified/two-branch-point Pryms.

This is a precise mathematical realization of Hodge structures growing in a tower while an involution splits inherited and new parts. It is natural for the chosen iteration construction, not established as selected by MTFT.

## 3. An explicit curve for the new Prym factor

The other quotient uses τ=σι:(x,y)↦(−x,−y). Put u=x²+c, v=xy. It gives

    D_(m,c): v²=(u−c)P_m(u,c),

of genus 2^(m−1)=g_m+1. The map C_(m+1)→D_m is degree two. Its pulled-back holomorphic forms are

    u^j du/v ↦ 2(x²+c)^j dx/y,

spanning exactly the σ-odd part. Consequently J(D_m) is isogenous to P_m, and

    J(C_(m+1)) ~ J(C_m) × J(D_m).

Beware degeneration multiplicities: Disc[(u−c)P_m(u)]=Disc(P_m)·p_(m+1)². Thus D_m may have a quadratic smoothing parameter where C_(m+1) has a transverse node. Isogenies preserve rational Hodge structures but can change integral monodromy normalization.

## 4. Primitive superstable periods are first nodal fibers

Suppose c_0 makes 0 have exact period L≥2. Then p_L(c_0)=0, all p_j(c_0)≠0 for 1≤j<L, and p_L'(c_0)≠0. At m=L,

    P_L(x,c_0)=A_L x²+O(x⁴),
    A_L=2^(L−1)∏_(j=1)^(L−1)p_j(c_0)≠0.

The only repeated root is x=0, with multiplicity two. Therefore C_L has exactly one ordinary node and normalization genus g_L−1. Locally the family is a transverse smoothing y²−A_Lx²≈p_L'(c_0)(c−c_0). It has the standard Picard–Lefschetz monodromy and one logarithmic entry after an appropriate choice of normalized period matrix. It is incorrect to say every period diverges: one cycle pinches, while an appropriately dual period has logarithmic behavior.

C_(L−1) remains smooth at c_0. Thus this first degeneration lies in the new anti-invariant/Prym Hodge sector.

Do not conflate map iteration m with period-doubling generation k. For period L=2^k, the first detecting fiber in this tower is C_(2^k), whose genus is

    2^[2^k−1]−1.

Periods 2,4,8,16 correspond to genera 1,7,127,32767. The explicit curve becomes large very quickly even though critical orbit recurrence remains cheap to evaluate.

At the nonperiodic Feigenbaum accumulation parameter c_infinity, every finite p_m(c_infinity) is nonzero, hence every finite C_m is smooth. The accumulation occurs only across an unbounded tower. More generally a generically smooth algebraic family of curves over an algebraic parameter curve has only finitely many bad fibers; a single finite algebraic family cannot encode infinitely many distinct singular parameters accumulating at an interior point.

## 5. Exact elliptic example and logarithmic period

C_2 is the elliptic quartic

    y²=x⁴+2cx²+c²+c.

Its binary quartic invariant I=4c(4c+3) and polynomial discriminant Δ=256c³(c+1) give

    j(C_2)=256 I³/Δ=64(4c+3)³/(c+1).

At c=−1 there is a simple j pole and the nodal fiber y²=x²(x²−2). At c=0 the quartic is singular but j has finite limit 1728. This is potentially good reduction, not a contradiction: use c=t⁴, x=tX, y=t²Y to obtain

    Y²=X⁴+2t²X²+t⁴+1,

which is smooth at t=0 with j=1728. A finite j value alone does not make the original equation smooth.

For a real period approach c=−1−ε with ε>0, let s=√(−c)>1, a²=s(s+1), b²=s(s−1). Then

    P_2=(x²−a²)(x²−b²),   0<b<a,

and the positive branch-to-branch integral is

    I(c)=∫_b^a dx/√[(a²−x²)(x²−b²)]
        = K(k)/a,    k²=1−b²/a²=2/(s+1).

Twice this integral, up to multiplication by i and orientation, is a closed period of dx/y. Since K(k)=log(4/k')+O(k'²log(1/k')) and k'=b/a≈√ε/2,

    I(−1−ε)=−log ε/(2√2)+log8/√2+o(1).

Equivalently I=π/[2 AGM(a,b)], a stable way to evaluate it. A direct four-branch-point Legendre cross ratio is λ=((a−b)/(a+b))². Beware: the k² in the simple integral representation is related by a degree-two Landen transformation and is not itself this λ. Indeed

    I=2 K(√λ)/(a+b).

## 6. Exact level-two modular correspondence

On C_2, τ(x,y)=(−x,−y) is fixed-free at good parameters. Its elliptic quotient, with (u,v)=(x²,xy), is

    E'_c: v²=u(u²+2cu+c²+c).

After choosing compatible origins the quotient is a degree-two isogeny. The Weierstrass invariants of E' give

    j(E')=64(3−c)³/(c+1)².

Set t=−64(c+1). Then

    j(C_2)=(t+16)³/t,
    j(E')=(t+256)³/t².

These are the standard rational parametrization of X_0(2). The Fricke involution t↦4096/t becomes c↦−c/(c+1). It exchanges the j invariants of the two-isogenous pair. This is an exact modular correspondence, but it is not Feigenbaum's renormalization operator or a parameter selection principle.

## 7. Sources actually retrieved

* Mumford, Prym Varieties I (1974), especially pages 3, 6–7: norm/pullback identity, odd/even decomposition, and finite-two-torsion isogeny splitting. https://www.dam.brown.edu/people/mumford/alg_geom/papers/1974a--PrymVar-I-NC.pdf
* Bloch, de Jong, Sertöz, Heights on curves and limits of Hodge structures (2023), section 3, formulas 3.14–3.18 and 3.29: explicit logarithmic nodal period behavior and limiting mixed Hodge structures. https://arxiv.org/abs/2206.01220
* Baker, Chen, Li, Qian, Necklaces, permutations, and periodic critical orbits for quadratic polynomials (2025): Gleason polynomials and primitive critical periods. https://arxiv.org/abs/2508.12924
* NIST DLMF 19.12: elliptic-integral logarithmic asymptotic. https://dlmf.nist.gov/19.12
* NIST DLMF 19.8: AGM and Landen transformations. https://dlmf.nist.gov/19.8

All concrete polynomial, tower, dimension, quotient, discriminant and j calculations above were derived independently here; the sources support the general established machinery. They are not claims of new mathematics or proof of a physical interpretation.
