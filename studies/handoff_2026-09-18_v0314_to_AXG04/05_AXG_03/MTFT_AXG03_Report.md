# MTFT AXG-03
## Tensor completion, integral flux, a candidate discrete remnant, and the Higgs parent

**Date:** 17 September 2026  
**Basis:** supplied MTFT 0.31.4 candidate; continuation of AXG-01 and AXG-02  
**Status:** exact conditional calculations and analytic obstructions; no completed parent theory

The proposed combination of a native axion and a six-dimensional tensor does not yet give a viable parent for the specified MTFT matter assignment. Four calculations identify what must change: a residual gauge-anomaly signature obstruction, a factor of three in the integral tensor charge, a nonzero fermionic anomaly for the resulting candidate discrete symmetry, and an elementary-Higgs mass and chirality obstruction. These constrain a precisely specified construction. They do not exclude every parent that could use the modular curve.

An isolated matter addition repairs the discrete fermion test, demonstrating that this particular obstruction is not absolute. It leaves the tensor signature obstruction unchanged and introduces twelve net new four-dimensional fermion modes. No Higgs vacuum, physical particle mass, curvature radius, or Einstein action is derived.

| Test | Exact result | Interpretation |
| --- | --- | --- |
| Residual tensor factorization | A restricted gauge Gram matrix has inertia \((2,1)\) | One ordinary nonchiral tensor product is insufficient; the specified conventional \((1,T)\) pairing cannot work through neutral additions alone |
| Integral tensor flux reduction | An explicit integral coupling induces \(3k_1\), not \(k_1\) | Primitive charge cannot be recovered by silently changing an axion's period |
| Candidate discrete symmetry | Smith invariants \((1,3)\); fermion cubic trace \(24\not\equiv0\pmod9\) | A classical \(\mathbb Z_3\) component requires additional quantum anomaly cancellation |
| Minimal elementary Higgs | \(\lambda_{\min}\geq1/(4R^2)>0\) | Minimal internal flux alone produces no Higgs instability |
| Parent scalar Yukawas | Same-chirality 6D Weyl pairs have no Lorentz-scalar bilinear | Four-dimensional gauge-invariant triangles do not establish six-dimensional scalar interactions |

The accompanying four scripts pass **149 exact checks**, with zero failures. These are algebraic and analytic-consequence checks under the declared assumptions; they are not experimental evidence for MTFT.

## 1. Fixed inputs and conventions

We work on a smooth, closed, compact genus-13 curve \(X=X_0(143)\) with its compact uniformizing metric, curvature \(-1/R^2\), and area

\[
\mathcal A=4\pi(g-1)R^2=48\pi R^2.
\]

This is the compact curve's metric. It must not be confused with the complete finite-area cusped modular metric on the punctured curve. The radius remains a free input.

The gauge group and central flux data are

\[
G=U(3)_c\times U(2)_L\times U(1)_a\times U(1)_b\times U(1)_d,
\qquad N=(3,2,1,1,1),
\]
\[
m=(0,-3,3,3,0)^T,\qquad f_i=F_i/(2\pi).
\]

Here \(m_i\) is the degree of each stack's central line bundle. For a central charge vector \(q\), the associated line degree is \(q\cdot m\). Character charges are normalized so a fundamental has central charge one; a determinant of \(U(N_i)\) has charge \(N_i\).

The continuous directions retained by the axion rows are

\[
p=(1,1,1,1,1)^T,\quad
Y=(1/6,0,-1/2,1/2,-1/2)^T,\quad
B-L=(1/3,0,0,0,-1)^T.
\]

The earlier rows and the native six-dimensional candidate are

\[
k_1=(0,-2,1,1,0),\quad k_2=(3,-4,0,0,1),
\]
\[
k_\Delta=k_2-k_1=(3,-2,-1,-1,1).
\]

Direct multiplication gives \(k_1m=k_2m=12\), while \(k_\Delta m=0\). A smooth fixed-modulus charged scalar on the closed curve requires zero associated degree. Thus \(k_\Delta\) passes this necessary native-scalar condition, whereas \(k_1\) and \(k_2\) do not. A tensor-derived axion obeys different patching rules and need not satisfy this scalar condition.

**The fermion parent being tested is deliberately specific.** It contains one complex 6D Weyl fermion in each unordered bifundamental pair. All have the same chirality except the \((c,\bar d)\) sector, which is reversed. This is the AXG-02 zero-index flip: \(m_c-m_d=0\), so the net four-dimensional chiral index is unchanged. It does not preserve an assumed 6D supersymmetric multiplet assignment. No charged gauginos, gravitinos, tensor multiplets, or ultraviolet completion are silently included.

## 2. Reconstructing the anomaly polynomial by hand

Write \(C=c_{2,c}\), \(L=c_{2,L}\), and \(c_{3,c}\) for the local traceless color and weak characteristic forms. Set

\[
W=\frac{7p_1^2-4p_2}{5760}.
\]

We define positive Weyl chirality to contribute \([\widehat A\,\mathrm{ch}]_8\), with

\[
\widehat A=1-\frac{p_1}{24}+W+\cdots.
\]

The needed traceless Chern characters through degree eight are

\[
\mathrm{ch}_c=3-C+\frac12c_{3,c}+\frac1{12}C^2,
\qquad
\mathrm{ch}_L=2-L+\frac1{12}L^2,
\qquad
\mathrm{ch}_a=\mathrm{ch}_b=\mathrm{ch}_d=1.
\]

For every \(i<j\), multiply

\[
\mathrm{ch}_{i\bar j}=e^{f_i-f_j}\mathrm{ch}_i\mathrm{ch}_{\bar j},
\qquad
I_8^{i\bar j}=[\mathrm{ch}_{i\bar j}]_8-\frac{p_1}{24}[\mathrm{ch}_{i\bar j}]_4+N_iN_jW.
\]

Dualization changes the sign of odd Chern-character degrees. Sum the ten sector contributions, placing a minus sign on \((c,\bar d)\). This recipe reconstructs every term; the script independently rebuilds it and compares it with the archived AXG-02 polynomial.

The net signed representation dimension is \(24-2(3)=18\). Consequently the pure gravitational part is \(18W\). The mixed color-cubic term is

\[
\frac12(k_\Delta f)c_{3,c}.
\]

A native scalar with shift \(k_\Delta\) can only cancel anomaly-polynomial terms divisible by \(k_\Delta f\). A decisive necessary test therefore sets that factor to zero:

\[
(f_c,f_L,f_a,f_b,f_d)=(x,0,y,z,y+z-3x).
\]

Setting \(f_L=0\) removes the matter-invisible common phase. Every scalar-cancellable term vanishes on this restriction. Any remaining anomaly must be handled by other specified fields or interactions.

## 3. A three-square certificate for the tensor obstruction

Restrict further to

\[
(f_c,f_L,f_a,f_b,f_d)=(x,0,0,0,-3x),
\qquad p_1=p_2=0,\qquad X=x^2.
\]

The remaining gauge polynomial is exactly

\[
\boxed{
Q_8=\frac14(C+2L+12X)^2-\frac12(L+18X)^2+108X^2.
}
\]

This is a short hand-checkable certificate. Expanding gives

\[
Q_8=\frac14C^2+CL+\frac12L^2+6CX-6LX-18X^2.
\]

In the independent four-form basis \((C,L,X)\), its matrix is

\[
G=\begin{pmatrix}
1/4&1/2&3\\
1/2&1/2&-3\\
3&-3&-18
\end{pmatrix},\qquad
\det G=-27/2.
\]

The triangular substitution

\[
\begin{pmatrix}u\\v\\w\end{pmatrix}
=\begin{pmatrix}1&2&12\\0&1&18\\0&0&1\end{pmatrix}
\begin{pmatrix}C\\L\\X\end{pmatrix}
\]

has determinant one, and \(Q_8=u^2/4-v^2/2+108w^2\). Sylvester's law of inertia therefore gives **two positive directions and one negative direction**. This is a statement about a coefficient pairing on invariant four-forms, not about spacetime signature or particle masses.

An ordinary single nonchiral tensor supplies one electric-magnetic product \(X_4Y_4\). Its symmetric Gram matrix has rank at most two and at most one direction of either sign. It cannot reproduce this restriction. More generally, a factorization

\[
I_8=+c\,\Omega_{\alpha\beta}X_4^\alpha X_4^\beta,\qquad c>0,\quad \mathrm{signature}(\Omega)=(1,T)
\]

cannot have a restriction with two positive directions. For this unchanged charged spectrum, neutral additions cannot fix the gauge block. The displayed sign is part of the statement: reversing the convention for all anomalies also reverses the matching factorization sign. The conventional 6D supergravity pairing and its separate kinetic-positivity requirements are discussed by [Park and Taylor, §§2.2–2.4 and 2.7](https://arxiv.org/pdf/1110.5916).

This does **not** exclude arbitrary nonsupersymmetric tensor sectors, additional charged matter, a different gauge parent, or a relative theory with specified inflow. An indefinite anomaly pairing alone is not a proof of negative physical kinetic energy. What is excluded is the small, explicitly specified cancellation ansatz.

### Tensor fields also contribute their own anomaly

One signed real chiral two-form contributes

\[
S=\frac{7p_2-p_1^2}{360}=\frac18L_8.
\]

Its opposite chirality reverses the sign. This normalization concerns one real chiral field; local cancellation is distinct from its global consistency. See [Monnier, §2.2, equation (2.2)](https://arxiv.org/pdf/1110.4639).

If \(n\) is the signed number of added neutral complex Weyls and \(t\) the number of tensors contributing \(+S\) minus those contributing \(-S\), the pure gravitational polynomial is \((18+n)W+tS\). Cancellation of irreducible \(p_2\) requires

\[
18+n=28t.
\]

The remaining coefficient of \(p_1^2\) is then \(d=t/32\). One cannot first cancel gravity with neutral fermions and subsequently add anomalous chiral tensors without updating this equation.

On the same Abelian restriction, the full matrix in basis \((C,L,p_1,X)\) is

\[
M(d)=\begin{pmatrix}
1/4&1/2&1/16&3\\
1/2&1/2&1/8&-3\\
1/16&1/8&d&0\\
3&-3&0&-18
\end{pmatrix}.
\]

Eliminating the first two variables leaves

\[
\begin{pmatrix}d-1/64&-3/4\\-3/4&108\end{pmatrix},
\qquad \det=108d-9/4.
\]

| Value of \(d\) | Inertia: positive, negative, zero |
| --- | --- |
| \(d<1/48\) | \((2,2,0)\) |
| \(d=1/48\) | \((2,1,1)\) |
| \(d>1/48\) | \((3,1,0)\) |

The purely illustrative addition of 18 opposite-chirality neutral Weyls gives \(d=0\), with exact diagonal pivots \((1/4,-1/2,-1/64,144)\). None of these choices removes the two-positive-direction obstruction. The supplied JSON also contains a full rational four-form Gram presentation. Its inertia is not used as an invariant: multivariable quartic monomials have relations that make such presentations nonunique.

## 4. Integral flux changes the axion charge lattice

For an ordinary tensor with normalized surface holonomy, write its Abelian characteristic class as

\[
X_4^{\rm ab}=\frac12\sum_{i,j}b_{ij}f_if_j,\qquad b_{ij}=b_{ji}.
\]

An axion obtained by integration over \(X\) has period \(2\pi\). Flux transgression gives

\[
q_i=\sum_jb_{ij}m_j,\qquad \delta\theta=-q_i\Lambda_i.
\]

The factor is \(bm\), not \(bm/2\): expanding the quadratic class supplies two equal external-internal cross terms. A patch-consistent reduction must retain both. An explicit treatment of this issue and of electric/dual axions appears in [Buchmüller et al., equations (55)–(60), (73)–(76)](https://arxiv.org/pdf/1506.05771).

Here

\[
m=3m_0,\qquad m_0=(0,-1,1,1,0)^T.
\]

The vector \(m_0\) is an allowed cocharacter of the actual unitary product group. If the tensor pairing takes integral values on this cocharacter lattice into the primitive tensor charge lattice, every induced charge is divisible by three. This excludes a primitive \(k_1\) tensor shift within that class. Quantization must be tested on the actual gauge group and string-charge lattice; see [Monnier, Moore and Park, equations (2.27), (3.32)](https://arxiv.org/pdf/1711.04777).

An explicit integral control is

\[
X_4=c_2(E_L)+c_1(L_a)c_1(L_b).
\]

On central curvatures its quadratic part is \(f_L^2+f_af_b\), with

\[
b=\begin{pmatrix}
0&0&0&0&0\\
0&2&0&0&0\\
0&0&0&1&0\\
0&0&1&0&0\\
0&0&0&0&0
\end{pmatrix},\qquad
bm=\begin{pmatrix}0\\-6\\3\\3\\0\end{pmatrix}=3k_1^T.
\]

Dividing this class by three is not an innocuous normalization. On spin \(S^2\times S^2\), take trivial \(E_L\), and unit first Chern classes for \(L_a,L_b\) pulled back from the two spheres. Then \(\int c_1(L_a)c_1(L_b)=1\), so \(X_4/3\) has period \(1/3\). Extra fractional or relative structure would have to be supplied to permit it.

Native mixing cannot restore the original primitive row: \(k_1-nk_\Delta\notin3\mathbb Z^5\) for every integer \(n\). The \(d\)-coordinate first forces \(n=0\pmod3\); the \(a\)-coordinate then remains one modulo three.

## 5. The candidate \(\mathbb Z_3\): selection rules and a failed fermion test

The integral control gives

\[
\widehat K=\begin{pmatrix}
3&-2&-1&-1&1\\
0&-6&3&3&0
\end{pmatrix}.
\]

It has the same continuous kernel as the earlier two-row matrix: \(p,Y,B-L\). Its Smith invariants are \((1,3)\). For a hand check, the greatest common divisor of entries is one, every two-column minor is divisible by three, and the minor using columns \(a,d\) equals \(-3\). Hence the second invariant is three.

This is a classical finite component of the unbroken gauge group, before other condensates. It survives the unitary-center bookkeeping: in determinant coordinates the rows are \((1,-1,-1,-1,1)\) and \((0,-3,3,3,0)\), again with invariants \((1,3)\).

Choose its generator

\[
\Lambda=\frac{2\pi}{3}g_3,\qquad g_3=(0,0,1,0,1)^T,\qquad \widehat K g_3=(0,3)^T.
\]

Both periodic axions return to the same point. A field with charge \(q\) transforms by \(\exp(2\pi i\,q\cdot g_3/3)\).

| Field | Central charge | \(\mathbb Z_3\) charge |
| --- | --- | --- |
| \(Q\) | \(e_c-e_L\) | 0 |
| \(u^c\) | \(e_a-e_c\) | 1 |
| \(d^c\) | \(e_b-e_c\) | 0 |
| \(\ell\) | \(e_d-e_L\) | 1 |
| \(e^c\) | \(e_b-e_d\) | 2 |
| \(\nu^c\) | \(e_a-e_d\) | 0 |
| \(H_u\) | \(e_L-e_a\) | 2 |
| \(H_d\) | \(e_L-e_b\) | 0 |

An operator of charge \(q\) can be dressed by integer powers of the two periodic axions only when \(q=\widehat K^T(n_1,n_2)^T\) with integer \(n_1,n_2\). Fractional axion powers are not single-valued under the declared periods.

| Operator | Required rational dressing \((n_1,n_2)\) | Result |
| --- | --- | --- |
| Four ordinary SM/Dirac-neutrino Yukawa charge triangles | \((0,0)\) | Gauge neutral |
| \(H_uH_d\) | \((0,-1/3)\) | No integer dressing |
| \(QQQ\ell\) | \((1,1/3)\) | No integer dressing |
| \(u^cu^cd^ce^c\) | \((-1,1/3)\) | No integer dressing |
| Cubes of the previous three charged operators | Three times the displayed pair | Integer dressing exists |
| \(\nu^c\nu^c\) | Outside the row span | Still carries \(B-L=2\) |

For example, the first permitted pure holomorphic power of the Higgs product is \(e^{-i\theta_T}(H_uH_d)^3\). This is a charge-selection statement, not a generated potential term or a supersymmetric superpotential. Additional charged condensates could alter the remnant. More fundamentally, the quantum anomaly must first cancel.

### The full chiral fermion ledger

The following includes the additional weak doublets in the bifundamental parent, not merely the usual three SM families. All modes are counted in a common four-dimensional chirality, conjugating representations when the index is negative. Integer representatives of discrete charges are used.

| Sector | Multiplicity | Representation dimension | Charge \(s\) | Contribution to \(\sum s\) and \(\sum s^3\) |
| --- | ---: | ---: | ---: | ---: |
| \(cL\) | 3 | 6 | 0 | 0 |
| \(ca\) | 3 | 3 | 1 | 9 |
| \(cb\) | 3 | 3 | 0 | 0 |
| \(La\) | 6 | 2 | 1 | 12 |
| \(Lb\) | 6 | 2 | 0 | 0 |
| \(Ld\) | 3 | 2 | 1 | 6 |
| \(ad\) | 3 | 1 | 0 | 0 |
| \(bd\) | 3 | 1 | −1 | −3 |

The index-zero \(cd\) and \(ab\) sectors contribute no net anomaly. The sums are

\[
S_1=S_3=24.
\]

For ordinary \(\mathrm{Spin}\times\mathbb Z_3\), the pure fermion anomaly test requires

\[
S_3=0\pmod9,\qquad S_1=0\pmod3.
\]

These are the \(n=3\) specialization of [Hsieh, equation (1.2)](https://arxiv.org/pdf/1808.02881). Here the linear condition passes, but \(S_3=6\pmod9\), corresponding to a nonzero anomaly invariant \(2/3\pmod1\) in that convention. Changing representatives or reversing the generator does not make it vanish.

Therefore the displayed fermion sector alone does not support a quantum-consistent gauged \(\mathbb Z_3\). A fully specified global Green–Schwarz, higher-form, topological, or inflow sector might cancel it; that contribution has not been constructed. The selection rules above cannot yet be promoted to physical proton protection.

### An isolated repair control

Add one positive-chirality 6D Weyl fermion in the genuine one-dimensional character

\[
(\det U_L)^{-1}z_a z_b,\qquad q_{\rm new}=k_1.
\]

It is neutral under \(SU(3)\), \(SU(2)\), \(Y\), \(B-L\), and the common phase. It has degree \(q_{\rm new}\cdot m=12\) and discrete charge one. Thus it supplies twelve net modes, giving \(S_1=S_3=36\), which pass this pure discrete fermion test.

This is an added model choice, not a prediction or a minimal complete repair. With \(r_1=k_1f\), its local anomaly changes by

\[
\Delta I_8=\frac{r_1^4}{24}-\frac{p_1r_1^2}{48}+W,\qquad
\Delta I_6=2r_1^3-\frac12p_1r_1.
\]

The net gravitational dimension becomes 19; a combined calculation must use \(19+n=28t\), rather than the earlier 18. After cancelling \(p_2\), \(d=t/32\) is unchanged. On the gauge-inertia witness, \(r_1(x,0,0,0,-3x)=0\). This field therefore leaves that obstruction intact.

## 6. Elementary Higgs fields: the operator determines the mass

Suppose elementary six-dimensional scalar doublets carry the charges of \(H_u,H_d\). Each associated central line has degree \(-6\). For positive kinetic energy and minimal coupling, their internal operator is the connection Laplacian \(D^\dagger D\).

On a closed compact curve,

\[
\langle\phi,D^\dagger D\phi\rangle=\|D\phi\|^2\geq0.
\]

A zero mode would be a parallel section. A nonzero parallel section of a line bundle never vanishes and forces its curvature to vanish, contradicting nonzero degree. Thus there is no scalar zero mode even before a numerical eigenvalue calculation.

For the central HYM connection, write

\[
D^2=-iB\,d\mathrm{vol},\qquad B=\frac{2\pi\deg\mathcal L}{\mathcal A}=-\frac1{4R^2}.
\]

The Kähler identities on scalar sections give

\[
D^\dagger D=2\bar\partial^\dagger\bar\partial+B
=2\partial^\dagger\partial-B.
\]

Using the second expression proves

\[
\boxed{\lambda_{\min}\geq\frac1{4R^2}>0.}
\]

The relevant operator identity is a specialization of the Bochner–Kodaira framework in [Demailly, Chapter VII, §1](https://www-fourier.univ-grenoble-alpes.fr/~demailly/manuscripts/agbook.pdf). Equality requires an appropriate holomorphic section of the dual line, and is not guaranteed by degree alone. In untwisted M1 the effective dual divisor supplies such a section; a flat twist can change that conclusion.

The eigenvalue above is in the internal/product metric. If the four-dimensional Einstein-frame ansatz is

\[
ds_6^2=\rho^{-2}ds_{4,E}^2+\rho^2ds_{X,0}^2,\qquad R=\rho R_0,
\]

then the scalar's four-dimensional kinetic normalization has no radius power, whereas its mass term gives

\[
m_E^2=\rho^{-2}\lambda(R)\geq\frac1{4R_0^2\rho^4}.
\]

Neither \(R_0\) nor \(\rho\) is determined here; this is no prediction in GeV.

### Why eighteen cohomology classes do not contradict this result

For a degree \(-6\) line bundle \(\mathcal L\) on genus 13,

\[
h^0(\mathcal L)=0,\qquad h^0-h^1=-6+1-13=-18,\qquad h^1=18.
\]

Serre duality identifies

\[
H^1(\mathcal L)\simeq H^0(K_X\otimes\mathcal L^{-1})^\vee,\qquad
\deg(K_X\otimes\mathcal L^{-1})=24+6=30.
\]

These are Dolbeault-harmonic forms or dual holomorphic sections. They are not eighteen zero eigenvectors of the minimally coupled scalar connection Laplacian on \(\mathcal L\). Even the degree-30 scalar Bochner operator has holomorphic level \(5/(4R^2)\), rather than zero. A physical twisted/form/gauge-field operator may cancel that shift, but it must be derived from the action.

Internal vector components can have gyromagnetic curvature terms absent from an elementary scalar. The distinction is explicit in the toroidal calculation of [Cremades, Ibáñez and Marchesano, equations (3.65)–(3.70)](https://arxiv.org/pdf/hep-th/0404229). That example does not itself establish the genus-13 spectrum. The present result preserves the cohomology calculation while testing a particular physical interpretation of it.

As a control, adding \(\mu^2+\xi R_6\) to the scalar operator in a flat-external unwarped product gives

\[
\lambda_{\min}^{\rm shifted}\geq\mu^2+\frac{1/4-2\xi}{R^2}.
\]

For a saturated lowest level and \(\mu^2=0\), \(\xi=1/8\) makes it zero. This shows which extra input can change the result; it does not derive that coefficient or a stable vacuum.

## 7. The six-dimensional Yukawa chirality test

The four familiar Yukawa charge sums vanish. A Lorentz contraction must also exist. Let

\[
\Gamma_7=\mathrm{diag}(I_4,-I_4),\qquad P_\pm=(I_8\pm\Gamma_7)/2.
\]

For a scalar Dirac bilinear,

\[
\overline{P_+\Psi}\,P_+\Chi=\bar\Psi P_-P_+\Chi=0.
\]

The representation-theoretic version uses the complexified spin group \(\mathrm{Spin}(6,\mathbb C)\simeq SL(4,\mathbb C)\):

\[
\mathbf4\otimes\mathbf4=\mathbf6\oplus\mathbf{10},\qquad
\mathbf4\otimes\mathbf4^*=\mathbf1\oplus\mathbf{15}.
\]

A scalar occurs for opposite chiralities, not two equal chiralities. The code also checks charge-conjugated bilinears, rather than presuming the Dirac contraction is the only possibility.

For hand replication, start with the Pauli matrices

\[
\sigma_1=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
\sigma_2=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
\sigma_3=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\]

Use the following 8-by-8 matrices, expressed compactly through Kronecker products:

\[
\begin{aligned}
\gamma_1&=\sigma_1\otimes I_2\otimes I_2,&
\gamma_2&=\sigma_2\otimes I_2\otimes I_2,\\
\gamma_3&=\sigma_3\otimes\sigma_1\otimes I_2,&
\gamma_4&=\sigma_3\otimes\sigma_2\otimes I_2,\\
\gamma_5&=\sigma_3\otimes\sigma_3\otimes\sigma_1,&
\gamma_6&=\sigma_3\otimes\sigma_3\otimes\sigma_2.
\end{aligned}
\]

They obey \(\{\gamma_i,\gamma_j\}=2\delta_{ij}I_8\). Define \(C_E=\gamma_2\gamma_4\gamma_6\), \(\Gamma_{7,E}=i\gamma_1\cdots\gamma_6\), and \(P_{+,E}=(I_8+\Gamma_{7,E})/2\). Direct multiplication gives \(P_{+,E}^TC_EP_{+,E}=0\). Complexification suffices for this Lorentz-invariant representation test.

All four proposed SM Yukawa triangles retain the original equal 6D chirality in the tested parent: only \(cd\) was flipped, and it appears in none of those triangles. Ordinary local nonderivative couplings to added elementary Higgs scalars therefore do not implement them.

A vector insertion \(\Gamma^m\) can connect the projectors and evade this scalar obstruction. However, the adjoint of the given direct-product gauge group does not contain the charged off-diagonal Higgs representations. A larger parent gauge group, new opposite-chirality matter, or a specified higher-derivative interaction would change the problem. None has been derived here. Likewise, the isolated chirality flip cannot be called an unchanged 6D \((1,0)\) hypermultiplet construction; the multiplet chiralities are fixed, as tabulated by [Park and Taylor, Table 1](https://arxiv.org/pdf/1110.5916).

## 8. What remains between this construction and gravity

A tensor field on the curve carries more data than one axion. Before flux interactions, ordinary harmonic reduction gives

| Six-dimensional field | Independent 4D real scalars | 4D vectors |
| --- | ---: | ---: |
| One nonchiral two-form | 2 | 26 |
| One chiral two-form | 1 | 13 |

The counts follow from \(b_0=b_2=1\), \(b_1=2g=26\), dualizing the external two-form, and imposing self-duality where appropriate. Flux mixing, constraints, projections and mass terms must be derived before these are called an interacting massless spectrum. Both electric and magnetic integral transgressions retain the factor-of-three issue under the declared lattice assumptions.

Anomaly-cancelling couplings involving curvature do not supply an Einstein kinetic term by themselves. If a six-dimensional Einstein action is independently specified,

\[
S_{\rm EH}=\frac{M_6^4}{2}\int d^6x\sqrt{-g_6}\,R_6,
\]

a fixed unwarped product gives

\[
M_4^2=M_6^4\mathcal A=48\pi M_6^4R^2.
\]

This follows by integrating the coefficient of \(R_4\) over \(X\). It relates two undetermined scales; it does not select either one. The extra Abelian factors and their axions can participate in interactions with a gravitational sector, but the calculations here do not generate its dynamical metric.

A productive next parent must specify, together, its charged fermions and chiralities, physical Higgs operator and Yukawa contraction, tensor lattice and global counterterms, and gravitational action. The unchanged charged spectrum with neutral repairs alone already fails the displayed conventional tensor condition. The smallest next search should therefore vary the charged parent structure explicitly and recompute all these necessary conditions before interpreting a scalar potential as an MTFT vacuum.

## 9. Reproduction and glossary

Run `python run_all.py` in an environment with the pinned requirements. The four audits report 46, 26, 33 and 44 checks, respectively. Source inputs, generated JSON, full matrices, the archived polynomial, and a SHA-256 manifest accompany this report. No original package file or previous investigation is changed.

| Term | Meaning in this report |
| --- | --- |
| Parent theory | A higher-dimensional action and field content whose reduction produces the proposed four-dimensional model |
| Native axion | A periodic scalar already present in six dimensions, with a declared gauge shift |
| Tensor axion | A periodic four-dimensional field obtained by integrating or dualizing a six-dimensional two-form |
| Anomaly polynomial | A characteristic-form expression encoding perturbative failure of a quantum symmetry |
| Green–Schwarz cancellation | Cancellation by specified transforming higher-form or scalar fields and their couplings |
| Inertia | The numbers of positive, negative and zero eigenvalues of a real symmetric form; preserved by invertible congruence |
| Cocharacter | A homomorphism from a circle into the gauge group; determines an allowed integral flux direction |
| Smith invariants | Integer diagonal factors of an integer matrix under invertible integer row/column operations |
| Transgression | Here, integration of a higher-dimensional characteristic class against internal flux to produce lower-dimensional charges |
| HYM connection | A Hermitian Yang–Mills connection; on a line over the curve it has constant central curvature in the chosen metric |
| Bochner Laplacian | The connection operator \(D^\dagger D\), whose zero modes are parallel sections |
| Dolbeault cohomology | Holomorphic-bundle cohomology represented by harmonic forms for the Dolbeault operator; not automatically a physical scalar spectrum |
| Index | Net chirality, which can remain unchanged while vectorlike pairs or massive modes change |
| Control | An explicitly added illustrative choice that tests an obstruction; not a derived prediction |

**Provenance.** The supplied source archive has SHA-256 `2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5`. The copied `smflux.py` source has SHA-256 `6858daddd8b8bf6e4b20dbf1572d69a08a6bbe65cf8d870d9dc31a980e1b2f9e`. Mathematical results in this report are either explicitly derived above or identified as standard framework inputs through the linked primary sources. Exact arithmetic establishes the stated conditional identities; it does not establish an anomaly-free ultraviolet completion or a theory of gravity.
