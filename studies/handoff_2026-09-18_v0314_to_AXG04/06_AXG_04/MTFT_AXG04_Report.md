# MTFT AXG-04
## An explicit charged six-dimensional parent candidate and its limits

**Date:** 18 September 2026  
**Continuation:** AXG-03; supplied MTFT 0.31.4 arithmetic data  
**Main candidate:** C3X — an investigation label, not an MTFT package release  
**Status:** explicit classical EFT, exact local anomaly cancellation and integral tensor data; global GS construction and a realistic vacuum remain open

We construct a concrete six-dimensional action with chiral charged matter, elementary Higgs Yukawas, a nonchiral Green–Schwarz tensor, and an Einstein term. It has an explicit flux background on the closed curve \(X_0(143)\) with exactly three chiral fermion zero modes in each Standard Model sector. Its complete local anomaly polynomial factors into integral characteristic classes.

This is a substantial change to the parent choice. Three families now come from **three six-dimensional copies times one internal mode**, rather than one parent field times three internal modes. The flavor matrices, gauge couplings, Higgs potential and gravitational scale are supplied inputs. The minimal unwarped gravitational solution is AdS, with no separation between the external and internal curvature radii. This is a testable EFT construction, not a derivation of our universe or a claim of unification.

Six accompanying audits reproduce the stated exact algebra, index consequences, matrix ranks and background equations. The final check counts are in the reproduction section.

## 1. What changes, and what is retained

| Ingredient | C3X choice | Status |
| --- | --- | --- |
| Internal complex curve | Closed \(X_0(143)\), genus 13 | Retained arithmetic background |
| Spin structure | Audited \(S_0^2=K_X\), \(h^0(S_0)=2\) | Retained input |
| Gauge group | \(SU(3)\times SU(2)\times U(1)_h\times U(1)_X\), as a direct product | New explicit parent choice |
| Hypercharge normalization | Integer charge \(q_h=6Y\) | Convention |
| Fermions | Three copies of six chiral 6D representations | New multiplicity input |
| Internal flux | One unit in \(U(1)_X\); no hypercharge/non-Abelian flux | New choice |
| Flux line | \(\mathcal L_X=\mathcal O(P+Q-R)\) | Explicit degree-one choice using CM points |
| Higgs | One elementary doublet, neutral under \(X\) | New physical field/operator choice |
| Tensor sector | One ordinary nonchiral two-form; hyperbolic charge lattice | Explicit local GS candidate |
| Gravity | Six-dimensional Einstein–Hilbert term | Supplied action, not derived |
| Flavor | Four arbitrary parent \(3\times3\) Yukawa matrices | Supplied parameters |

The previous independent \(B-L\) and common-phase gauge directions are not retained as additional gauge fields. Nor do we silently retain the old all-bifundamental fermion spectrum or its three-mode Yukawa tensors.

The choice of **direct-product gauge global form** matters. No quotient by a shared Standard Model center is assumed. Changing that global form changes allowed bundles and requires another charge-quantization test.

## 2. The charged spectrum

All entries in the following table occur **three times**. The signs are six-dimensional Weyl chiralities. \(U,D,E,N\) are the physical representations associated with right-handed four-dimensional fields, rather than their left-handed conjugates.

| Field | \(SU(3)\) | \(SU(2)\) | \(q_h=6Y\) | \(q_X\) | 6D chirality |
| --- | --- | --- | ---: | ---: | ---: |
| \(Q\) | \(\mathbf3\) | \(\mathbf2\) | 1 | 1 | \(+\) |
| \(U\) | \(\mathbf3\) | \(\mathbf1\) | 4 | 1 | \(-\) |
| \(D\) | \(\mathbf3\) | \(\mathbf1\) | −2 | 1 | \(-\) |
| \(\ell\) | \(\mathbf1\) | \(\mathbf2\) | −3 | 1 | \(+\) |
| \(E\) | \(\mathbf1\) | \(\mathbf1\) | −6 | 1 | \(-\) |
| \(N\) | \(\mathbf1\) | \(\mathbf1\) | 0 | 1 | \(-\) |
| \(H\), scalar | \(\mathbf1\) | \(\mathbf2\) | 3 | 0 | — |

There is one Higgs field, not three copies. Define
\[
\widetilde H=i\sigma_2H^*,\qquad
\sigma_2=\begin{pmatrix}0&-i\\i&0\end{pmatrix}.
\]
It has charges \((-3,0)\).

The fermion chirality assignment is a known nonsupersymmetric 6D Standard Model pattern. The distinction between its six-dimensional multiplicity and four-dimensional compactification modes is essential; see [Dobrescu–Poppitz, equations (2)–(4)](https://arxiv.org/pdf/hep-ph/0102010). We are testing a particular arithmetic compactification and tensor realization of that pattern, not claiming to invent the assignment.

Each copy has signed representation dimension
\[
6+2-3-3-1-1=0.
\]
There are twelve positive-chirality weak doublets when color and the three parent copies are counted, and zero net signed color fundamentals. These counts satisfy the familiar necessary non-Abelian six-dimensional conditions.

### Connection with the earlier five-stack notation

The formal curvature embedding is
\[
\begin{pmatrix}f_c\\f_L\\f_a\\f_b\\f_d\end{pmatrix}
=
\begin{pmatrix}
1&0\\
0&-1\\
-3&-1\\
3&-1\\
-3&0
\end{pmatrix}
\begin{pmatrix}h\\x\end{pmatrix},
\qquad h=\frac{F_h}{2\pi},\quad x=\frac{F_X}{2\pi}.
\]
The first column is \(6Y\); the second is the unit-flux direction \((0,-1,-1,-1,0)^T\). It annihilates \(k_1=(0,-2,1,1,0)\) identically. Thus the earlier native axion is unnecessary in this subgroup construction. This is a genuine homomorphism into the unitary product, while the allowed backgrounds are now those of the chosen direct-product covering group.

## 3. An explicit line bundle giving exactly one mode per copy

An index alone would not suffice: degree one guarantees \(h^0-h^1=1\), and can still leave vectorlike zero-mode pairs.

The supplied arithmetic-spin construction gives
\[
S_0=\mathcal O\!\left(6(0)+6(1/11)\right),\qquad
S_0^2\simeq K_X,\qquad
H^0(S_0)=\operatorname{span}\{1,u\},
\]
where
\[
u(\tau)=\frac{\eta(13\tau)\eta(143\tau)}
{\eta(\tau)\eta(11\tau)}.
\]
These are inherited arithmetic inputs, recorded in the bundled v0.31.4 source. Their previous proofs are not being replaced by the new tests.

Choose three distinct \(W_{13}\)-fixed CM points \(P,Q,R\), with
\[
u(P)=+\frac{i}{\sqrt{13}},\qquad
u(Q)=-\frac{i}{\sqrt{13}},
\]
and set
\[
\boxed{\mathcal L_X=\mathcal O(P+Q-R).}
\]
There are four such fixed points, with two of each sign, so this choice is available. The negative divisor coefficient is a presentation of a smooth holomorphic line bundle; it is not a singular magnetic source.

Here is the hand proof.

Evaluation of \(1,u\) at \(R\) has rank one. Therefore
\[
h^0(S_0(-R))=1.
\]
Serre duality and Riemann–Roch give
\[
h^0(S_0(R))-h^0(S_0(-R))=1,
\]
hence \(h^0(S_0(R))=2\). The natural inclusion \(H^0(S_0)\to H^0(S_0(R))\) is an isomorphism.

Now impose vanishing at \(P,Q\). Because neither equals \(R\), the evaluation matrix, up to nonzero row scalings, is
\[
E=\begin{pmatrix}
1&i/\sqrt{13}\\
1&-i/\sqrt{13}
\end{pmatrix},
\qquad \det E=-\frac{2i}{\sqrt{13}}\ne0.
\]
Thus
\[
h^0(S_0(R-P-Q))=0.
\]
Since \(K_XS_0^{-1}\simeq S_0\), this is \(h^1(S_0\mathcal L_X)=0\). Finally
\[
\deg(S_0\mathcal L_X)=13,\qquad
h^0-h^1=13+1-13=1,
\]
so
\[
\boxed{(h^0,h^1)(S_0\mathcal L_X)=(1,0).}
\]

All six physical 6D representations have \(X=1\), and therefore share this internal line. Their opposite 6D chiralities give the desired opposite 4D chiralities. Before electroweak breaking, three parent copies yield exactly three chiral zero modes per representation, with no opposite zero-mode partners from these fermion sectors.

For comparison, \(\mathcal L_X=\mathcal O(P)\) gives \((h^0,h^1)=(2,1)\) at these points. That superficially simpler unit-flux choice would introduce a vectorlike pair per copy. The subtraction in \(P+Q-R\) has a concrete purpose.

A degree-one holomorphic line admits a smooth constant-curvature unitary connection. No explicit pointwise metric or eigenfunction normalization is required for this cohomological count; those would be needed for detailed nonconstant overlap computations.

## 4. The six-dimensional action

Use mostly-plus Lorentzian signature and positive gauge/tensor kinetic coefficients. Let \(\Psi_{f,r}\) denote each charged Weyl field, with parent-copy label \(r=1,2,3\). A local action is
\[
\begin{aligned}
S_{\mathrm{bulk}}=\int d^6x\sqrt{-g}\,\Bigg[
&\frac{M_6^4}{2}R_6-\Lambda_6
-\sum_{A=3,2}\frac{1}{2g_{6A}^2}
\operatorname{tr}(F_{A\,MN}F_A^{MN})
-\sum_{A=h,X}\frac{F_{A\,MN}F_A^{MN}}{4g_{6A}^2}\\
&-(D_MH)^\dagger D^MH-V_6(H)
+\sum_{f,r}i\bar\Psi_{f,r}\Gamma^MD_M\Psi_{f,r}
-\big(\mathcal Y_6+\mathrm{h.c.}\big)\Bigg],
\end{aligned}
\]
where the non-Abelian fundamental trace satisfies \(\operatorname{tr}(T^aT^b)=\delta^{ab}/2\), and
\[
V_6(H)=\mu_6^2H^\dagger H+\lambda_6(H^\dagger H)^2,\qquad \lambda_6>0,
\]
\[
\begin{aligned}
\mathcal Y_6={}&
\bar Q_r\,\widetilde H\,(\mathbf y_{6u})_{rs}U_s
+\bar Q_r\,H\,(\mathbf y_{6d})_{rs}D_s\\
&+\bar\ell_r\,H\,(\mathbf y_{6e})_{rs}E_s
+\bar\ell_r\,\widetilde H\,(\mathbf y_{6\nu})_{rs}N_s.
\end{aligned}
\]
Repeated parent labels are summed. Gauge indices contract in the usual fundamental/antifundamental pairs.

Every term is gauge invariant. For example, the up vertex has \(h\)-charge \(-1-3+4=0\) and \(X\)-charge \(-1+0+1=0\). Each scalar bilinear joins opposite 6D chiralities, removing the AXG-03 projector obstruction.

The engineering dimensions are
\[
[H]=2,\quad[\Psi]=5/2,\quad[g_6]=[\mathbf y_6]=-1,\quad[\lambda_6]=-2.
\]
This is a nonrenormalizable effective field theory. Its coefficients and cutoff require physical input; an ultraviolet completion has not been supplied.

Add a normalized nonchiral two-form \(B\), whose surface holonomy is \(\exp(2\pi i\int B)\):
\[
S_B=-\frac{1}{2g_B^2}\int\mathcal H_3\wedge *\mathcal H_3
+2\pi\int B\wedge\mathcal Y_4,
\]
\[
\mathcal H_3=dB+\omega_3(\mathcal X_4),\qquad
d\omega_3(\mathcal X_4)=\mathcal X_4.
\]
The next section supplies both characteristic forms and the exact anomaly cancellation.

This formula specifies local patches and local gauge transformations. When \([\mathcal X_4]\ne0\), it is not a globally smooth ordinary \(B\)-field on an arbitrary bundle. A global differential-cohomology/string-source definition remains part of completing the quantum theory.

## 5. Exact anomaly cancellation and the tensor lattice

Define
\[
C=c_2(SU(3)),\qquad W=c_2(SU(2)),\qquad \eta=\frac{p_1(TM_6)}2.
\]
On a spin manifold \(\eta\) is integral. With positive Weyl contribution \([\widehat A\,\mathrm{ch}]_8\), the entire three-copy matter anomaly is
\[
\boxed{
I_8=(W+9h^2)(3C+W+\eta-27h^2-6x^2).
}
\]
Thus choose
\[
\mathcal X_4=W+9h^2,\qquad
\mathcal Y_4=3C+W+\eta-27h^2-6x^2.
\]
Both are integral characteristic classes for the declared direct-product group on spin backgrounds. There are no residual irreducible color-cubic or pure gravitational terms.

### Hand reconstruction

For each field in the spectrum table, multiply its color, weak and Abelian Chern characters:
\[
\mathrm{ch}_{\mathbf3}=3-C+\frac12c_3+\frac1{12}C^2,\qquad
\mathrm{ch}_{\mathbf2}=2-W+\frac1{12}W^2,
\]
\[
\mathrm{ch}_{\mathrm{Abelian}}=\exp(q_hh+x),\qquad
\widehat A=1-\frac{p_1}{24}+\frac{7p_1^2-4p_2}{5760}+\cdots.
\]
Take the degree-eight part, multiply by the chirality sign and by three copies, then sum. For a quick expansion check:
\[
I_8=
3CW+W^2+\frac12p_1W
+27Ch^2-18Wh^2+\frac92p_1h^2
-243h^4-6Wx^2-54h^2x^2.
\]
Two scripts reconstruct the result through physical representations and the original quiver embedding; an additional Cartan-weight calculation checks the character formula.

The tensor source vector uses the hyperbolic lattice
\[
\Omega=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
I_8=\frac12
\begin{pmatrix}\mathcal X_4&\mathcal Y_4\end{pmatrix}
\Omega
\begin{pmatrix}\mathcal X_4\\\mathcal Y_4\end{pmatrix}.
\]
It is even, unimodular, and has signature \((1,1)\). This is the electric/magnetic charge pairing of a nonchiral tensor, not a negative physical kinetic-energy eigenvalue.

Writing the source convention as
\[
\binom{\mathcal X_4}{\mathcal Y_4}
=a\,\frac{p_1}{4}-b_C C-b_W W+\frac12b_{hh}h^2+\frac12b_{XX}x^2+b_{hX}hx,
\]
the coefficient vectors are
\[
a=\binom02,\quad b_C=\binom0{-3},\quad b_W=\binom{-1}{-1},
\]
\[
b_{hh}=\binom{18}{-54},\quad b_{XX}=\binom0{-12},\quad b_{hX}=\binom00.
\]
The Abelian diagonal vectors are even. For \(v=(v_1,v_2)^T\),
\[
v^T\Omega v=2v_1v_2,\qquad
a^T\Omega v-v^T\Omega v=2v_1(1-v_2),
\]
so \(a\) is characteristic. These explicit integral checks address the obstruction that rational factorization alone misses. The role of actual cocharacters and tensor-charge quantization is explained in [Monnier–Moore–Park](https://arxiv.org/abs/1711.04777); this report does not assume their supersymmetric dynamics.

Under descent, choose \(I_7=\omega_3(\mathcal X_4)\mathcal Y_4\). If
\[
\delta B=-\omega_2^{(1)}(\mathcal X_4),
\]
the variation of \(2\pi\int B\mathcal Y_4\) cancels the fermion descent in this convention. Local anomaly cancellation is therefore explicit.

### A separate global-topology check

For the ordinary spin background category and the chosen direct-product group, the accompanying mod-two calculation finds
\[
\boxed{\Omega_7^{\mathrm{Spin}}(BG)=0.}
\]
A short reproducible proof uses
\[
H^*(BG;\mathbb F_2)=\mathbb F_2[h,x,C,W,D],\qquad
\deg(h,x,C,W,D)=(2,2,4,4,6),
\]
where \(D=c_3(SU(3))\). The Steenrod-square rules are
\[
Sq^2h=h^2,\quad Sq^2x=x^2,\quad Sq^2C=D,\quad Sq^2W=Sq^2D=0.
\]
The sequence \(H^4\to H^6\to H^8\) has dimensions \(5,9,16\), with matrix ranks \(2,7\). Its middle kernel equals the image:
\[
\operatorname{span}_{\mathbb F_2}\{h^2x+hx^2,D\}.
\]
The dual differential removes the only possible total-degree-seven term in the spin Atiyah–Hirzebruch spectral sequence. Full matrices and exhaustive finite-field checks are included.

This leaves no additional torsion bordism anomaly in the ordinary spin-\(BG\) category **once globally defined cancellation is constructed**. It does not construct that cancellation. A category containing a chosen trivialization \(d\mathcal H_3=\mathcal X_4\) can have different bordism. These distinctions and the differential rules are treated in [Lee–Tachikawa, Appendix B and §3.3](https://arxiv.org/pdf/2012.11622). We do not promote this result to a completed GS partition function.

## 6. Flux reduction and the extra Abelian gauge field

Take
\[
\int_X x=1,\qquad \int_Xh=0.
\]
The two characteristic forms transgress as
\[
\int_X\mathcal X_4=0,\qquad
\int_X\mathcal Y_4=-12x_4.
\]
The second factor of two comes from the two cross terms in \((x_{\mathrm{internal}}+x_4)^2\).

The four-dimensional topological coupling is therefore
\[
S_{\mathrm{BF},4}=-24\pi\int_{M_4}B_2\wedge x_4.
\]
Its integer BF level is \(-12\). Dualizing the external two-form gives a periodic axion whose charge has magnitude twelve under \(U(1)_X\). With finite positive kinetic coefficients, this supplies a Stückelberg mass for the \(X\) vector while preserving hypercharge. The mass depends on those coefficients and the radius; no value in GeV follows.

For primitive periods and the chosen covering group, the classical remaining finite component is \(\mathbb Z_{12}\). In common left-handed 4D notation, \(Q,\ell\) have \(X=1\), while \(u^c,d^c,e^c,\nu^c\) have \(X=-1\).

The exact reduced anomaly is
\[
I_6=\frac{\partial I_8}{\partial x}=-12x(W+9h^2).
\]
For the complete three-family fermion ledger,
\[
\sum s=\sum s^3=0,\quad
\sum q_h^2s=-216,\quad
\sum q_hs^2=0.
\]
The weak-instanton \(X\)-charge is twelve and the color-instanton charge is zero. The displayed pure \(\mathrm{Spin}\times\mathbb Z_{12}\) congruences and mixed instanton tests pass; the pure congruences use [Hsieh, equation (1.2)](https://arxiv.org/pdf/1808.02881). This is stronger than the failed \(\mathbb Z_3\) fermion test of AXG-03, while still conditional on completing the full global theory.

| Operator | \(X\)-charge in left-handed notation | Candidate \(\mathbb Z_{12}\) rule |
| --- | ---: | --- |
| Standard Yukawas | 0 | Allowed |
| \(QQQ\ell\) | 4 | Forbidden |
| \(u^cu^cd^ce^c\) | −4 | Forbidden |
| \(\nu^c\nu^c\) | −2 | Forbidden |
| \((QQQ\ell)^3\) | 12 | Can have an integer axion dressing |

These are selection rules, not computed rates or a proton-lifetime prediction. New charged condensates can change the remnant. In particular, the displayed parent naturally has Dirac neutrinos; adding a Majorana sector requires more structure.

## 7. How fermion masses arise in this action

The Higgs is neutral under the only internal flux. Its minimal kinetic operator has a constant internal profile, normalized as
\[
H_6(x,y)=\frac{H_4(x)}{\sqrt{\mathcal A}},
\qquad \mathcal A=48\pi R^2.
\]
Let \(\chi(y)\) be the single normalized internal fermion mode,
\[
\int_X\chi^\dagger\chi\,d\mathrm{vol}=1.
\]
All copies and species share the same \(X\)-line and internal mode; their internal hypercharge and non-Abelian connections are trivial. The Yukawa integral gives
\[
\boxed{\mathbf Y_{4f}=\frac{\mathbf y_{6f}}{\sqrt{\mathcal A}}.}
\]
For a fixed product metric, the potential reduces to
\[
\mu_4^2=\mu_6^2,\qquad
\lambda_{H,4}=\frac{\lambda_6}{\mathcal A}.
\]
If the chosen coefficient \(\mu_6^2\) is negative,
\[
H_4=\frac1{\sqrt2}\binom0{v_4},\qquad
v_4^2=-\frac{\mu_6^2\mathcal A}{\lambda_6},
\]
and
\[
\boxed{\mathbf M_f=\frac{v_4}{\sqrt2}\mathbf Y_{4f}.}
\]
This is an explicit mechanism connecting a 6D scalar interaction, normalized internal modes and 4D fermion masses. It does not derive the tachyonic coefficient, the Yukawa matrices or their singular values. Indeed, substituting the chosen potential minimum cancels the area from the tree-level mass:
\[
\mathbf M_f=
\frac{\mathbf y_{6f}}{\sqrt2}
\sqrt{-\frac{\mu_6^2}{\lambda_6}}.
\]
The arithmetic geometry controls the existence/count of modes here; this simple common-profile construction supplies no flavor hierarchy by itself.

Gauge masses arise through the Higgs covariant derivative and the separate \(X\) Stückelberg mechanism. Their couplings are also input parameters. No common-coupling relation or weak-mixing-angle prediction has been imposed.

In an Einstein-frame radius variation the appropriate canonical rescalings must be included. For \(ds_6^2=\rho^{-2}ds_{4,E}^2+\rho^2ds_{X,0}^2\), the physical mass in that Einstein frame is \(M_{f,E}=\rho^{-1}M_{f,\mathrm{product}}\). At a fixed radius the expressions above are unambiguous; the dimensionless Yukawa scales as \(1/\sqrt{\mathcal A}\). These formulas do not silently reuse the older three-internal-mode Yukawa tensors.

## 8. Does the action support its proposed spacetime?

An Einstein term makes this question dynamical. Use the **smooth closed** genus-13 uniformizing metric, Gaussian curvature \(-1/R^2\). Its area is \(48\pi R^2\). The cusped modular metric has a different area and is not substituted into these equations.

Assume an unwarped maximally symmetric product, a covariantly constant Higgs at a critical point of its potential, internal \(X\)-flux, and zero other background fields. Set
\[
M=M_6^4,\quad
R_{\mu\nu}=\lambda_4g_{\mu\nu},\quad
R_{mn}=kg_{mn},\quad
U=\Lambda_6+V_6(H_0).
\]
For positive magnetic energy density \(\rho_F\),
\[
T_{\mu\nu}=-(U+\rho_F)g_{\mu\nu},\qquad
T_{mn}=(\rho_F-U)g_{mn}.
\]
The six-dimensional Einstein equations give
\[
\lambda_4=\frac{U-\rho_F}{2M},\qquad
k=\frac{U+3\rho_F}{2M}.
\]
Setting \(k=-1/R^2\) yields
\[
U=-\frac{2M}{R^2}-3\rho_F,\qquad
\boxed{\lambda_4=-\frac1{R^2}-\frac{2\rho_F}{M}<0.}
\]

**The stated action and ansatz admit neither a Minkowski nor a de Sitter product on this negatively curved surface.** This is a bounded result; warping, gradients, localized sources and corrections change the equations and must be assessed explicitly. The flux/radion framework can be compared with [Navarro–Santiago, §§2–3](https://arxiv.org/pdf/hep-th/0405173).

For \(R_{\mu\nu}=-3g_{\mu\nu}/\ell_4^2\),
\[
\ell_4^2=\frac{3}{R^{-2}+2\rho_F/M}<3R^2.
\]
There is no parametric external/internal curvature-radius separation in this product. The zero-mode reduction is a well-defined calculation of coefficients, but this background does not establish a controlled near-flat four-dimensional regime.

With unit \(X\)-flux and Maxwell term \(-F_X^2/(4g_{6X}^2)\),
\[
\rho_F=\frac{2\pi^2}{g_{6X}^2\mathcal A^2}
=\frac{\beta}{R^4},\qquad
\beta=\frac1{1152g_{6X}^2}.
\]
If \(U=-u_0<0\), the unique positive radius solution is
\[
R^2=\frac{M+\sqrt{M^2+3u_0\beta}}{u_0}.
\]
The scale is fixed only after supplying action parameters.

For
\[
ds_6^2=e^{-2\sigma}ds_{4,E}^2+e^{2\sigma}ds_{X,0}^2,
\qquad \phi=2M_P\sigma,
\]
the radius potential has the form
\[
V_E=a e^{-4\sigma}+b e^{-6\sigma}+c e^{-2\sigma},\quad a,b>0.
\]
At a stationary point chosen as \(\sigma=0\), \(c=-2a-3b\) and \(V_*=-a-2b<0\). The canonical radion has
\[
\boxed{m_\phi^2=\frac2{R^2}+\frac{6\rho_F}{M}>0.}
\]
This is stability of the radius in the specified truncation, not full KK or quantum stability.

### An explicit classical control

In arbitrary consistent mass units choose
\[
M_6^4=1,\quad g_{6X}^2=\frac1{1152},\quad
\mu_6^2=-1,\quad\lambda_6=1,\quad\Lambda_6=-\frac{19}{4}.
\]
Then
\[
H_0^\dagger H_0=\frac12,\quad V_6(H_0)=-\frac14,\quad U=-5,\quad\beta=1,
\]
\[
R^2=\ell_4^2=1,\qquad m_\phi^2=8.
\]
The radial Higgs mass squared at that reference radius is \(2\). On the maximally symmetric product, with the other gauge backgrounds zero, both GS four-form sources vanish as differential forms: the internal two-form flux squares to zero and the Pontryagin form of the constant-curvature product is zero. A zero tensor background is therefore compatible with its local equations.

This exhibits an action with an on-shell classical product background. It is not a fit to observation, a GeV prediction, a fully stable vacuum, or a weak-coupling/UV-control certificate.

### Remaining geometric fields

The classical curvature, constant potential and central magnetic energy depend on area, not the hyperbolic shape. The \(6g-6=72\) real shape directions are not selected by this background potential. Thus choosing the arithmetic complex structure \(X_0(143)\) is still an input.

An ordinary nonchiral two-form on this curve also has, before interactions and projections, two scalar modes and \(2g=26\) vector modes. One scalar participates in the \(X\) Stückelberg mechanism. The remaining harmonic sector requires a physical spectrum and lifting analysis. It cannot be omitted merely because the charged fermion zero modes match three SM families.

The four-dimensional Planck coefficient is
\[
M_P^2=M_6^4\mathcal A.
\]
It follows from the supplied Einstein action; it is not a derivation of gravity from the extra Abelian factor.

## 9. Why the closer alternatives were insufficient

The calculation retained unsuccessful candidates so that the changes leading to C3X are reviewable.

### One parent and index-three flux

Keep six bifundamental Weyls with signs
\[
(cL,ca,cb,Ld,ad,bd)=(+,-,-,+,-,-).
\]
The desired SM index signs force
\[
m'=(0,-3,-3,-3,0)
\]
up to a common shift. This makes elementary Higgs line degrees zero and permits all four scalar Yukawas.

Modulo the native \(k_1\) scalar direction, define
\[
a=f_c-f_L,\quad b=f_a-f_L,\quad d=f_d-f_L,\quad f_b-f_L=-b.
\]
Its anomaly factors locally:
\[
I_8=(W+b^2)
\left(C+\frac W3+\frac{p_1}{6}-\frac32a^2-\frac12d^2-\frac13b^2\right).
\]
But the genuine primitive cocharacter \((0,0,1,-1,0)\) gives a quartic coefficient \(-1/3\). In the convention \(X_4=\tfrac12 b_{\mathrm{diag}}f^2\) and \(I_8=\tfrac12\Omega(X_4,X_4)\), the required string norm is \(-8/3\), which is not integral. On spin \(S^2\times S^2\), a line with \(c_1=u+v\) has \(\int c_1^2=2\), making the string-charge test explicit.

Adding two index-zero positive weak bifundamentals and four negative neutral Weyls repairs the conventional weak-doublet count but still requires norm \(-4/3\). The extra matter also restores a two-positive-direction gauge pairing obstruction to the small conventional tensor ansatz. These repairs are insufficient.

Furthermore, identifying the line connections to obtain constant Higgs modes makes the three internal family spaces identical. The overlap becomes
\[
Y_f=\frac{y_{6f}}{\sqrt{\mathcal A}}U_f^\dagger U_Q,\qquad
Y_f^\dagger Y_f=\frac{|y_{6f}|^2}{\mathcal A}I_3.
\]
All three singular values in a sector coincide. A CKM matrix can be chosen identity within the degeneracy, but it is then unobservable and basis-undetermined, not a prediction of the measured CKM matrix. C3X escapes this particular tree-level degeneracy by supplying parent-generation matrices; that is a cost in predictivity.

### A larger \(U(8)\) gauge parent

A six-dimensional \(U(8)\) action with one adjoint Dirac fermion supplies off-diagonal charged gauge components and the group factors needed for gauge-origin Yukawa vertices. It also cancels paired fermion anomalies.

However, the two six-dimensional chiralities give opposite index matrices and retain mirrors under the same smooth flux. The original block flux alone has centralizer \(U(4)\times U(2)\times U(2)\); extra holonomies are needed to distinguish all five blocks.

Its unequal-slope split connection is a saddle of the fixed-metric Yang–Mills functional. The complex negative normal dimension is
\[
8(12+3)+4(12+6)+8(12+3)=312,
\]
or real Morse index \(624\), with \(180\) complex colored negative directions. This uses [Atiyah–Bott, Proposition 5.4 and equation (5.10)](https://www.uvm.edu/~cvincen1/files/teaching/spring2019-math382/atiyahbott.pdf). The degree-\(-6\) Higgs cohomology does acquire an interpretation as negative gauge-Hessian directions in this different action. It does not thereby become a selected electroweak vacuum. Coupled AdS stability is a further question, not established by the Morse count alone.

## 10. What has been achieved, and the next bounded problem

C3X supplies an explicit charged spectrum, scalar and tensor interactions, an exact three-family compactification of those fermions, integral local anomaly-cancellation data, a finite-field bordism calculation, and an on-shell classical AdS control. It closes several gaps that were previously only proposals.

It does not derive the parent group, the number of parent copies, flavor matrices, an absolute scale, the arithmetic complex structure, or a realistic spacetime vacuum. Nor does it provide the complete global GS partition function or a UV completion. It is not a supersymmetric construction.

The next work should keep these questions separate:

1. Construct the global GS/differential-cohomology sector for the stated direct-product group, with its tensor strings and allowed bundle backgrounds.
2. Derive a mechanism that changes the AdS/no-scale-separation result while keeping the charged spectrum and anomaly data consistent.
3. Determine the remaining tensor, Wilson-line and shape spectrum, and calculate how any proposed lifting sector changes the action.
4. If arithmetic flavor prediction remains the goal, seek an interaction or parent structure that derives the generation matrices rather than inserting them.

A further arbitrary potential or unlisted spectator cannot be treated as a solution to these questions. Each proposed change must return to the same anomaly, mode-count and field-equation tests.

## 11. Reproduction, provenance and terminology

Run the bundled run_all.py after installing requirements.txt. All calculations use exact rational, symbolic or finite-field arithmetic; no numerical spectrum or observed-mass fit is used.

| Audit | Checks |
| --- | ---: |
| Chiral-parent controls and subgroup factorization | 84 |
| Explicit charged action and reductions | 62 |
| Higgs action, divisor purity and flavor | 68 |
| Einstein background and radion | 47 |
| \(U(8)\) comparison | 31 |
| Spin-bordism calculation | 30 |
| **Total** | **322** |

A passing assertion establishes the stated identity or conditional consequence. The underlying standard index, Kähler and anomaly frameworks remain mathematical assumptions of the construction; check counts are not a measure of physical confirmation.

The supplied MTFT source archive has SHA-256
**2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5**.
The original source files used as arithmetic/comparison inputs, their license, the generated result ledgers, full matrices, supplementary derivations and a file-hash manifest are included. Previous reports and package files are unchanged.

| Term | Meaning here |
| --- | --- |
| Parent copy | A distinct six-dimensional field with its own generation label, before compactification |
| Internal mode | A section solving an internal kinetic equation; several may arise from one parent field |
| Pure chiral zero-mode spectrum | One chirality has zero modes while the opposite kernel vanishes; stronger than a net index |
| Green–Schwarz tensor | A two-form with transforming field strength and topological couplings that cancel gauge/Lorentz anomalies |
| Integral tensor lattice | The lattice of electric/magnetic string charges with integral Dirac pairings |
| Characteristic vector | A lattice vector \(a\) satisfying \(a\cdot v=v\cdot v\pmod2\) for every lattice vector |
| BF coupling | A four-dimensional \(B_2\wedge F_2\) interaction; dual to a periodic-axion Stückelberg coupling |
| Radion | The scalar degree of freedom describing the internal radius |
| AdS product | A spacetime with a negatively curved maximally symmetric four-dimensional factor |
| Bordism check | A global-topology anomaly test in a precisely declared category of background fields |
| EFT | A local effective field theory with a cutoff and supplied couplings, not automatically a UV-complete theory |
