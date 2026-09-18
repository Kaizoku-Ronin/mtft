# MTFT SC7-01: testing the spin-circle route to seven dimensions

**Investigation date:** 16 September 2026  
**Input:** supplied MTFT 0.31.4 candidate  
**Status:** exact topology and zero-mode calculations; conditional classical gravity calculation  
**Outcome:** the proposed geometry exists, but its simplest smooth, commuting, twisted 7D gauge extension does not reproduce the M1 three-family spectrum.

## 1. What was tested

The candidate is

$$
M^{1,3}\times P,\qquad
P=S(S_0)\xrightarrow{\pi}X,\qquad
X=X_0(143),\quad S_0^{\otimes2}\simeq K_X.
$$

Here X is the **compact** genus-13 curve, and S(S0) is the unit-circle bundle of its degree-12 theta characteristic. Thus the internal space has three real dimensions and the gauge theory has seven spacetime dimensions.

The question was whether the existing arithmetic spin geometry and M1 flux bundles can provide the internal space of a standard twisted 7D supersymmetric gauge construction while retaining the observed type of chiral matter.

We fixed the following scope before running the experiment: smooth closed P, commuting stack backgrounds, the standard partial topological twist, and no defects, boundaries, singular sources, or additional projections. The original six-dimensional index was not transplanted into the new theory. The complete preregistration is included.

The standard 7D Higgs-bundle framework uses a gauge connection A and an adjoint one-form phi. In anti-Hermitian conventions its source-free BPS equations can be written

$$
F_A-\phi\wedge\phi=0,\qquad d_A\phi=0,\qquad d_A^*\phi=0.
$$

The complex connection A+i phi is flat. Its de Rham cohomology controls the massless fields; matter modes occur in degrees one and two. These are the equations and mode dictionary adopted here. See [Braun, Cizel, Hübner and Schäfer-Nameki, sections 2.1–2.3](https://arxiv.org/html/1812.06072v2).

**Parent-theory qualification.** Pure SYM with only the product gauge group has adjoints of its separate factors; it does not automatically contain the M1 bifundamentals. A natural way to place the candidate sectors in the standard 7D adjoint complex is a provisional U(8) bundle

$$
\mathcal E=\bigoplus_{x=c,L,a,b,d}
  \bigl(\pi^*L_x\bigr)\otimes\mathbb C^{N_x},
\qquad (N_x)=(3,2,1,1,1).
$$

Its off-diagonal adjoint bundles include the required Hom bundles. This is a possible bookkeeping realization, **not a derived parent gauge group or an established ultraviolet completion**. The calculation below grants the candidate sectors this realization and still finds that their required zero modes disappear. Equivalently, it evaluates each proposed sector conditional on the standard twisted adjoint mode complex. It does not assume a new unspecified 7D matter action.

## 2. Results at a glance

| Question | Result | Meaning |
|---|---|---|
| Does the spin circle define a closed 3-manifold? | Yes; Euler class 12 | The geometric proposal is well defined. |
| What happens to degree-three flux? | It becomes an order-four torsion class on P | It survives topologically, but no longer as a nonzero real cohomology class. |
| Can its pullback support a commuting BPS connection? | Yes, after changing to a flat connection | The naive curved pullback itself fails the flatness equation. |
| What is the resulting fiber holonomy? | Minus i for oriented degree +3 | The torsion is detectable by parallel transport. |
| Do the original three-family sectors have massless modes? | No; all their twisted cohomology groups vanish | The specified 7D uplift fails the matter gate. |
| Do all sectors vanish? | No; degree-zero controls survive | The failure is specific to the required sectors, not an empty or broken calculation. |
| Does the spin-compatible Q8 algebra survive geometrically? | Yes | This does not establish Q8 symmetry of a flux vacuum. |
| Is P a hyperbolic 3-manifold? | No | A hyperbolic base and a circle bundle do not produce constant negative curvature in three dimensions. |
| Does the restricted gravity ansatz fix the overall size? | No | Its classical Einstein-frame scale potential is strictly decreasing. |
| Has a Hecke-equivariant physical parent been obtained? | No | That construction remains open. |

## 3. Integral topology: flux becomes torsion

Write g=13 and e=deg(S0)=g−1=12. The homotopy sequence of the circle bundle, using the asphericity of the genus-13 base, gives an infinite cyclic fiber subgroup. One presentation is

$$
\pi_1(P)=
\left\langle a_1,b_1,\ldots,a_{13},b_{13},t\ \middle|\
[a_j,t]=[b_j,t]=1,\ 
\prod_{j=1}^{13}[a_j,b_j]=t^{12}
\right\rangle.
$$

The sign of the exponent depends on the orientation convention. All orders and dimensions below are unchanged if that convention is reversed.

In the abelianization all commutators vanish, leaving 12t=0. Smith normal form gives one nonzero invariant factor, 12:

$$
H_1(P;\mathbb Z)=\mathbb Z^{26}\oplus\mathbb Z/12.
$$

Independently, the Gysin sequence contains multiplication by the Euler class:

$$
H^0(X;\mathbb Z)\xrightarrow{\ \times12\ }H^2(X;\mathbb Z)
\xrightarrow{\ \pi^*\ }H^2(P;\mathbb Z)
\longrightarrow H^1(X;\mathbb Z)\longrightarrow0.
$$

It follows that

$$
\begin{aligned}
H^0(P;\mathbb Z)&=\mathbb Z,\\
H^1(P;\mathbb Z)&=\mathbb Z^{26},\\
H^2(P;\mathbb Z)&=\mathbb Z^{26}\oplus\mathbb Z/12,\\
H^3(P;\mathbb Z)&=\mathbb Z.
\end{aligned}
$$

The direct-sum splitting of H2 cohomology is not canonical. The image of the base degree class is its canonical torsion subgroup. In particular,

$$
c_1(\pi^*L_m)=m\bmod12,\qquad
\operatorname{ord}\bigl(c_1(\pi^*L_m)\bigr)
=\frac{12}{\gcd(12,m)}.
$$

Degree three becomes order four; degree six becomes order two; degree twelve becomes topologically trivial. The fiber is infinite in the fundamental group, but has order twelve in first homology. Those are different statements.

## 4. The connection and the holonomy

Use a real unitary connection convention with c1(L) represented by F/(2 pi) and holonomy exp(i integral A). Let Theta be a principal connection on P, normalized by

$$
\oint_t\Theta=2\pi,\qquad d\Theta=\pi^*F_{S_0}.
$$

On the curve choose constant-curvature line-bundle connections with

$$
F_{L_m}=\frac{m}{12}F_{S_0}.
$$

This is possible for every degree m; flat twists on the base remain available. The pullback connection pi*A_m is generally curved. On the **same pulled-back line bundle**, change it by a globally defined one-form:

$$
\widehat A_m=\pi^*A_m-\frac{m}{12}\Theta.
$$

Then

$$
F_{\widehat A_m}=0,\qquad
\lambda_m:=\operatorname{Hol}_t(\widehat A_m)
=\exp\!\left(-\frac{2\pi i m}{12}\right).
$$

This gives a concrete flat representative of the torsion bundle. A flat connection can have nonzero torsion first Chern class: vanishing curvature detects real cohomology, which forgets torsion.

For the degrees relevant here:

| Oriented degree m | Torsion order | Fiber holonomy |
|---|---:|---|
| 0 | 1 | 1 |
| +3 | 4 | −i |
| −3 | 4 | +i |
| +6 or −6 | 2 | −1 |
| 12 | 1 | 1 |

**The holonomy cannot be tuned away while keeping the same line bundle and a flat connection.** Two unitary flat connections on that fixed bundle differ by a closed one-form. Its integral around t is zero because t is torsion in integral first homology. A base Wilson line therefore cannot change lambda.

For commuting smooth Higgs backgrounds, each diagonal phi is also a closed one-form, so its fiber period vanishes for the same reason. Passing to the flat complex connection does not remove the fiber monodromy. This explains why the conclusion applies beyond the special choice phi=0.

## 5. The decisive theorem: nontrivial fiber holonomy removes all bulk zero modes

**Fiber-acyclicity statement.** Let a rank-one local system on P have fiber holonomy lambda different from one. Then

$$
H^k(P;\mathcal L)=0\qquad\text{for every }k.
$$

**Proof.** On a circle, the local-system cochain complex is

$$
\mathbb C\xrightarrow{\ \lambda-1\ }\mathbb C.
$$

When lambda is not one, this map is invertible, so both fiber cohomology groups vanish. The fiber-cohomology spectral sequence for P over X therefore has a zero second page and zero total cohomology. This argument allows arbitrary compatible holonomy around the base cycles.

For the declared twisted 7D construction, harmonic representatives of these cohomology groups are the massless internal modes. Thus this result removes the modes themselves; it says more than merely that their net chiral index is zero.

The factor S0 used in the original six-dimensional spinor bundle does not restore them. The pullback pi*S0 is topologically trivial, and its degree twelve gives fiber holonomy one in this flat construction. More fundamentally, the standard 7D twist uses a complex of forms with gauge coefficients; it is not the old Dolbeault complex with one extra coordinate attached.

### Independent exact computation

The included program evaluates Fox derivatives of the fundamental-group presentation. With 27 generators and 27 displayed relators it constructs

$$
C^0\xrightarrow{d_0}C^1\xrightarrow{d_1}C^2,
\qquad
\dim C^0=1,\quad \dim C^1=\dim C^2=27.
$$

For a character chi, d0 has entries chi(gj)−1 and d1 has the evaluated Fox derivatives of each relator. The program checks every relator and d1 d0=0 exactly.

This presentation complex determines H0 and H1. The missing top-dimensional information is supplied by **Poincare duality on the actual closed oriented 3-manifold**, using the inverse character:

$$
\dim H^2(P;\mathcal L_\chi)=\dim H^1(P;\mathcal L_{\chi^{-1}}),
\qquad
\dim H^3(P;\mathcal L_\chi)=\dim H^0(P;\mathcal L_{\chi^{-1}}).
$$

We do not treat the two-dimensional presentation complex as a complete three-dimensional cellular model.

All twelve residues were checked with three base-character choices: trivial, unitary nontrivial, and complex nonunitary. These are 36 primary configurations, with inverse characters also evaluated for duality. Computation is exact over the field generated by i and square root of three; no floating-point rank tolerance is used.

| Fiber / base character | Rank d0 | Rank d1 | Dimensions H0, H1, H2, H3 |
|---|---:|---:|---|
| Fiber 1; trivial base | 0 | 1 | (1,26,26,1) |
| Fiber 1; nontrivial base | 1 | 2 | (0,24,24,0) |
| Any nontrivial twelfth root on fiber | 1 | 26 | (0,0,0,0) |

The two controls agree with the untwisted and local-system Gysin calculations. In the nontrivial-base case, H0 and H2 of the surface vanish and its Euler characteristic gives dim H1=24.

## 6. Application to the actual M1 sectors

The program reads the stack data and representation table from the supplied source:

$$
(m_c,m_L,m_a,m_b,m_d)=(0,-3,3,3,0).
$$

For a sector with charges qx, the oriented degree is the sum of qx mx. Negative-degree sectors in the source table are conjugated to present their original chiral representation. The code accounts for that conjugation when assigning holonomy.

| Original sector | Original representation (SU3, SU2, Y) | Oriented degree / original multiplicity | Fiber holonomy | 7D H1 and H2 |
|---|---|---:|---|---|
| (c, L-bar) | (3,2,1/6) | 3 | −i | 0,0 |
| (c-bar, a) | (3-bar,1,−2/3) | 3 | −i | 0,0 |
| (c-bar, b) | (3-bar,1,1/3) | 3 | −i | 0,0 |
| (L-bar, a) | (1,2,−1/2) | 6 | −1 | 0,0 |
| (L-bar, b) | (1,2,+1/2) | 6 | −1 | 0,0 |
| (L-bar, d) | (1,2,−1/2) | 3 | −i | 0,0 |
| (a, d-bar) | (1,1,0) | 3 | −i | 0,0 |
| (b, d-bar) | (1,1,+1) | 3 | −i | 0,0 |

Thus every nonzero-degree M1 sector is acyclic. Neither its three-family matter nor its degree-six doublet sectors survive as massless bulk modes of this candidate.

There are also two degree-zero off-diagonal pairs: (c,d) and (a,b). Their cohomology is generally nonzero, as shown by the controls. These are not necessarily gauge-neutral: their displayed orientations have representations (3,1,2/3) and (1,1,−1). If they belong to the proposed adjoint decomposition, they can supply unwanted bulk matter. Adjoint sectors also survive. Therefore the result is **not** that the entire theory becomes empty or all charged matter vanishes. It fails to reproduce the required spectrum.

The original six-dimensional equality

$$
\chi\!\left(X,S_0\otimes L_x\otimes L_y^{-1}\right)=m_x-m_y
$$

remains valid in its original problem. The new calculation has a different space and a different mode operator. There is no contradiction between them.

## 7. What survives of the arithmetic spin structure

Using the theta-compatible normalization from the preceding audit, the half-form matrices are

$$
\widetilde A=
\begin{pmatrix}i&0\\0&-i\end{pmatrix},\qquad
B=
\begin{pmatrix}0&-1/\sqrt{13}\\ \sqrt{13}&0\end{pmatrix}.
$$

They satisfy

$$
\widetilde A^2=B^2=(\widetilde AB)^2=-I,\qquad
\widetilde AB\widetilde A^{-1}B^{-1}=-I.
$$

Exact group closure gives eight elements, with one identity, one element of order two, and six of order four: Q8. Their symmetric squares match the previously established action on the three products of half-forms. This retains the correction distinguishing an arbitrary line-bundle lift from a lift compatible with the fixed map S0 squared to K. The present study rechecks those matrices and their square-map algebra; it does not independently repeat the earlier eta-multiplier derivation.

Averaging a Hermitian metric over this finite group realizes the action on the unit-circle bundle. Its central minus identity is a half-turn of the circle, isotopic to the identity. Consequently it acts trivially on ordinary real cohomology. Since pi* identifies H1(X;R) with H1(P;R), the pulled-back Atkin–Lehner traces are 2, −2 and −18 for W11, W13 and W143. These were checked against the supplied exact weight-two matrices.

This establishes a symmetry of the spin geometry. The original odd-degree M1 CM-divisor background breaks W11. Full Q8 equivariance of a newly chosen flat 7D background has **not** been established, nor does geometric Q8 by itself restore the original flux-divisor data.

Ordinary Fourier functions do retain useful bundle information:

$$
f(p\,e^{i\theta})=e^{-in\theta}f(p)
\quad\Longleftrightarrow\quad
f\text{ represents a smooth section of }S_0^{\otimes n}.
$$

Holomorphicity is a further condition on the base. This realizes the familiar half-forms kinematically, but a Fourier section need not be a zero mode of the physical 7D operator. Confusing these two spaces would incorrectly preserve the three-family count.

A full lift of the Hecke correspondences, with compatible pull-push maps and the relevant operator domains, remains unconstructed in this study.

## 8. The metric and gravity test

This is a **separate conditional test**: add an unwarped seven-dimensional Einstein action and the stated ordinary nonnegative classical energy terms. We have not embedded the gauge construction into a complete supergravity solution.

Let the compact base have Gaussian curvature −1/R squared and let r be the circle radius. Its area is 48 pi R squared. Use a uniform-curvature connection:

$$
ds_P^2=ds_X^2+r^2\Theta^2,\qquad
F_\Theta=\frac{e}{24R^2}\operatorname{vol}_X,\qquad
e=12.
$$

The volume is

$$
\operatorname{Vol}(P)=96\pi^2R^2r.
$$

Locally in upper-half-plane coordinates one may write

$$
ds_P^2=
R^2\frac{dx^2+dy^2}{y^2}
r^2\left(d\theta+\frac12\frac{dx}{y}\right)^2.
$$

The program computes the Christoffel symbols and Ricci tensor directly. The scalar curvature is

$$
\boxed{\mathcal R_P=-\frac{2}{R^2}-\frac{r^2}{8R^4}.}
$$

It agrees with the submersion formula. In an orthonormal frame, the Ricci eigenvalues are

$$
\left(
-\frac1{R^2}-\frac{r^2}{8R^4},\
-\frac1{R^2}-\frac{r^2}{8R^4},\
+\frac{r^2}{8R^4}
\right).
$$

Hence this connection metric is neither Einstein nor constant-negative-curvature. More generally, P cannot be a closed hyperbolic 3-manifold: its fundamental group has a nontrivial infinite center coming from the fiber. This spin-circle route is not the same construction as the Kleinian/Bers route to hyperbolic three-dimensional geometry.

**Compactness matters.** The metric here uniformizes the compact curve after filling its four cusps. It is not the usual complete finite-area metric on the punctured modular curve Y0(143), whose area at curvature −1 is 56 pi. Removing the cusps changes the topology and the mode problem; it is not an innocuous substitution into this calculation.

If the Einstein term is normalized by M7 to the fifth power divided by two, dimensional reduction gives

$$
M_4^2=M_7^5\,96\pi^2R^2r,\qquad
\frac1{g_4^2}=\frac{96\pi^2R^2r}{g_7^2}
$$

for a common bulk gauge normalization. These are reduction identities, not predictions of either coupling or an absolute scale.

### Overall-size obstruction in the declared classical ansatz

Scale both radii by a dimensionless rho while keeping their ratio fixed. After the four-dimensional Weyl rescaling to Einstein frame, the contributions from negative internal curvature, fixed internal two-form flux and nonnegative bulk vacuum energy scale as

$$
V_E(\rho)=\frac{a}{\rho^5}
          +\frac{b}{\rho^7}
          +\frac{c}{\rho^3},
\qquad a>0,\quad b,c\ge0.
$$

The volume scales as rho cubed, curvature as rho to minus two and the squared norm of a fixed two-form as rho to minus four. The Einstein-frame Weyl factor supplies rho to minus six, producing these exponents.

Therefore

$$
\frac{dV_E}{d\rho}
=-\frac{5a}{\rho^6}-\frac{7b}{\rho^8}
-\frac{3c}{\rho^4}<0.
$$

There is no stationary overall size in this restricted ansatz. For the actual flat gauge backgrounds above, their Yang–Mills flux contribution has b=0; torsion holonomy itself has no curvature energy. The nonzero circle-bundle curvature is already included in the geometric coefficient a.

This is an obstruction to this particular classical reduction. It is not a theorem about all seven-dimensional theories, all supergravity potentials, or theories with additional sources, warping, quantum corrections, or other fields.

## 9. What this decides and what it does not

**Decided within the declared model:** a smooth spin-circle construction with the original independent line bundles and the standard commuting 7D twisted mode complex does not preserve M1's three families. Changing radii or base Wilson lines cannot change its fiber-acyclicity result. The separate restricted gravity ansatz cannot stabilize the overall size.

**Still open:** defect or boundary sectors, different compactifications, and different parent actions. Noncommuting backgrounds can change individual cohomology groups, but do not automatically solve chirality. For any finite-rank flat local system on a closed 3-manifold,

$$
\chi(P;\mathcal L)=\operatorname{rank}(\mathcal L)\chi(P)=0,
\qquad h^1-h^2=h^0-h^3.
$$

In a matter sector with H0 and H3 absent this already forces equal chiral and conjugate-chiral dimensions. A proposed noncommuting alternative must address this broader index constraint as well as the specific fiber obstruction. No such alternative was inserted after the failed gate.

The source-free calculation identifies the next worthwhile target: a mechanism that changes the physical mode complex. Defects and boundary conditions are established ingredients in 7D Higgs-bundle constructions, but their spectrum and consistency must be calculated afresh; a singularity is not itself a prediction of three families. See [Pantev and Wijnholt](https://arxiv.org/abs/0905.1968) and [the relative-cohomology treatment in section 4.1 of Braun et al.](https://arxiv.org/html/1812.06072v2).

The four cusps of the arithmetic curve provide distinguished places to investigate possible defect data. A follow-up would have to specify their sources, bundles, boundary conditions and mode operator before counting states. It would also have to handle the missing bulk modes; merely adding a local calculation to the old six-dimensional count would not do so.

A useful next study would choose one such mechanism, derive its appropriate index and anomaly conditions, and require the complete M1 representation content with net multiplicity three. Its boundary data must be justified rather than chosen to insert the answer.

The Fano-plane and exceptional-group connections can remain algebraic tools. They do not select this spacetime dimension or cure the calculated zero-mode obstruction. Likewise, a 7D gauge worldvolume here is distinct from the seven internal dimensions of an eleven-dimensional M-theory compactification.

## 10. Reproducibility and evidence status

The supplied archive has SHA-256

~~~
2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5
~~~

The bundle contains the preregistration, runnable script, result ledger, execution log, dependency versions, and unchanged input archive. The script hashes the imported source files before and after calculation.

Exact checks include:

- genus from the modular index and elliptic/cusp counts, checked against the supplied weight-two basis dimension;
- integral Smith normal form and the independent Gysin calculation;
- all twelve holonomy residues under three base-character choices;
- every nonzero-degree M1 representation and all ten unordered stack pairs;
- the theta-compatible group algebra and symmetric-square identities;
- direct scalar/Ricci curvature computation against the submersion expression;
- strict negativity of the declared Einstein-frame scale derivative.

The fiber-acyclicity theorem supplies the result for arbitrary compatible base holonomies; the finite scan is its independent computational check, not a claim to have sampled every connection.

Two intermediate development runs stopped because assertions compared different factorizations of the same symbolic expression. The assertions were changed to simplify their differences to zero. The curvature formulas and registered physical criteria were unchanged. The completed execution is the one recorded in the bundled result ledger and run log.

No finite-element estimate, measured particle mass, complete gravitational solution, or ultraviolet completion is certified. The exact statements concern the specified mathematics; their physical interpretation is conditional on the chosen 7D field theory.

## 11. Short glossary

| Term | Meaning in this investigation |
|---|---|
| Theta characteristic | A holomorphic line bundle S with S squared isomorphic to the canonical bundle K of a curve. |
| Unit-circle bundle | The manifold of norm-one vectors in a complex line bundle; each point of the base has a circle above it. |
| Euler class | The integral class measuring the twisting of an oriented circle bundle; here its degree is 12. |
| Torsion class | A nonzero integral cohomology class killed by an integer; invisible to real differential-form cohomology. |
| Flat connection | A connection with zero curvature, which can nevertheless have nontrivial parallel transport around loops. |
| Holonomy | The transformation acquired by parallel transport around a closed loop; here a complex phase on the fiber. |
| Local system | A vector space transported consistently along paths, specified by a representation of the fundamental group. |
| Acyclic | Having zero cohomology in every degree. |
| Partial topological twist | A reorganization using an internal rotation group and R-symmetry that puts the supersymmetric fields into differential-form bundles. |
| BPS equations | Background equations required by the specified preserved supersymmetry. |
| Higgs one-form phi | An adjoint-valued internal one-form in the twisted theory; not automatically the observed electroweak Higgs. |
| Fox derivative | A group-ring derivative of a word, used here to construct a twisted cochain differential from group relations. |
| Gysin sequence | An exact sequence connecting the cohomology of an oriented circle bundle and its base. |
| Chiral index | A signed count of chiral and conjugate-chiral modes; it is different from the total number of modes. |
| Radion | A four-dimensional field describing an internal size, such as the overall scale rho. |
| Einstein frame | The four-dimensional metric convention in which the Einstein curvature term has a constant coefficient. |
| Hecke correspondence | An arithmetic correspondence inducing pull-push operators; extending it to a new physical construction requires compatible bundle and operator data. |
