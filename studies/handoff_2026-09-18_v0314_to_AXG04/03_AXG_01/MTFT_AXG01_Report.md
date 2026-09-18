# MTFT AXG-01: two axions, anomaly cancellation, and the mass–gravity question

**Investigation date:** 17 September 2026  
**Input:** supplied MTFT 0.31.4 candidate, representative five-stack model M1  
**Status:** exact charge and anomaly algebra; conditional four-dimensional effective theory; no derived vacuum or absolute mass scale

## Abstract

We investigate whether the two anomalous Abelian directions of the declared MTFT five-stack model can support an axion sector connecting gauge consistency, Higgs interactions and internal-radius dynamics. The complete four-dimensional perturbative anomaly polynomial factors through a primitive rank-two integer charge matrix. This supplies a local Green–Schwarz/Stückelberg cancellation construction, including a gravitational Pontryagin coupling. Canonical normalization gives two massive vector combinations and three massless combinations before additional symmetry breaking. An axion-dressed Higgs bilinear is gauge invariant, but its flavor coefficients are not determined. The same charge lattice permits dressed baryon-number-violating operators. The minimal kinetic extension contributes no classical radius potential and leaves the previously identified radius runaway unchanged. These results narrow the next parent-theory calculation; they do not establish a six-dimensional completion, a stable compactification, electroweak symmetry breaking, or gravity unification.

The accompanying scripts execute **159 exact assertions: 106 anomaly/lattice/mass checks, 21 operator checks and 32 vacuum checks**. Their count measures reproducibility of the stated algebra, not physical evidence for the model.

## 1. What this investigation establishes

| Question | Result | Qualification |
|---|---|---|
| Can two axions cancel the declared 4D local anomalies? | The full anomaly polynomial factors through the two shift charges | Requires the displayed Wess–Zumino terms and an appropriate local counterterm prescription |
| Do two vector combinations acquire mass? | Yes, for positive definite gauge and axion kinetic matrices | Their absolute masses depend on undetermined kinetic data |
| Does this introduce a gravity coupling? | The cancellation includes an axion coupling to the gravitational Pontryagin four-form | It does not generate an Einstein action or determine Newton's constant |
| Is an axion-dependent Higgs bilinear allowed? | Yes | Its coefficient, flavor structure and dynamical origin are additional data |
| Does the minimal extension stabilize the radius? | No, within the stated classical ansatz | A more complete scalar sector may change the result |
| Do the extra Abelian D-terms select a Higgs flavor direction? | Not with identical charges and flavor-symmetric quadratic kinetic terms alone | Family-sensitive interactions can lift those directions |
| Is a 6D parent or quantum-gravity completion established? | No | 4D local anomaly cancellation is only one necessary consistency calculation |

The internal space remains the real two-dimensional curve $X_0(143)$. Axion gauge fibers are not added physical spacetime coordinates. The effective theory below is four-dimensional.

## 2. Inputs, charge conventions and source audit

Use stack order $(c,L,a,b,d)$ and fundamental central charge $+1$. The supplied package declares

$$
G=U(3)_c\times U(2)_L\times U(1)_a\times U(1)_b\times U(1)_d,
\qquad N=(3,2,1,1,1),\qquad m=(0,-3,3,3,0).
$$

The index of the oriented bifundamental $(i,\bar j)$ is $m_i-m_j$. A negative index is represented by left-handed fermions in the conjugate representation. These indices, including the extra doublet sectors, determine the anomaly sums. They are not replaced by the Standard Model spectrum alone.

Three distinguished generator vectors are

$$
p=\begin{pmatrix}1\\1\\1\\1\\1\end{pmatrix},\qquad
Y=\begin{pmatrix}1/6\\0\\-1/2\\1/2\\-1/2\end{pmatrix},\qquad
B-L=\begin{pmatrix}1/3\\0\\0\\0\\-1\end{pmatrix}.
$$

For a field with charge vector $q_\Phi$, its hypercharge is $q_\Phi^{\mathsf T}Y$. The common phase $p$ acts trivially on all these bifundamentals; its gauge field is not thereby automatically removed from an action.

The candidate axion charge map is

$$
\boxed{K=\begin{pmatrix}0&-2&1&1&0\\3&-4&0&0&1\end{pmatrix}.}
$$

Introduce dimensionless, $2\pi$-periodic axions $a^I$, $I=1,2$, with

$$
A\mapsto A+d\lambda,\qquad a\mapsto a-K\lambda,
\qquad Da=da+KA.
$$

Direct multiplication gives $Kp=KY=K(B-L)=0$. Their independence and $\operatorname{rank}K=2$ imply

$$\ker_{\mathbb R}K=\operatorname{span}_{\mathbb R}\{p,Y,B-L\}.$$

This specifies an effective extension. It does not derive axions from a higher-dimensional tensor field.

**What is actually in 0.31.4.** `src/mtft/surface/smflux.py` supplies the stack data, anomaly sums and Majorana obstruction. `docs/SM/SM01_PREREGISTRATION.md` and `SM01_REPORT.md` treat Green–Schwarz cancellability as a parent assumption. They do not supply the axion kinetic matrix or a complete potential. The inspected archive mentions CW-01 in the review response but does not include a corresponding derivation sufficient to reconstruct a supersymmetric parent. We therefore do not assume that every Higgs direction is protected merely from an $N=1$ label.

The older `src/mtft/cosmology.py` is a separate phenomenological construction: it takes Newton's constant and temporal-axion parameters as inputs. It does not derive the cancellation or radius potential studied here. Likewise, `surface/condensation.py` concerns a separate bundle-extension calculation and does not by itself provide this five-stack potential. Retracted SM-14 flavor fits are not inputs to AXG-01.

## 3. Full anomaly factorization, reproducible by hand

### 3.1 The cubic polynomial

Let $q=(q_c,q_L,q_a,q_b,q_d)$ denote an arbitrary generator. Define

$$x=q_c-q_L,\quad y=q_a-q_L,\quad z=q_b-q_L,\quad w=q_d-q_L.$$

The cubic trace over the net left-handed spectrum is

$$P(q)=\sum_{i<j}(m_i-m_j)N_iN_j(q_i-q_j)^3.$$

For these inputs it can be written without any hidden data as

$$
\begin{aligned}
P={}&18x^3-9(x-y)^3-9(x-z)^3\\
&+12y^3+12z^3+6w^3+3(y-w)^3+3(z-w)^3.
\end{aligned}
$$

Set

$$r_1=y+z=(Kq)_1,\qquad r_2=3x+w=(Kq)_2.$$

Expanding the cubes and grouping terms gives the exact identity

$$\boxed{P=r_1Q_1+r_2Q_2,}$$

$$
Q_1=24(y^2-yz+z^2)+27x^2+9w^2,
\qquad Q_2=-9(y^2+z^2).
$$

For example, use $(y+z)(y^2-yz+z^2)=y^3+z^3$. The remaining terms are
$27x^2(y+z)+9w^2(y+z)-9(3x+w)(y^2+z^2)$.
This gives a short hand check of the complete cubic identity.

### 3.2 Mixed anomaly matrix

Normalize the fundamental non-Abelian Dynkin index to $T(\mathrm{fund})=1/2$. Then

$$
\begin{pmatrix}A_3(q)\\A_2(q)\\A_{\mathrm{grav}}(q)\end{pmatrix}
=
\begin{pmatrix}
0&-3&3/2&3/2&0\\
9/2&-12&3&3&3/2\\
0&-48&24&24&0
\end{pmatrix}q
=
\begin{pmatrix}
\tfrac32r_1\\3r_1+\tfrac32r_2\\24r_1
\end{pmatrix}.
$$

Here $A_3$ and $A_2$ are the $SU(3)^2U(1)$ and $SU(2)^2U(1)$ traces; $A_{\mathrm{grav}}$ is the linear charge trace. The package's `mixed_nonabelian_anomalies` uses twice the first two coefficients, whereas `anomaly_ledger` uses the convention displayed here. The scripts check both conventions.

The pure $SU(3)^3$ anomaly vanishes. The net chiral spectrum contains 24 Weyl $SU(2)$ doublets when color multiplicity is included, so its Witten parity is even. Additional genuinely vector-like pairs do not change either conclusion.

### 3.3 Surviving currents still require careful anomaly allocation

The polynomial vanishes on $\ker K$, but mixed cubic traces involving two surviving generators and an arbitrary third one need not vanish. With
$A(u,v,w)=\operatorname{tr}_{\mathrm{chiral}}(Q_uQ_vQ_w)$,

$$
A(Y,Y,q)=7r_1-\frac32r_2,\qquad
A(Y,B-L,q)=2r_1,\qquad
A(B-L,B-L,q)=4r_1.
$$

Thus “the cubic vanishes on the surviving span” must not be read as “all its mixed current anomalies vanish in every prescription.” A local Bardeen/generalized Chern–Simons counterterm, or an equivalent regulator prescription, redistributes the consistent anomaly so that the surviving gauge Ward identities are preserved. Such terms are part of the general multi-axion effective action; anomaly cancellation does not fix every allowed interaction. See [Anastasopoulos et al., §§1–2](https://arxiv.org/abs/hep-th/0605225).

There is an explicit hand-checkable descent identity. Temporarily take $A_a$ to be normalized one-forms with $dA_a=f_a$. For the unnormalized cubic polynomial $P(f)$, define

$$
I^{\mathrm{sym}}_5=\frac13\sum_a A_a\,\partial_aP(f),
\qquad I^{\mathrm{fac}}_5=\sum_I(KA)_I Q_I(f),
$$

$$C_4=\frac13\sum_{I,a}(KA)_I\wedge A_a\,\partial_aQ_I(f).$$

Euler's identity for homogeneous polynomials and the graded product rule give

$$dI^{\mathrm{sym}}_5=dI^{\mathrm{fac}}_5=P(f),\qquad
I^{\mathrm{sym}}_5-I^{\mathrm{fac}}_5=dC_4.$$

For the physical $P/6$ convention below, divide these three expressions by six. This is a local construction; it is not a uniqueness claim for counterterms or a global extension theorem.

## 4. The local cancellation action and its gravity coupling

Use normalized Abelian curvatures $f_a=F_a/(2\pi)$ and characteristic classes with
$\mathrm{ch}_2(\text{SU fundamental})=-c_2$ and
$\widehat A=1-p_1/24+\cdots$. The four-dimensional chiral anomaly is encoded by the formal six-form

$$
I_6=\frac16P(f)-2A_3(f)c_{2,3}-2A_2(f)c_{2,2}
-\frac1{24}A_{\mathrm{grav}}(f)p_1(T).
$$

Substituting the identities above yields

$$\boxed{I_6=(Kf)_1X_4^1+(Kf)_2X_4^2,}$$

$$
X_4^1=\frac16Q_1(f)-3c_{2,3}-6c_{2,2}-p_1(T),
\qquad
X_4^2=\frac16Q_2(f)-3c_{2,2}.
$$

In a factorized descent convention, take the fermionic anomalous variation of the Euclidean effective action to be
$2\pi i\int(K\lambda/2\pi)_I X_4^I$. Then the local Wess–Zumino term

$$S^{E}_{\mathrm{WZ}}=2\pi i\int_{M_4}\sum_I\frac{a^I}{2\pi}X_4^I$$

has the opposite variation because $\delta a=-K\lambda$. The counterterm described in §3 relates this convention to the symmetric Abelian descent. The common overall sign changes if the chirality/anomaly convention changes.

### 4.1 A concrete connection to curvature

The first axion has a term proportional to

$$\boxed{a^1p_1(T)\ \propto\ a^1\operatorname{tr}(\mathcal R\wedge\mathcal R).}$$

This follows here from the mixed gauge–gravitational anomaly. It is a coupling to the gravitational Pontryagin density, often expressed as $a^1 R\widetilde R$. It can be defined with a background metric. Making that metric dynamical requires a gravitational action and consistency conditions of its own. A Pontryagin coupling and the Einstein–Hilbert term $\int\sqrt{-g}\,\mathcal R$ are different terms; the former does not supply the latter's kinetic normalization or Newton's constant. For an example of studying Pontryagin terms as additions to an already specified gravitational theory, see [Jackiw and Pi](https://arxiv.org/abs/gr-qc/0308071).

Both axions in this minimal construction are eaten locally. A unitary-gauge description must retain the equivalent anomaly-canceling interactions; it does not turn them into two free observable pseudoscalars or into a spin-two graviton.

### 4.2 What is and is not globally checked

On a closed spin four-manifold, the square of an integral degree-two class has even integral. For split $SU(3)\times SU(2)\times U(1)^5$ cover bundles, the half-integral diagonal coefficients in $Q_I/6$ therefore have integral periods; mixed coefficients and the displayed characteristic-class coefficients are integral. The proposed $2\pi$ axion periods pass this restricted test.

General bundles for the quotient groups underlying $U(3)$ and $U(2)$ need a separate global analysis. The restricted period calculation does not establish all large-gauge transformations, torsion effects, the six-dimensional anomaly polynomial, or a UV completion.

## 5. Vector masses: use the kinetic metric

Write the Abelian kinetic and Stückelberg terms as

$$
\mathcal L=-\frac14H_{ab}F^a_{\mu\nu}F^{b\mu\nu}
-\frac12G_{IJ}(Da^I)_\mu(Da^J)^\mu,
$$

where $H$ and $G$ are positive definite. $H$ is dimensionless and $G$ has mass dimension two in these conventions. In canonical gauge coordinates the squared-mass matrix is

$$\boxed{\mathcal M^2=H^{-1/2}K^{\mathsf T}GK H^{-1/2}.}$$

It has rank two. The remaining three zero modes correspond to $p,Y,B-L$ in the original connection coordinates; canonical-coordinate vectors are obtained with the appropriate $H^{1/2}$ factor. Kinetic mixing can change their normalization and couplings without changing this rank statement.

For a common stack coupling and fundamental central charges $+1$, the package's normalization gives

$$H=\frac{H_0}{g_4^2},\qquad H_0=\operatorname{diag}(6,4,2,2,2).$$

Choose, **only as an illustrative extra input**, $G=f_{\mathrm{ax}}^2I_2$. The nonzero eigenvalues reduce to those of

$$
g_4^2f_{\mathrm{ax}}^2KH_0^{-1}K^{\mathsf T}
=g_4^2f_{\mathrm{ax}}^2\begin{pmatrix}2&2\\2&6\end{pmatrix}.
$$

The dimensionless characteristic polynomial is $\lambda^2-8\lambda+8$, hence

$$\boxed{m_\pm^2=g_4^2f_{\mathrm{ax}}^2(4\pm2\sqrt2).}$$

If one instead sets both kinetic metrics to the identity in unit scales, the control matrix is

$$KK^{\mathsf T}=\begin{pmatrix}6&8\\8&26\end{pmatrix},$$

whose eigenvalues are $16\pm2\sqrt{41}$. This is a different kinetic assumption, not a competing physical prediction. Neither calculation fixes $g_4$, $f_{\mathrm{ax}}$, their radius dependence, or a mass in GeV.

### 5.1 Integral charges and possible remnants

Columns $a,d$ of $K$ form the identity matrix. Consequently its Smith invariant factors are $(1,1)$, so on the primitive $U(1)^5$ torus the kernel is connected: there is no extra finite Stückelberg remnant forced by this charge map. An integral kernel basis is

$$
\begin{pmatrix}
1&0&0\\0&1&0\\0&0&1\\0&2&-1\\-3&4&0
\end{pmatrix}.
$$

The shifts also respect the stated unitary centers at character level:

$$\chi_1=(\det U_L)^{-1}z_a z_b,\qquad
\chi_2=(\det U_c)(\det U_L)^{-2}z_d.$$

Thus their color and weak central charges have the required multiples of three and two. This checks the shift map, not the unresolved global quantization of the full anomaly action.

## 6. A Higgs interaction is allowed—and the same test finds other operators

Let $e_i$ be the standard charge basis. The matter and candidate Higgs scalar charges used by the triangle Yukawa couplings are

| Field | Charge vector expression | Hypercharge | $B-L$ |
|---|---|---:|---:|
| $Q$ | $e_c-e_L$ | $1/6$ | $1/3$ |
| $u^c$ | $e_a-e_c$ | $-2/3$ | $-1/3$ |
| $d^c$ | $e_b-e_c$ | $1/3$ | $-1/3$ |
| $L$ | $e_d-e_L$ | $-1/2$ | $-1$ |
| $e^c$ | $e_b-e_d$ | $1$ | $1$ |
| $\nu^c$ | $e_a-e_d$ | $0$ | $1$ |
| $H_u$ | $e_L-e_a$ | $1/2$ | $0$ |
| $H_d$ | $e_L-e_b$ | $-1/2$ | $0$ |

The Higgs rows specify scalar charges required by these Yukawas; this is not a derivation of a supersymmetric pairing with the net chiral fermion table. All four operators $Qu^cH_u$, $Qd^cH_d$, $Le^cH_d$ and $L\nu^cH_u$ have zero total stack charge.

For an operator $\mathcal O$ and integer $n\in\mathbb Z^2$,

$$\boxed{e^{in_Ia^I}\mathcal O\text{ is gauge invariant}\iff q_{\mathcal O}=K^{\mathsf T}n.}$$

Write $k_1,k_2$ for the row charges of $K$. Direct addition gives

| Operator | Stack charge | Gauge-invariant dressing |
|---|---|---|
| $H_u\cdot H_d$ | $-k_1$ | $e^{-ia^1}H_u\cdot H_d$ |
| $QQQL$ | $k_2$ | $e^{ia^2}QQQL$ |
| $u^cu^cd^ce^c$ | $2k_1-k_2$ | $e^{i(2a^1-a^2)}u^cu^cd^ce^c$ |
| $\nu^c\nu^c$ | Outside the row span; $B-L=2$ | These axions cannot neutralize it |

Non-Abelian and Lorentz indices in the four-fermion operators are contracted to singlets; suitable flavor choices avoid antisymmetry zeros.

The potential can therefore contain

$$V_B=\sum_{ij}B_{ij}(\rho,\ldots)e^{-ia^1}H_{u,i}\cdot H_{d,j}+\mathrm{h.c.},$$

where $B_{ij}$ has mass dimension two. This is a concrete candidate interaction to derive from a parent. Gauge charges alone do not determine $B_{ij}$, select a Higgs direction, or ensure a bounded potential. If a supersymmetric completion supplies a chiral field $T_1=s_1+ia^1$, then a holomorphic term $\mu_{ij}e^{-T_1}H_{u,i}\cdot H_{d,j}$ is allowed by the same charges. Its existence, coefficients and saxion dependence require that completion.

The two dressed four-fermion operators are allowed with coefficients proportional to $1/\Lambda^2$ in a nonsupersymmetric effective Lagrangian. They violate baryon number while preserving $B-L$. This is a selection-rule result, **not** a prediction of a nonzero coefficient or proton lifetime. A parent that generates the Higgs interaction must also determine which other allowed amplitudes it generates and their suppression.

A bare Majorana mass remains forbidden while $B-L$ is unbroken. These two axions leave that gauge direction unshifted; they do not replace a $B-L$-breaking ingredient.

There is also a source correction worth retaining: the older SM-01 discussion describes an $(a,\bar b)$ singlet as hypercharge-neutral. The declared charges instead give $Y_a-Y_b=-1$. AXG-01 does not use that proposed neutral-singlet mechanism.

## 7. Radius and Higgs vacuum tests

### 7.1 A minimal action makes the missing input explicit

With signature $(-,+,+,+)$, a possible four-dimensional Einstein-frame ansatz is

$$
\begin{aligned}
S_4={}&\int d^4x\sqrt{-g}\bigg[
\frac{M_P^2}{2}\mathcal R-\frac12Z_\rho(\rho)(\partial\rho)^2
-\frac14H_{ab}(\rho)F^aF^b\\
&\hspace{15mm}-\frac12G_{IJ}(\rho)Da^IDa^J
-V_0(\rho)+\mathcal L_{\mathrm{matter}}\bigg]
+S_{\mathrm{WZ}}+S_{\mathrm{ct}}.
\end{aligned}
$$

Here the Einstein term, $Z_\rho$, $G$, and any additional scalar interactions are assumptions or functions to derive. Displaying this action does not claim that the package already implements it. The Lorentzian Wess–Zumino term is the continuation of §4.

In a Lorentz-invariant vacuum with $A_\mu=0$ and constant axions,

$$Da=0\quad\Longrightarrow\quad \mathcal L_{\mathrm{St}}=0,
\qquad \partial_\rho\mathcal L_{\mathrm{St}}=0.$$

Radius-dependent vector masses are not by themselves a classical potential. The anomaly terms likewise do not supply the missing ordinary scalar potential in this background; the gravitational Pontryagin density vanishes on maximally symmetric spacetime.

### 7.2 The inherited restricted radius potential remains a runaway

For the previously considered unwarped product compactification, take a genus-13 internal metric, a two-derivative Einstein term, positive gauge kinetic energy with fixed quantized flux, and nonnegative bulk vacuum energy. Exclude localized sources and further corrections for this test. With $\rho=R/R_0>0$, Weyl rescaling to four-dimensional Einstein frame gives

$$V_0(\rho)=a\rho^{-4}+b\rho^{-6}+c\rho^{-2},
\qquad a>0,\quad b,c\ge0.$$

Here $a$ is the contribution from negative internal curvature, $b$ from flux energy and $c$ from bulk vacuum energy. The coefficients have mass dimension four. Differentiation gives

$$\boxed{V_0'(\rho)=-4a\rho^{-5}-6b\rho^{-7}-2c\rho^{-3}<0.}$$

The minimal kinetic axion extension leaves this result unchanged. It has no stationary finite radius under these assumptions. This does not exclude saxion potentials, warping, other sources, quantum terms or a different parent.

For orientation, at a chosen constant radius the unwarped reduction gives
$\mathcal A_X=48\pi R^2$ and $M_{P,\mathrm{red}}^2=M_6^4\mathcal A_X$ in the corresponding constant-radius metric convention. Those relations still require $M_6$ and $R$. For a fluctuating radius one must perform the Einstein-frame rescaling; an area formula does not fix the physical radius.

### 7.3 Two useful distinctions about axion and Higgs potentials

**Pure axions.** Gauge invariance of a local potential depending only on the two axions implies

$$K^{\mathsf T}\nabla_aV=0.$$

Because $K^{\mathsf T}$ is injective, $\nabla_aV=0$. An isolated nonconstant $\cos a^I$ potential is therefore incompatible with these gauged shifts. Charged matter dressings such as §6 evade this restriction through their simultaneous gauge transformations.

**Flavor directions.** Under a flavor-symmetric quadratic Kähler metric, an extra-Abelian D-term has the form

$$D_a=\xi_a(\rho)+\sum_s q_{a,s}\|H_s\|^2,
\qquad V_D=\frac12D_a(H^{-1})^{ab}D_b.$$

Fields within each identical-charge sector enter through norms, so this interaction alone cannot choose their family orientation. For a hand example, take three real components and

$$V_D=\frac{g^2}{2}\left[\xi+q(h_1^2+h_2^2+h_3^2)\right]^2.$$

At $(h_1,h_2,h_3)=(v,0,0)$ and $\xi=-qv^2$,

$$\operatorname{Hess}V_D=\operatorname{diag}(4g^2q^2v^2,0,0).$$

The norm direction has curvature while the two tangent directions remain flat. This example does not exclude family-sensitive Kähler terms, Yukawa-induced effects, the matrix $B_{ij}$, soft terms or nonperturbative contributions.

A full supersymmetric model may also contain saxions and gauged moment maps that generate a potential already at tree level. It is therefore essential to distinguish the minimal bosonic kinetic extension tested here from a complete supersymmetric action. Supersymmetric scalar potentials require actual Kähler, superpotential and gauge data; see [Martin, §§3.4 and 4.10](https://arxiv.org/abs/hep-ph/9709356).

### 7.4 A control that shows exactly where the sign restriction matters

As a mathematical control only, change $c=-C<0$. This is not a source derived for MTFT. Set $x=\rho^{-2}$, obtaining

$$V=-Cx+ax^2+bx^3,\qquad C,a>0,\ b\ge0.$$

The unique positive stationary point is

$$x_* = \frac{C}{a+\sqrt{a^2+3bC}},\qquad \rho_*=x_*^{-1/2}.$$

This includes the $b=0$ limit. At the stationary point,

$$V_*=-ax_*^2-2bx_*^3<0,\qquad
V''(\rho_*)=\frac{8a}{\rho_*^6}+\frac{24b}{\rho_*^8}>0.$$

Thus a different sign assumption can produce a radial minimum. With positive kinetic energy, Einstein gravity and all other fields stationary, it is an AdS vacuum. This does not establish stability in other field directions or a consistent negative-$c$ compactification source. For $a=b=C=1$ in arbitrary units, $\rho_*^2=3$ and $V_*=-5/27$; these are a hand-check example, not masses or measured scales.

Loop effects require a separate calculation. A vector-only Coleman–Weinberg expression is insufficient: one needs the complete bosonic and fermionic spectrum, gauge fixing and ghosts, interactions, regulator and renormalized local counterterms. Existing KK diagnostics do not automatically provide that complete input. Controlled radius stabilization by competing loop contributions is possible in other specified models, as illustrated by [von Gersdorff and Hebecker](https://arxiv.org/abs/hep-th/0504002); their result is not a potential for this curve.

## 8. How this bears on particle masses and gravity

There are now three distinct mechanisms to keep separate.

In the mass formula below, define the neutral Higgs expectation values by $\langle H^0_{u,\alpha}\rangle=v_{u,\alpha}/\sqrt2$.

| Mechanism | Object determining the mass/coupling | What AXG-01 supplies |
|---|---|---|
| Extra vector mass | Canonical eigenvalues of $H^{-1/2}K^{\mathsf T}GKH^{-1/2}$ | Exact rank and conditional two-by-two spectrum |
| Fermion mass after Higgs condensation | $M^u_{ij}=\sum_\alpha Y^u_{ij\alpha}(\rho)v_{u,\alpha}/\sqrt2$, and analogous sectors | An allowed interaction that could help select Higgs combinations; no derived $v_\alpha$ |
| Gravitational dynamics | Metric action, its coupling scale and a consistent compactification vacuum | A local anomaly-related Pontryagin coupling; no derivation of Einstein dynamics or $M_P$ |

The exact charge algebra determines which interactions can occur. To turn it into a fermion hierarchy, the parent must determine their coefficients and a stable vacuum. To turn the internal radius into a physical mass scale, it must determine the radius potential and kinetic normalization. The additional $U(1)$ fields can participate in these processes, but their spin-one degrees of freedom do not supply the spin-two gravitational field simply by mixing.

The next investigation should therefore require the following concrete inputs and tests:

1. **Parent origin:** derive the two axions, their periods, kinetic matrix and complete anomaly action from a specified parent field content. Check the parent's own local and global anomalies.
2. **Scalar action:** derive the radius/saxion kinetic terms and potential. If supersymmetry is imposed, provide its Kähler potential, superpotential, moment maps and breaking assumptions.
3. **Interaction coefficients:** calculate whether $B_{ij}$ or its supersymmetric counterpart is generated, its flavor structure, and the coefficients of the other operators permitted by the same charges.
4. **Vacuum gate:** solve all stationarity conditions modulo gauge transformations; test the physical scalar Hessian, flux constraints and scale separation. A radial minimum alone is insufficient.
5. **Mass gate:** canonically normalize the surviving fields and compute masses from that vacuum. Keep $B-L$ breaking and the common-phase gauge field explicit.

Until these data exist, an arbitrary potential fitted to a desired minimum would add assumptions without testing why MTFT selects them.

## 9. Reproduction and provenance

The bundle contains three independent-purpose scripts, their machine-readable results, `run_all.py`, the copied input module and its MIT license. Run from the extracted folder with Python, SymPy and NumPy:

```bash
python run_all.py
```

Expected total: **159 passed, zero failed**. The source snapshot lets the script compare its independent polynomial construction with the supplied package's own anomaly functions; no installed `mtft` package, network access or proprietary algebra system is needed. This is a targeted extension audit, not a rerun of the package's entire test suite.

The supplied archive SHA-256 is:

```text
2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5
```

The source snapshot hash, runtime versions and individual check counts are in `results_summary.json`. `manifest.json` records deliverable file hashes. The original package, old registers and prior reports are unchanged.

## 10. Short glossary and hand exercises

| Term | Meaning in this calculation |
|---|---|
| Anomaly polynomial | A formal characteristic-class polynomial whose descent encodes the local quantum gauge variation |
| Stückelberg field | A scalar shifting under a gauge transformation so its covariant derivative is gauge invariant; it can supply a vector's longitudinal mode |
| Axion dressing | A periodic scalar exponential multiplying a charged operator so the combined gauge charge is zero |
| Wess–Zumino term | Here an axion times a four-form whose variation cancels a fermion anomaly |
| Generalized Chern–Simons counterterm | A local term shifting the distribution of a consistent anomaly among gauge currents |
| Pontryagin density | The curvature four-form $\operatorname{tr}(\mathcal R\wedge\mathcal R)$, distinct from scalar curvature in the Einstein action |
| Radion | The four-dimensional scalar describing changes of internal size |
| Saxion | A real scalar partner of an axion in a supersymmetric chiral field |
| Smith normal form | Integer diagonalization revealing the primitive charge lattice and possible finite remnants |
| Canonical normalization | Rescaling/mixing fields so their kinetic matrix becomes the standard identity before interpreting masses |
| Wilson coefficient | The strength multiplying an allowed effective operator; permission by symmetry does not calculate its value |

For practice, reproduce these four short calculations without the scripts:

1. Multiply $K$ by $Y$ and $B-L$, and solve $Kq=0$ using $q_c,q_L,q_a$ as free parameters.
2. Expand the eight cubes in §3.1 and regroup them into $r_1Q_1+r_2Q_2$.
3. Form $KH_0^{-1}K^{\mathsf T}$ and solve its quadratic characteristic equation.
4. Add the charge vectors of $H_uH_d$ and $QQQL$, verify their axion dressings, then explain why the same construction fails for $\nu^c\nu^c$.

The algebraic answers are displayed in the preceding sections. The remaining physical questions require new parent-theory information, not a longer manipulation of the same matrices.
