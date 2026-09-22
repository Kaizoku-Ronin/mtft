# Hodge dynamics and explicit period transport

Research continuation · 19 September 2026

## Result

The proposed “gauge realm” now has an explicit mathematical realization at the level of connections and monodromy. The elliptic family has a verified rational period connection. The next curve has an exact connection that splits into inherited elliptic and new genus-two components. Closed parameter paths can generate exponential transport growth, and that growth survives changes of frame.

Two exact obstructions also emerged. Chaotic parameter motion confined to a compact nonsingular real interval produces zero asymptotic growth under ordinary Gauss–Manin transport. Separately, holomorphic pullback along the tower, or along elliptic isogenies, is an isometry after dividing by the square root of its degree. Neither operation alone supplies Feigenbaum scaling.

This is a mathematical bundle and transport construction over parameter space. It does not yet specify a physical gauge field, a spacetime action, or an MTFT-selected parameter evolution.

## 1. The exact elliptic connection

Start with the previous investigation’s first elliptic member,

\[
C_{2,c}:y^2=P(x,c)=x^4+2cx^2+c^2+c,
\qquad c\in B=\mathbb C\setminus\{0,-1\}.
\]

Use the residue-free de Rham basis

\[
\omega_0=dx/y,\qquad \omega_2=x^2dx/y.
\]

The first differential is holomorphic; the second is of the second kind. An arbitrary list of powers times dx/y on an even-degree hyperelliptic model is not automatically a compact-curve de Rham basis: residues at the two points at infinity must be checked.

For a flatly transported cycle γ, let Π=(∮γω₀,∮γω₂)ᵗ. Then

\[
\frac{d\Pi}{dc}=A_E(c)\Pi,
\qquad
\boxed{A_E(c)=
\begin{pmatrix}
-1/(4c)&1/[4c(c+1)]\\
-1/4&1/(4c)
\end{pmatrix}.}
\]

This is an exact identity in cohomology. The discarded exact differentials have primitives R₀/y and R₂/y, where

\[
R_0=-\frac{x(x^2+3c+1)}{4c(c+1)},
\qquad R_2=\frac{x}{4}-\frac{x^3}{4c}.
\]

The script `root_checks/elliptic_exact.py` clears denominators and verifies both polynomial identities over Q[c,x], without specializing c.

The period equation defines a flat connection ∇=d−A_Edc. Under the frame change Π̃=G(c)Π, its matrix becomes

\[
\widetilde A_E=G'G^{-1}+GA_EG^{-1}.
\]

For a closed loop based at c*, monodromy changes by conjugation, M̃=G(c*)MG(c*)⁻¹. Its trace and eigenvalues therefore survive a change of frame. These are the gauge-geometric features being used here. The full flat connection and the generally curved connection induced on the holomorphic Hodge subbundle are different objects.[1]

Eliminating the second period gives a scalar Picard–Fuchs equation:

\[
I''+\left(\frac1c+\frac1{c+1}\right)I'
+\frac{4c-1}{16c^2(c+1)}I=0.
\]

Its regular singular points are 0, −1, and infinity. The local exponents are respectively {−1/4,1/4}, {0,0}, and {1/2,1/2} in the local coordinate 1/c at infinity.

## 2. An exact obstruction: a chaotic base can have zero transport growth

Let K be a compact real interval avoiding the singular parameters, and Y(c) an invertible fundamental matrix of the period equation on an interval containing K. For any base map b:K→K, direct real-interval transport is

\[
T(c)=Y(b(c))Y(c)^{-1}.
\]

Products telescope:

\[
T(b^{n-1}c)\cdots T(c)=Y(b^n(c))Y(c)^{-1}.
\]

Both Y and its inverse are bounded on K. Consequently all singular values of the product stay bounded above and away from zero, uniformly in n. Every transport Lyapunov exponent is zero. The same conclusion holds in the Hodge norm, since continuous positive-definite norms are uniformly equivalent on this compact nonsingular set.

This proof does not assume the base dynamics is regular. It applies equally to a chaotic base. In dynamical terminology the transport is a bounded coboundary.

**Chosen numerical control.** We deliberately imposed the fully chaotic logistic map on K=[−3,−2]:

\[
x=c+3,\qquad x_{n+1}=4x_n(1-x_n),\qquad c_n=x_n-3.
\]

This parameter rule was chosen to test the obstruction; it is not derived from MTFT. Three initial conditions were followed for 8192 iterations, with two ODE tolerances. The base finite exponents were 0.69315809, 0.69313784, and 0.69313264, close to the usual almost-everywhere value log 2. The transport finite exponents were approximately 1.19×10⁻⁵, 4.81×10⁻⁶, and 7.75×10⁻⁶. Their exact asymptotic limit is zero by the telescoping proof, not by extrapolation of those finite values.

Direct matrix products and the telescoped expression agreed within 1.9×10⁻¹⁴ at the final checkpoints. These numbers are floating-point diagnostics. The logistic trajectories were computed in float64 and are not certified trajectories of specified exact real initial conditions.

## 3. Winding around singularities produces nontrivial growth

The compact-interval argument does not apply to path protocols with nontrivial winding in B. A global single-valued Y is then unavailable.

Take standard counterclockwise based loops S around c=0 and U around c=−1. The code specifies the paths and uses c*=1. The scalar equation and its local exponents give

\[
S^2=-I,\qquad U=I+N,\quad N^2=0,\ N\ne0,
\qquad \operatorname{tr}(US)=-2.
\]

The last identity follows from the loop relation at infinity. Since tr A_E=0, determinants are one. Therefore

\[
\operatorname{tr}(U^kS)=-2k.
\]

This determines several exact growth regimes:

| Repeated loop protocol | Asymptotic operator-norm growth |
|---|---|
| Sⁿ | Bounded; S⁴=I |
| Uⁿ | Linear; Uⁿ=I+nN |
| (US)ⁿ | Linear; a nontrivial parabolic matrix up to sign |
| (U²S)ⁿ | Exponential, with factor 2+√3 per repetition |

For the final word,

\[
\det(\lambda I-U^2S)=\lambda^2+4\lambda+1,
\qquad
\rho(U^2S)=2+\sqrt3,
\]

and the rate per repeated word is log(2+√3)≈1.3169578969. More generally, k≥2 gives spectral radius k+√(k²−1). Thus the protocol must be specified; the geometry does not select a unique growth constant.

Numerical monodromy integration found ρ(U²S)=3.732050807568958, maximum determinant defect below 5.5×10⁻¹⁵, and a nonconstant-frame conjugacy defect of 3.6×10⁻¹⁴. Three solver tolerances were compared. These are consistency diagnostics, not interval bounds. Defective parabolic eigenvalues show apparent splitting near 10⁻⁷ from errors near 10⁻¹⁴; the exact trace relations determine the correct classification.

The repeated word is periodic parameter driving. It creates exponential linear transport without requiring chaotic driving. Its rate is not a physical rate until a traversal-time convention is supplied. No identification with Feigenbaum δ is established. The computed matrices are in a complex period frame; an integral cycle basis was not reconstructed.

## 4. The inherited and new components remain separated under transport

At the next level,

\[
C_{3,c}:y^2=f_c^{\circ3}(x),\qquad f_c(x)=x^2+c,
\]

the genus is three. The two quotient curves are C₂ and

\[
D_{2,c}:v^2=Q(u,c)=(u-c)((u^2+c)^2+c),
\]

which has genus two. The maps use u=x²+c, with y unchanged for the C₂ quotient and v=xy for the D₂ quotient.

For D₂, the compact-curve de Rham basis is ηⱼ=uʲdu/v, j=0,1,2,3. The first two forms are holomorphic. The rational 4×4 connection A_D(c) was solved exactly over Q(c). Every row is certified by an identity

\[
-\tfrac12 u^j Q_c
=\left(\sum_{k=0}^3 A_{D,jk}u^k\right)Q
+R_j'Q-\tfrac12R_jQ_u.
\]

All coefficients are included in the JSON output, together with the exact primitives Rⱼ. Its possible finite poles lie at c=0,−1 and the roots of c³+2c²+c+1, precisely among the critical-return singular parameters for these models.

The six pulled-back forms give the exact block connection

\[
\boxed{A_{C_3}=A_E\oplus A_D.}
\]

This was checked directly on C₃, including the c-dependence of the quotient maps and basis. It is not inferred only from dimensions. In that basis the reflection has matrix diag(I₂,−I₄), and its two projectors commute with the connection.

The general reason is functoriality: the fiberwise reflection is a morphism of the family and commutes with Gauss–Manin transport. Thus the inherited/new splitting is a splitting of local systems as well as a fiberwise Hodge decomposition. Ordinary parameter transport cannot mix these two sectors. A mixing interaction would require additional structure that does not preserve the same projectors.

**A loop that distinguishes the sectors.** Let a≈−1.754877666 be the real root of c³+2c²+c+1. This is a primitive period-three critical-return parameter. The inherited C₂ remains smooth there, whereas the new D₂ degenerates. We integrated a counterclockwise loop of radius 0.05 around a, with a common basepoint c=−2.

The inherited monodromy was identity to a norm defect of 7.9×10⁻¹⁶. The new monodromy was nontrivial: M_D−I had singular values approximately 19.13, 4.3×10⁻¹⁵, 3.0×10⁻¹⁵, and 3.8×10⁻¹⁶, while ‖(M_D−I)²‖ was 2.5×10⁻¹³. These are numerical signatures of the expected rank-one unipotent response to the node. The large first singular value depends on the frame; the distinction between identity and nontrivial unipotent monodromy does not.

The exact classification comes from the local geometry: the inherited connection extends across the disk, while the new curve develops a nonseparating node. Thus M_E=I and M_D=I+N with N≠0, N²=0, and rank N=1. The companion `monodromy/NEW_COMPONENT_NOTE.md` gives the local equation and explains the quadratic smoothing; numerical singular values are supporting evidence rather than the rank proof.

The D₂ discriminant contains the factor (c³+2c²+c+1)²: its local node smoothing has quadratic parameter dependence. The related C₃ model has a transverse node at the same critical return. Integral cycle normalizations should therefore not be identified across their isogeny without further work. This experiment tests a genuinely new sector; it does not test a period-doubling parameter, whose next primitive period is four.

## 5. Isogeny scaling is controlled by degree

For a degree-d finite holomorphic map φ:C→D and a cohomology class v on D,

\[
\|\phi^*v\|_{H,C}^2=d\|v\|_{H,D}^2.
\]

Holomorphic pullback preserves harmonicity of one-forms, and integration of their wedge product multiplies by d. Therefore d⁻¹ᐟ²φ* is an isometric embedding. It is an isometry for elliptic isogenies, where the dimensions agree.

For the degree-two tower maps, the inherited Hodge directions have exactly this normalization. Their raw norm factor √2 records covering degree; normalized iteration cannot supply an expanding Feigenbaum mode.

The elliptic family also admits the previously found Fricke parameter involution

\[
\iota(c)=-\frac{c}{c+1},\qquad \iota^2(c)=c.
\]

An explicit lift requires a square-root choice. The quartic is birational to

\[
W_c:Y^2=X^3-4cX^2-4cX,
\quad X=2(y+x^2+c),\quad Y=2xX.
\]

Here dX/Y=dx/y. Set s²=2(c+1) and c′=ι(c). The standard two-isogeny followed by scaling is

\[
X'=\frac{X-4c-4c/X}{s^2},\qquad
Y'=\frac{Y(1+4c/X^2)}{s^3}.
\]

It maps W_c to W_c′ after adjoining s. On the quartic cohomology basis its pullback matrix is

\[
B(c)=\begin{pmatrix}s&0\\2c/s&2/s\end{pmatrix},
\qquad\det B=2.
\]

With compatible choices s(c′)s(c)=2, the round trip satisfies B(c′)B(c)=2I. Reversing a square-root branch can introduce a central minus sign. Degree normalization gives the identity, or its negative, on the round trip. This lift is therefore explicit but does not supply an expanding parameter renormalization. The underlying equality of j-invariants alone would not have fixed these lift choices.

The lift also intertwines the connection exactly:

\[
\frac{dB}{dc}+BA_E(c)=\iota'(c)A_E(\iota(c))B.
\]

Writing B=sD removes the square root from this verification. The note `dynamics/DYNAMICS_AND_FRICKE.md` displays both sides as the same rational matrix; its script supplies 80 supporting exact rational specialization checks.

## 6. A scaling ratio in period samples can come from the input

At every fixed tower level, the Feigenbaum accumulation parameter c∞ is smooth. A locally chosen period P(c) is therefore analytic near that point. Suppose an already chosen parameter sequence obeys

\[
c_n-c_\infty\sim K\delta^{-n},\qquad K\ne0.
\]

If P′(c∞)≠0, Taylor expansion immediately gives

\[
\frac{P(c_{n-1})-P(c_n)}{P(c_n)-P(c_{n+1})}\longrightarrow\delta.
\]

If the first nonzero derivative has order k, the limit instead is δᵏ. Consequently, finding the Feigenbaum ratio in smooth period samples at preselected Feigenbaum parameters would not independently derive it. The ratio can be inherited from the input sequence. Cumulative local transport along that sequence converges to a finite endpoint transport.

## 7. What now qualifies as the proposed framework

| Ingredient | Status |
|---|---|
| Hodge/cohomology bundle over a punctured parameter space | Established construction |
| Explicit elliptic and new genus-two connections | EXACT algebraic reduction |
| Gauge covariance and frame-independent monodromy characters | EXACT; numerical checks included |
| Inherited/new projectors preserved by transport | EXACT |
| Exponential growth under specified winding protocols | EXACT spectral calculation; numerical monodromy support |
| Positive transport exponent from arbitrary chaotic interval forcing | Ruled out under the stated bounded-flat-transport hypotheses |
| Expansion from degree-normalized holomorphic pullback | Ruled out for the inherited directions |
| Feigenbaum renormalization acting on this tower | Open; no such operator constructed here |
| MTFT selects the base evolution or a physical time | Not established |
| A physical gauge field or gravity model | Not established by these calculations |

Kontsevich–Zorich theory remains the relevant established comparison: Teichmüller flow supplies the base motion, Gauss–Manin supplies transport, and the Hodge norm measures growth. Its Lyapunov spectrum has proven relations to Hodge-bundle geometry.[1,2] Applying that mechanism here would require specifying the relevant flow and showing that the chosen curve/differential locus is preserved, or enlarging the family. It is not enough to relabel logistic iteration as a flow on moduli space.

The next genuinely different target is a renormalization correspondence that survives the interval and degree-normalization obstructions. It must specify how it acts on parameters, curves, cohomology, and norms, and why it is selected. The present work supplies explicit connections and projectors against which such a proposal can be tested.

## Reproduction and scope

Read `README.md` for commands. Exact reduction scripts use Python’s standard library. Numerical ODE experiments use NumPy and SciPy; the figure uses Matplotlib. No repository files were changed. This bundle extends the preceding triangle/Feigenbaum/Hodge report; it does not rerun or certify the full MTFT suite.

EXACT labels apply to algebraic identities or proofs with stated hypotheses. Floating-point solver agreement, determinants, and finite growth rates are DIAGNOSTIC. No claim of historical novelty is made.

## References

1. Forni, Matheus and Zorich, *Lyapunov spectrum of invariant subbundles of the Hodge bundle*, especially §§2–3: [author-hosted paper](https://webusers.imj-prg.fr/~anton.zorich/Papers/MatheusForniZorich%20Lyapunov%20Spectrum%20ETDS%202012.pdf).
2. Eskin, Kontsevich and Zorich, *Sum of Lyapunov exponents of the Hodge bundle with respect to the Teichmüller geodesic flow*: [research paper](https://arxiv.org/abs/1112.5872).

The rational reductions, compact-interval proof, and monodromy character calculation are given explicitly above and in the accompanying programs. The references establish the surrounding Hodge-dynamical framework; they are not being cited as proofs that our family realizes Feigenbaum renormalization.
