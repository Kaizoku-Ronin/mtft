# MTFT v0.31.4: the half-form lift, vacuum assumptions, and gravity

Investigation date: 16 September 2026. This is an additive audit of the supplied candidate, not an edit of the release or a replacement for earlier registers.

**Main finding.** The new matrices generate a valid order-eight dihedral group on the section space. However, the extra compatibility required of a lift of the theta characteristic changes the $W_{11}$ lift by a factor of $i$. The resulting group is $Q_8$. Both actions have commutator $-I$; the new geometric check determines which action respects $S_0^{\otimes2}\simeq K_X$.

The anomaly extension and restored M2 normalization provenance check out. The assertion that all Higgs directions are protected by four-dimensional $N=1$ supersymmetry remains unestablished in the supplied archive. This distinction changes the recommended order of work on the vacuum.

| Claim | Audit result | Status |
|---|---|---|
| The displayed $A,B$ generate $D_8$ | Correct for the stated normalization | EXACT |
| The same $A,B$ preserve the theta-square identification | $A$ fails; $iA,B$ satisfy it | EXACT correction |
| The projective AL commutator is $-I$ | Confirmed with either normalization | EXACT |
| This answers the spacetime Pin/Lorentz question | It establishes an internal curve lift; a spacetime identification is additional work | Not established |
| Cubic Abelian anomaly vanishes on span{phase,$Y$,$B-L$} | Entire symbolic polynomial is zero; mixed and gravitational coefficients also vanish | EXACT for the declared spectrum |
| A bare $\nu^c\nu^c$ term respects unbroken $B-L$ | Its charge is $+2$, so it is forbidden | EXACT selection rule |
| Raw M2 data reproduce the corrected normalized tensors | Both up and down reconstructed with zero numerical difference in this run | DIAGNOSTIC |
| $N=1$ makes every Higgs direction an exact modulus | Needs a parent, preserved supercharges, and explicit F/D-flatness conditions | Unestablished from supplied materials |
| Adding an Einstein term fixes the compactification radius | The minimal positive-energy extension has a radius runaway | EXACT within the restricted added ansatz |

## 1. What was inspected and reproduced

Input archive: `mtft-0.31.4.tar.gz`.

SHA-256:

```text
2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5
```

The archive contains 20 Markdown files under `docs/SM/`, including the normalization retraction and the response to the review. Important inspected files are:

- `src/mtft/surface/arithspin.py`: the added AL lift and the existing theta-characteristic machinery.
- `src/mtft/surface/smflux.py`: stack charges and the added anomaly/Majorana functions.
- `src/mtft/surface/hym.py`: corrected normalization and M2 data loading.
- `src/mtft/surface/_data/x0143_weight2_basis.json`: exact rational Atkin–Lehner matrices and q-expansions.
- `src/mtft/surface/_data/x0143_m2_tensors_h02.npz`: raw tensors, Gram matrices, and corrected normalized tensors.
- `docs/SM/SM01_REPORT.md`, `REVIEW_V0311_RESPONSE.md`, and the later SM/KK registers.

Seven relevant shipped tests passed: four in the review/arithspin modules and three in the M2 module. This audit did **not** rerun or independently certify Claude's reported full 810-test release suite. The separate audit program verifies the new algebra symbolically and records input hashes, intermediate matrices, and numerical diagnostics in `results.json`.

## 2. The section-space $D_8$ calculation is correct

Let

$$
X=X_0(143),\quad D=6(0)+6(1/11),\quad S_0=\mathcal O_X(D),
\qquad u=\frac{\eta(13\tau)\eta(143\tau)}{\eta(\tau)\eta(11\tau)}.
$$

The cusp order vector of $u$, in the order $0,1/11,1/13,\infty$, is

$$
\operatorname{div}(u)=(-6,-6,6,6).
$$

In the basis $(1,u)$ of $H^0(S_0)$, the section actions are

$$
A(s)=s\circ W_{11},\qquad
B(s)=\sqrt{13}\,u\,(s\circ W_{13}),
$$

$$
A=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad
B=\begin{pmatrix}0&-1/\sqrt{13}\\\sqrt{13}&0\end{pmatrix}.
$$

Their identities follow from

$$
u\circ W_{11}=-u,\qquad u\circ W_{13}=-\frac{1}{13u}:
$$

$$
A^2=I,\quad B^2=-I,\quad AB=-BA,\quad (AB)^2=I.
$$

The generated group has eight elements: one of order one, five of order two, and two of order four. Here $D_8$ means the dihedral group of order eight; this avoids the conflicting $D_4$ naming conventions from the earlier conversation.

The normalization matters. Without the factor $\sqrt{13}$, the raw pullback operator $B_0(s)=u(s\circ W_{13})$ has $B_0^2=-I/13$, not $-I$. The source docstring does disclose its normalization.

The audit independently derives both unit identities using exact rational Dedekind sums in the eta transformation formula. It does not promote a floating-point `allclose` test into an exact proof. The transformation formula is given by [NIST DLMF §23.18, equations 5–7](https://dlmf.nist.gov/23.18).

## 3. The theta-square condition selects $Q_8$

A theta characteristic comprises a line bundle together with an isomorphism

$$
\Phi:S_0^{\otimes2}\longrightarrow K_X.
$$

To lift an automorphism $w$ as a symmetry of this spin structure, its action $L_w$ must satisfy

$$
\Phi(L_ws\otimes L_wt)=w^*\Phi(s\otimes t).
\tag{1}
$$

Acting on the line bundle alone does not check (1).

For this particular $S_0$, the isomorphism is explicit. Set

$$
f(\tau)=\eta(\tau)^2\eta(11\tau)^2,\qquad
\omega=f(\tau)\,d\tau.
$$

The modular form has cusp orders $(13,13,1,1)$, so the differential has divisor

$$
\operatorname{div}(\omega)=(12,12,0,0)=2D.
$$

Thus

$$
\Phi(s\otimes t)=st\,\omega.
$$

Exact eta transformation laws give

$$
W_{11}^*\omega=-\omega,\qquad
W_{13}^*\omega=13u^2\omega.
\tag{2}
$$

These identities were also checked against the shipped rational $13\times13$ Atkin–Lehner matrices. The differential is the first frozen basis vector. Its q-expansion and the transformed products agree exactly through $q^{130}$; the weight-two Sturm bound at level 143 is 28. The eta calculation separately supplies the global identities and divisor information.

Equation (2) exposes the missing phase. For $A$, the left side of (1) is $(st)\circ W_{11}\,\omega$, while the right side is its negative. For $B$, the two sides already agree. Consequently take

$$
\widetilde A=iA,\qquad \widetilde B=B.
$$

Then

$$
\boxed{\widetilde A^2=\widetilde B^2=(\widetilde A\widetilde B)^2=-I,
\qquad [\widetilde A,\widetilde B]=-I.}
$$

The group is $Q_8$: the identity, the central involution $-I$, and six elements of order four. These are the full two sign choices above each element of the AL Klein group when compatibility with $\Phi$ is required. Multiplying $\Phi$ by a nonzero constant does not change condition (1).

Another exact way to see the correction is to use the product basis $(\omega,u\omega,u^2\omega)$ in $H^0(K_X)$. The natural actions are

$$
W_{11}=\operatorname{diag}(-1,1,-1),\qquad
W_{13}=\begin{pmatrix}0&0&1/13\\0&-1&0\\13&0&0\end{pmatrix}.
$$

These are $\operatorname{Sym}^2(\widetilde A)$ and $\operatorname{Sym}^2(B)$. They are **not** $\operatorname{Sym}^2(A)$ and $\operatorname{Sym}^2(B)$. The supplied audit asserts this intertwining identity explicitly.

Recommended additive correction: retain the valid $D_8$ result as the specified line-bundle normalization; add a separately named theta-compatible lift returning $iA,B$, and a test of (1). Replace the unqualified phrase “D8, not Q8, answers the Pin-lift question” with the distinction established here.

This does not automatically retract the existing Yukawa selection rules. Their section-space phase conventions and product maps must be compared consistently; the present correction concerns the additional canonical-square condition.

## 4. What survives, and what this says about triality

The robust shared datum is

$$
\kappa=ABA^{-1}B^{-1}=-I.
$$

It is unchanged by multiplying either generator by a nonzero scalar. The two groups induce the same projective action on $\mathbb P H^0(S_0)$. Their extension signs differ when the allowed scalar kernel is restricted to $\{\pm I\}$; specifying the canonical square map resolves that ambiguity for the spin structure.

The geometry gives a further constraint:

| Involution | Trace on $H^0(K_X)$ | Quotient genus | Fixed points |
|---|---:|---:|---:|
| $W_{11}$ | 1 | 7 | 0 |
| $W_{13}$ | -1 | 6 | 4 |
| $W_{143}$ | -9 | 2 | 20 |

For genus 13, the invariant differential space has dimension $g'=(13+\operatorname{tr}w)/2$. Riemann–Hurwitz then gives $r=28-4g'=2-2\operatorname{tr}w$. These are exact consequences of the frozen matrices, not approximate CM-point counts.

At a fixed point of a nonidentity holomorphic involution, the derivative is $-1$. A compatible square root of its cotangent action therefore squares to $-1$. This independently explains why the $W_{13}$ and $W_{143}$ spin lifts have order four.

Abstractly, $\operatorname{Out}(Q_8)\simeq S_3$, acting on its three noncentral pairs. But an automorphism of $X$ conjugating the AL involutions must preserve their numbers of fixed points. Since $0,4,20$ differ, no automorphism normalizing this $V_4\subset\operatorname{Aut}(X)$ can nontrivially permute its three involutions. Thus this proposed geometric realization of the abstract permutation symmetry fails an exact test. This statement concerns automorphisms of this fixed curve; it does not rule out a separately defined duality or correspondence on a larger construction.

The earlier TRI-01 obstruction also remains: the centered exterior-square construction from the quartic Hecke factor has cubic trace zero, whereas the observed sextic factor has cubic trace $-6$. The present internal spin lift does not supply a missing Hecke intertwiner, identify $8_v\oplus8_s\oplus8_c$ with the MTFT tensors, or identify AL involutions with spacetime parity/time reversal. In quantum theory time reversal also requires attention to antiunitarity.

## 5. Anomalies, Majorana masses, and normalization

For the supplied M1 stack spectrum, define

$$
q_x=\alpha+\beta Y_x+\gamma(B-L)_x,
\quad
\mathcal A(q)=\sum_{x<y}(m_x-m_y)N_xN_y(q_x-q_y)^3.
$$

Symbolic expansion gives $\mathcal A(q)=0$ as a polynomial in all three variables. The mixed non-Abelian coefficients and the linear gravitational coefficient also vanish. Polarization therefore establishes all mixed Abelian cubic coefficients within this span, rather than just cancellation at a few sampled charge vectors. The supplied mixed-anomaly helper omits the conventional overall fundamental Dynkin-index factor $1/2$; that convention does not affect any vanishing claim.

The common phase acts trivially on bifundamentals: adding the same number to every $q_x$ leaves every charge difference unchanged. The three-dimensional stack-charge span consequently contains two nontrivial matter-charge directions, represented by $Y$ and $B-L$.

For left-handed $\nu^c$, $B-L=+1$, hence the bilinear has charge $+2$. A bare Majorana mass is forbidden while this symmetry is unbroken. A renormalizable coupling to a scalar of charge $-2$ is one possible completion; another symmetry-breaking mechanism could generate an effective operator. Its scale is not derived by this selection rule. These are four-dimensional anomaly/charge checks of the chosen spectrum, not a six-dimensional quantum-consistency calculation.

M2 reconstruction used the supplied raw up/down tensors and all required Gram matrices with `hym.normalise_yukawa`. Both reconstructed normalized arrays agree with the frozen arrays with relative difference 0.0 in this runtime. Both family-basis invariance checks pass. With 10,000 isotropic complex Gaussian Higgs directions per sector and seed 20260916:

| Sector | Median $m_1/m_3$ | Median $m_2/m_3$ |
|---|---:|---:|
| M2 up | 0.19051 | 0.55760 |
| M2 down | 0.26073 | 0.61917 |

These confirm the corrected broad distribution and dependence on a Higgs-direction choice. Finite-sample variations explain small differences from earlier quoted medians. They are numerical diagnostics, not observational predictions or rigorous error bounds. The withdrawn fit vectors are absent from the frozen data.

Some historical register text still presents the retracted hierarchy results, and the `load_m2_tensors` docstring still mentions fit vectors. Preserve the historical entries, but add prominent supersession notices pointing to CC-26 and update the loader description. The phrase “corrected everywhere” is stronger than the shipped prose supports.

## 6. The supersymmetry assertion needs a parent and a flatness calculation

The response invokes CW-01, but no CW-01 derivation is included in the supplied source/archive registers. More concretely, SM-01's recorded caveat distinguishes the implemented fermion index from a proposed supersymmetric parent. The index identities themselves are

$$
\chi(S_0\otimes E)=\deg E,\qquad
\chi(E)=\deg E+1-g=\deg E-12.
$$

This is not a theorem forbidding a supersymmetric completion. Different multiplets, R-charges, and twists can produce different bundle assignments. It does mean that the actual parent and zero-mode bundles must be supplied before transferring the three-family result to its chiral multiplets.

Likewise a holomorphic section of $S_0$, or a Hermitian–Yang–Mills metric on a flux line bundle, is not by itself a preserved supercharge. The curved compactification requires appropriate Killing-spinor/twist conditions and vanishing fermion supersymmetry variations. A concrete six-to-four-dimensional construction exhibits this extra step explicitly in [Razamat, Vafa and Zafrir, §2 and §7](https://arxiv.org/pdf/1610.09178).

Even after $N=1$ is established, it does not imply every scalar direction is flat. For example, with canonical kinetic terms and

$$
W=\mu H_u\cdot H_d,\qquad \mu\ne0,
$$

an unbroken supersymmetric vacuum can sit at the origin while

$$
V_F=|\mu|^2\left(|H_u|^2+|H_d|^2\right)
$$

gives both Higgs fields a potential. This is a counterexample to the general implication, not an assertion that MTFT contains this particular term. In rigid $N=1$, the full scalar potential depends on the superpotential, Kähler metric and D-terms; a moduli branch must satisfy all F- and D-flatness conditions modulo gauge equivalence. See [Martin, A Supersymmetry Primer, §§3.4 and 4.10](https://arxiv.org/pdf/hep-ph/9709356).

Once an actual branch of supersymmetric vacua is established, perturbative protection can indeed apply. The correction is to establish that branch and its assumptions first. Nonperturbative effects and a supergravity completion require their own analysis. A statement of perturbative protection alone does not fix the radius or establish a cosmological constant.

The next parent/vacuum calculation should therefore specify the six-dimensional multiplets and chirality, the R-symmetry background or twist, flux and auxiliary-field equations, and the surviving four-dimensional bundles. Then write $W,K,D_a$, solve the vacuum equations, and identify the surviving Higgs directions. Only after that calculation can one decide what breaking or additional dynamics is needed.

## 7. What adding gravity would actually require

The nested-spaces image is a hierarchy of mathematical structures and induced topologies. It is not a hierarchy of physical scales. A Lorentzian metric is nondegenerate but indefinite; it does not induce a positive norm or a distance by the usual inner-product construction. MTFT's positive Hilbert/Gram inner products can coexist with an external Lorentzian spacetime metric because they act on different objects.

For a gravity model one must specify a spacetime manifold, a dynamical metric, fields/bundles, an action, and a background satisfying the field equations. The modular curve can be an internal factor of such a construction; its fixed metric does not itself supply the spacetime Einstein equations.

For an **assumed**, unwarped product with constant radius $R$, write the six-dimensional Einstein term as $M_6^4\mathcal R_6/2$. If the compact genus-13 curve has Gaussian curvature $-1/R^2$, Gauss–Bonnet gives

$$
\mathcal A_X=48\pi R^2,\qquad
M_4^2=M_6^4\mathcal A_X,
\qquad \frac1{g_4^2}=\frac{\mathcal A_X}{g_6^2},
\qquad m_{\mathrm{KK},n}^2=\frac{\lambda_n}{R^2}.
$$

The Planck and gauge relations follow by integrating the corresponding kinetic terms over the internal factor. The KK relation refers to the relevant dimensionless internal operator without additional mass terms. None fixes $R$, $M_6$, or a mass in GeV. For a varying radius the four-dimensional metric must be rescaled to Einstein frame.

There is already a useful conditional exclusion test. Take only two-derivative Einstein gravity, nonnegative bulk vacuum energy, and ordinary positive gauge kinetic terms with fixed quantized internal flux. Assume no warping, localized sources, or other corrections. In four-dimensional Einstein frame, with $\rho=R/R_0$, the radius potential has the form

$$
V_E(\rho)=\frac a{\rho^4}+\frac b{\rho^6}+\frac c{\rho^2},
\qquad a>0,\ b\ge0,\ c\ge0.
$$

Here negative internal curvature produces the positive $a$ term, flux energy the $b$ term, and bulk vacuum energy the $c$ term. Its derivative is

$$
\frac{dV_E}{d\rho}=-\frac{4a}{\rho^5}-\frac{6b}{\rho^7}-\frac{2c}{\rho^3}<0.
$$

Thus this restricted extension has no stationary finite radius. This is a constraint on a candidate gravity parent, not a no-go theorem for all MTFT completions. Stabilization requires a contribution or background outside these assumptions. The next useful gravity experiment is to derive the full allowed radius/Higgs potential from a specified parent and test its extrema and stability.

## 8. Recommended next work and glossary

1. Append the theta-compatible $Q_8$ result and its square-map test to the AL lift record. Preserve the existing $D_8$ normalization as a distinct valid action.
2. Supply or reconstruct CW-01's actual parent, twist and F/D-flatness calculation. Verify that its spectrum is the one used in the numerical Yukawa construction.
3. Use that parent to derive the coupled radius/Higgs potential and six-dimensional anomaly conditions. Test whether the compactification is a solution before treating its radius or Higgs direction as a physical input to a prediction.

| Term | Meaning in this audit |
|---|---|
| Theta characteristic | A line bundle $S$ with a chosen isomorphism $S^2\simeq K_X$ |
| Half-form | A section of such a square root of the canonical bundle |
| Section lift | An action covering a curve automorphism on a line bundle or its sections |
| Theta-compatible lift | A section lift whose tensor square agrees with natural pullback of differentials |
| Projective action | A linear action considered up to nonzero scalar factors |
| $D_8$ | The dihedral group of order eight, with two elements of order four |
| $Q_8$ | The quaternion group, with six elements of order four |
| Modulus | A coordinate on a family of physically inequivalent vacua; an unspecified parameter is not automatically one |
| F/D-flatness | Vanishing auxiliary-field conditions for a supersymmetric vacuum in the stated rigid theory |
| Radion | The four-dimensional scalar describing variation of the internal size |
| EXACT | An identity proved algebraically under stated mathematical/model assumptions |
| DIAGNOSTIC | Numerical evidence without a rigorous total error bound |

The concrete new result is the compatible $Q_8\to V_4$ lift and its exact product-map verification. The next physical question is which parent equations make the chosen curve, fluxes, radius, and Higgs configuration a vacuum.
