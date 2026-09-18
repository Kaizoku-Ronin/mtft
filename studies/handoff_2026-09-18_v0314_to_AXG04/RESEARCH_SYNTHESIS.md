# From v0.31.4 to an explicit charged parent candidate

**Synthesis date:** 18 September 2026. This document summarizes the included studies; their reports contain the complete proofs, conventions, references and machine ledgers. New module names elsewhere in the handoff are implementation proposals.

## 1. Where the research now stands

The work moved from identifying internal arithmetic symmetries to testing explicit higher-dimensional field theories. It produced reusable exact tools, several sharp obstructions, and a concrete local six-dimensional EFT candidate. It did not derive a unique parent theory, realistic gravitational vacuum, observed flavor hierarchy, or absolute mass scale.

Three different notions must remain explicit throughout:

| Object | What its mathematics controls | What is not implied |
|---|---|---|
| Modular curve and its bundles | Divisors, cohomology, indices, Hecke actions and internal mode spaces | A physical field equation or spacetime metric |
| Chosen parent and mode operator | Which fields propagate, which sections are zero modes, allowed interactions and anomalies | That the chosen background solves the equations |
| Action plus a vacuum | Kinetic normalization, masses, curvature and stability within that action | Uniqueness, UV completion or agreement with observations |

The positive Hodge/Gram metric acts on internal state spaces. A Lorentzian spacetime metric is indefinite and does not produce a positive norm by the same construction. The nested norm/metric/topology diagram is a hierarchy of mathematical structures, not a derivation of spacetime dynamics.

## 2. Baseline, conventions and the first correction

The fixed input is `baseline/mtft-0.31.4.tar.gz`, SHA-256

```text
2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5
```

Its relevant modules are `mtft.surface.arithspin`, `smflux` and `hym`. The supplied candidate includes the corrected CC-26 contraction, raw M2 tensors and Gram matrices, the extended anomaly check and Majorana obstruction. The initial audit independently reconstructs M2 normalization and records numerical agreement; it does not turn finite-element or sampling results into error-certified predictions.

Let

\[
X=X_0(143),\quad g=13,\quad D=6c_0+6c_{11},\quad S_0=\mathcal O_X(D),
\quad u=\frac{\eta(13\tau)\eta(143\tau)}{\eta(\tau)\eta(11\tau)}.
\]

In cusp order \((0,1/11,1/13,\infty)\),

\[
\operatorname{div}(u)=(-6,-6,6,6),\qquad H^0(S_0)=\operatorname{span}(1,u).
\]

The section-space matrices are

\[
A=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad
B=\begin{pmatrix}0&-1/\sqrt{13}\\\sqrt{13}&0\end{pmatrix}.
\]

They satisfy \(A^2=I\), \(B^2=-I\), \(AB=-BA\), \((AB)^2=I\): this is **D8, the dihedral group of order eight**. The unnormalized pullback \(B_0=uW_{13}^*\) instead has \(B_0^2=-I/13\).

The spin object includes a specified square map \(\Phi:S_0^2\to K_X\), not just the line bundle. Set

\[
\omega=\eta(\tau)^2\eta(11\tau)^2d\tau,\qquad
\Phi(s\otimes t)=st\omega.
\]

Exact eta transformations give

\[
W_{11}^*\omega=-\omega,\qquad W_{13}^*\omega=13u^2\omega.
\]

Consequently the compatible lift is \(\widetilde A=iA\), \(\widetilde B=B\), with

\[
\widetilde A^2=\widetilde B^2=(\widetilde A\widetilde B)^2=-I,
\qquad[\widetilde A,\widetilde B]=-I.
\]

This is **Q8**. On the product basis \((\omega,u\omega,u^2\omega)\), the natural matrices are

\[
W_{11}=\operatorname{diag}(-1,1,-1),\qquad
W_{13}=\begin{pmatrix}0&0&1/13\\0&-1&0\\13&0&0\end{pmatrix}.
\]

They equal \(\operatorname{Sym}^2(\widetilde A)\) and \(\operatorname{Sym}^2(B)\). This intertwining identity is the decisive regression test. The projective commutator survives both normalizations; the actual central extension requires the square-map convention.

The involutions have fixed-point counts \(0,4,20\) for \(W_{11},W_{13},W_{143}\). Therefore a curve automorphism normalizing this AL group cannot permute its three nonidentity elements. Abstract \(\operatorname{Out}(Q_8)\simeq S_3\) does not provide that geometric action. Nor does this internal lift identify AL involutions with spacetime parity/time reversal or establish Spin(8) triality for the MTFT tensors.

**Source:** [initial audit §§2–6](studies/00_v0314_audit/MTFT_v0314_Spin_and_Vacuum_Audit.md).

## 3. The seven-dimensional spin-circle test

SC7-01 tests a specified smooth compactification

\[
M^{1,3}\times P,\qquad P=S(S_0)\xrightarrow{\pi}X,
\]

where the circle bundle has Euler degree \(e=12\). A standard commuting twisted 7D gauge mode complex is assumed, with no defects or boundaries. A provisional U(8) adjoint decomposition supplies the proposed bifundamental sectors; pure product-group SYM does not provide them automatically.

Its integral topology is

\[
H_1(P;\mathbb Z)=\mathbb Z^{26}\oplus\mathbb Z/12,\qquad
H^2(P;\mathbb Z)=\mathbb Z^{26}\oplus\mathbb Z/12.
\]

A degree-\(m\) line bundle pulls back to a torsion class of order
\(12/\gcd(12,m)\). For a circle connection \(\Theta\) with period \(2\pi\),

\[
\widehat A_m=\pi^*A_m-\frac m{12}\Theta,
\qquad F_{\widehat A_m}=0,
\qquad\lambda_m=e^{-2\pi im/12}.
\]

If \(\lambda_m\ne1\), the fiber cochain differential is multiplication by the invertible scalar \(\lambda_m-1\). The fiber-cohomology spectral sequence gives

\[
H^k(P;\mathcal L_m)=0\quad\text{for all }k.
\]

Thus the original degree-three and degree-six M1 sectors have **no massless bulk modes in this chosen 7D complex**. Their holonomies are respectively \(-i\) and \(-1\) in the displayed orientation. Degree-zero controls survive: the cohomology dimensions are \((1,26,26,1)\) for the trivial character and \((0,24,24,0)\) for a nontrivial base character. The finite exact scan checks 36 primary characters; the theorem covers arbitrary compatible base characters.

The connection metric has

\[
ds_P^2=ds_X^2+r^2\Theta^2,\qquad
\mathcal R_P=-\frac2{R^2}-\frac{r^2}{8R^4},\qquad
\operatorname{Vol}(P)=96\pi^2R^2r.
\]

It is not Einstein or a hyperbolic 3-manifold. Its separately assumed Einstein reduction has

\[
V_E(\rho)=a\rho^{-5}+b\rho^{-7}+c\rho^{-3},
\qquad a>0,\ b,c\ge0,
\]

so no stationary overall scale. Defects, boundaries, other mode operators and different parents require a new spectrum/anomaly calculation; they cannot inherit the old 6D family count unchanged.

**Source:** [SC7-01](studies/01_SC7_01/MTFT_SC7_01_Spin_Circle_Investigation.md).

## 4. Hopf geometry, Hecke multiplicities and the number 143

The actual level factorization is

\[
\mathbb Z/143\simeq\mathbb F_{11}\times\mathbb F_{13},\qquad
(a,b)\mapsto78a+66b\pmod{143}.
\]

It yields \(\mathbb P^1(\mathbb Z/143)\simeq\mathbb P^1(\mathbb F_{11})\times\mathbb P^1(\mathbb F_{13})\), a set of \(12\cdot14=168\) elements. The modular matrices

\[
S=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
T=\begin{pmatrix}1&1\\0&1\end{pmatrix}
\]

act on primitive rows; T has cusp cycles of lengths \(1,11,13,143\). The level-forgetting maps to \(X_0(11)\) and \(X_0(13)\) have generic degrees 14 and 12. Their moduli relation is over the common j-line, with normalization caveats for compact coarse curves. It is not a product of independent modular curves.

Standard \(S^{143}\) denotes a real sphere of dimension 143, not a Shimura variety assigned to level 143. The quaternionic Shimura curve sometimes denoted \(X^{143}\) is another complex curve, with 143 as a quaternion-algebra discriminant. Neither notation turns \(X_0(143)\) into a high-dimensional sphere.

There is a precise sphere statement in the **state space**:

\[
H^1(X;\mathbb R)=V_{12}\oplus V_{14},\qquad
S(H^1)=S(V_{12})*S(V_{14})\cong S^{11}*S^{13}=S^{25}.
\]

Here \(V_{12}\) is the sextic Hecke block and \(V_{14}\) its elliptic + old + quartic complement. The join follows from Hodge orthogonality; the factors are not assigned to the level primes individually. HOPF-02 records rational 26×26 complementary projectors. Its matrices use the package homology basis; transpose when acting on the dual cohomology basis.

The half-form pencil

\[
f:X\to\mathbb{CP}^1,\qquad p\mapsto[1:\sqrt{13}u(p)]
\]

is base-point-free of degree 12. Pulling the quaternionic Hopf connection back along the embedded \(\mathbb{CP}^1\subset\mathbb{HP}^1\) gives

\[
E\simeq S_0^{-1}\oplus S_0,\qquad\deg E_\pm=\mp12.
\]

This connection is reducible to U(1); the underlying SU(2) bundle on the curve is topologically trivial. Arithmetic left Q8 actions moving the Hopf base and right fiber actions fixing it must be distinguished. The map has total ramification 48, so the pulled-back round metric degenerates at branch points; the bundle connection remains smooth.

The quartic Hecke block has real dimension eight, but four distinct real T2 eigenvalues each have multiplicity two. A commuting quaternionic structure would have to act on each real plane. For

\[
I=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
IJ=-JI\Rightarrow J=\begin{pmatrix}a&b\\b&-a\end{pmatrix},
\quad J^2=(a^2+b^2)I_2,
\]

the required \(J^2=-I_2\) is impossible over the reals. Doubling multiplicity gives an explicit 16-dimensional positive control. It enlarges the state space and supplies a control metric; it does not identify the original modes or their Petersson metric with that control.

**Source:** [HOPF-02](studies/02_HOPF_02/MTFT_HOPF02_Report.md).

## 5. AXG-01: the two extra Abelian directions

M1 uses stack order \((c,L,a,b,d)\), ranks \((3,2,1,1,1)\) and flux
\(m=(0,-3,3,3,0)^T\). The primitive shift matrix is

\[
K=\begin{pmatrix}0&-2&1&1&0\\3&-4&0&0&1\end{pmatrix}
=\begin{pmatrix}k_1\\k_2\end{pmatrix}.
\]

Its kernel is spanned by common phase \(p\), hypercharge Y and B−L:

\[
p=(1,1,1,1,1)^T,\quad
Y=(1/6,0,-1/2,1/2,-1/2)^T,\quad
B-L=(1/3,0,0,0,-1)^T.
\]

The common phase acts trivially on these bifundamentals; this alone does not remove its gauge field from an action. The two rows are primitive (Smith factors 1,1).

For generator differences \(x=q_c-q_L, y=q_a-q_L,z=q_b-q_L,w=q_d-q_L\), set

\[
r_1=y+z,\quad r_2=3x+w,
\]
\[
Q_1=24(y^2-yz+z^2)+27x^2+9w^2,\qquad Q_2=-9(y^2+z^2).
\]

The complete cubic trace is \(P=r_1Q_1+r_2Q_2\). With normalized curvatures and characteristic-class conventions from the report,

\[
I_6=(Kf)_1X_4^1+(Kf)_2X_4^2,
\]
\[
X_4^1=Q_1/6-3c_{2,3}-6c_{2,2}-p_1,
\quad X_4^2=Q_2/6-3c_{2,2}.
\]

Periodic axions with \(Da=da+KA\) supply a **local** cancellation action, including the needed anomaly-allocation counterterm. The gravitational term is of Pontryagin type \(a^1p_1(T)\), not the Einstein kinetic term.

For positive gauge and axion kinetic matrices H and G,

\[
\mathcal M^2=H^{-1/2}K^TGKH^{-1/2},\qquad\operatorname{rank}\mathcal M^2=2.
\]

Choosing \(H=g_4^{-2}\operatorname{diag}(6,4,2,2,2)\) and \(G=f_{\rm ax}^2I\) gives the conditional masses \(g_4^2f_{\rm ax}^2(4\pm2\sqrt2)\). These are not absolute predictions.

The integer dressing criterion is \(q_{\mathcal O}=K^Tn\). It allows dressings of \(H_uH_d\), \(QQQL\), and \(u^cu^cd^ce^c\), while \(\nu^c\nu^c\) remains charged under unbroken B−L. Permission by symmetry does not determine a coefficient or proton lifetime. Minimal Stückelberg kinetic terms vanish on constant-axion, zero-vector vacua and add no classical radius potential.

**Source:** [AXG-01](studies/03_AXG_01/MTFT_AXG01_Report.md).

## 6. AXG-02 and AXG-03: why a parent is a stronger requirement

The first explicit lift put a 6D Weyl fermion in each unordered M1 bifundamental. Its anomaly pushforward exactly reproduces the 4D polynomial, yet its I8 contains irreducible \(p_2\) and mixed Abelian–color-cubic terms not canceled by ordinary products of four-forms. Changing the chirality of the degree-zero \((c,\bar d)\) sector changes the parent anomaly while leaving net 4D indices unchanged. Zero index is not permission to discard a parent field.

Native smooth fixed-modulus periodic scalars require their charge row to annihilate flux:

\[
Km=(12,12)^T,\qquad k_\Delta=k_2-k_1=(3,-2,-1,-1,1),\qquad k_\Delta m=0.
\]

Only the latter line in this shift span passes this necessary topology condition. It is not a sufficient flatness or global-holonomy condition. Tensor-derived axions have different patching rules.

For the chirality-flipped parent, AXG-03 finds a restricted gauge anomaly quadratic form with inertia \((2,1)\). This rules out a single ordinary nonchiral tensor product and the specified conventional \((1,T)\) pairing using neutral additions alone. It is not a no-go theorem for all charged spectra or all nonsupersymmetric tensor constructions.

Integral tensor flux reduction induces \(3k_1\), not \(k_1\), in the explicit construction:

\[
\widehat K=\begin{pmatrix}k_\Delta\\3k_1\end{pmatrix},\qquad
\operatorname{SNF}(\widehat K)=(1,3).
\]

The candidate Z3 fermion ledger has \(S_1=S_3=24\), failing the relevant cubic mod-9 test. An isolated charged singlet repair gives 12 additional net modes and does not fix the tensor obstruction. The full GS/topological contribution was not constructed, so the fermion-only failure is not a theorem about every completion.

The operator audit adds two crucial distinctions:

- A minimally coupled degree-−6 elementary scalar has a positive Bochner bound \(\lambda_{\min}\ge1/(4R^2)\). Eighteen Dolbeault H1 classes are not eighteen zero modes of that scalar Laplacian.
- Same-chirality 6D Weyl pairs do not admit the ordinary Lorentz-scalar Yukawa bilinear under test. A 4D charge triangle alone is insufficient to construct its 6D interaction.

With negative curvature, positive flux energy and nonnegative bulk potential, the 6D Einstein-frame radius potential is

\[
V_E=a\rho^{-4}+b\rho^{-6}+c\rho^{-2},\quad a>0,\ b,c\ge0,
\]

and has no stationary point. Negative bulk terms can yield conditional AdS minima. AXG-02's two-field exponential control requires \(\alpha/\beta>3\) under its stated definitions; the potential and exponents are inputs, not arithmetic outputs.

**Sources:** [AXG-02](studies/04_AXG_02/MTFT_AXG02_Report.md), [AXG-03](studies/05_AXG_03/MTFT_AXG03_Report.md).

## 7. AXG-04: the C3X candidate

C3X deliberately changes the parent to the **direct product**

\[
G=SU(3)\times SU(2)\times U(1)_h\times U(1)_X,\qquad h=6Y.
\]

Each fermion row occurs three times; the scalar occurs once. U,D,E,N below are physical right-handed-type representations, not their 4D left-handed conjugates.

| Field | SU(3) | SU(2) | h | X | 6D chirality |
|---|---|---|---:|---:|---:|
| Q | 3 | 2 | 1 | 1 | + |
| U | 3 | 1 | 4 | 1 | − |
| D | 3 | 1 | −2 | 1 | − |
| L | 1 | 2 | −3 | 1 | + |
| E | 1 | 1 | −6 | 1 | − |
| N | 1 | 1 | 0 | 1 | − |
| H | 1 | 2 | 3 | 0 | scalar |

There is one unit of X flux. Choose three distinct W13-fixed CM points P,Q,R with opposite u-values at P and Q, and

\[
\mathcal L_X=\mathcal O(P+Q-R).
\]

Using the inherited S0 pencil and the nonzero evaluation determinant

\[
\det\begin{pmatrix}1&i/\sqrt{13}\\1&-i/\sqrt{13}\end{pmatrix}
=-2i/\sqrt{13},
\]

Serre duality and Riemann–Roch give
\((h^0,h^1)(S_0\mathcal L_X)=(1,0)\). This is stronger than index one. The alternative \(\mathcal O(P)\) gives \((2,1)\) and retains a mirror pair. Arithmetic inputs about S0 and the CM points remain inherited; the new divisor consequence is the audited step.

### Local action and anomaly data

The report writes the Einstein, Yang–Mills, Weyl, elementary-Higgs and nonchiral two-form terms explicitly. Its scalar potential is
\(V_6=\mu_6^2H^\dagger H+\lambda_6(H^\dagger H)^2\), \(\lambda_6>0\). Opposite parent chiralities allow
\(\bar Q\widetilde H y_{6u}U\), \(\bar QH y_{6d}D\), \(\bar LH y_{6e}E\), \(\bar L\widetilde H y_{6\nu}N\).

Let \(C=c_2(SU3), W=c_2(SU2),\eta=p_1/2\), with h,x normalized Abelian curvatures. Then

\[
\boxed{I_8=(W+9h^2)(3C+W+\eta-27h^2-6x^2)=X_4Y_4.}
\]

The factors are integral on spin manifolds for the specified global gauge group. The tensor lattice is the even unimodular hyperbolic lattice

\[
\Omega=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]

In the report's source convention,

\[
a=(0,2),\quad b_C=(0,-3),\quad b_W=(-1,-1),
\quad b_{hh}=(18,-54),\quad b_{XX}=(0,-12),\quad b_{hX}=0.
\]

These are anomaly/source-lattice data, not a negative kinetic-energy prescription. A local nonchiral tensor realization uses
\(H_3=dB+\omega_3(X_4)\), \(\delta B=-\omega_2^1(X_4)\), a positive ordinary kinetic term, and \(2\pi\int B\wedge Y_4\).

The ordinary spin-bordism calculation gives
\(\Omega_7^{\mathrm{Spin}}(B(SU3\times SU2\times U1^2))=0\). Its exact mod-two matrices are supplied. This does not construct a globally defined GS partition function: nonzero integral X4, tensor patching, string sources and the relevant field-configuration category still require work.

Flux reduction gives

\[
I_6=-12x(W+9h^2),
\]

and BF/Stückelberg level −12. A candidate Z12 remains on the specified cover. The included fermionic and listed mixed discrete checks pass; full global consistency remains conditional. The old independent B−L/common-phase gauge fields have not been retained.

### Masses and the gravitational background

For the constant normalized Higgs and a common normalized internal fermion mode,

\[
\mathcal A=48\pi R^2,\quad
Y_4=\frac{y_6}{\sqrt{\mathcal A}},\quad
\lambda_{H4}=\frac{\lambda_6}{\mathcal A},\quad
v_4^2=-\frac{\mu_6^2\mathcal A}{\lambda_6}.
\]

When \(\mu_6^2<0\) is supplied, the fixed-product mass matrix is

\[
M_f=\frac{y_{6f}}{\sqrt2}\sqrt{-\frac{\mu_6^2}{\lambda_6}}.
\]

The four parent 3×3 matrices remain free inputs. Einstein-frame masses gain the appropriate radion Weyl factor when the radius varies. The old three-internal-mode M1 Yukawa tensors cannot be inserted into C3X without a new derivation.

Write \(M=M_6^4\), \(U=\Lambda_6+V_6(H_0)\), \(\rho_F>0\), and define external Ricci curvature by \(R_{\mu\nu}=\lambda_4g_{\mu\nu}\). The unwarped constant-field Einstein equations give

\[
\lambda_4=\frac{U-\rho_F}{2M},\qquad
k=\frac{U+3\rho_F}{2M}=-\frac1{R^2},
\]

so

\[
\lambda_4=-\frac1{R^2}-\frac{2\rho_F}{M}<0,
\qquad \ell_4^2=\frac3{R^{-2}+2\rho_F/M}<3R^2.
\]

The minimal solution is AdS and lacks external/internal curvature scale separation. An explicit arbitrary-unit control has \(R^2=\ell_4^2=1\), a positive radion mass squared 8 and radial Higgs mass squared 2. It is a classical control, not a prediction or proof of full KK stability.

The 72 real complex-structure/shape moduli are not selected by this truncation. A free nonchiral two-form on genus 13 supplies 26 harmonic vectors and two scalars before the specified gauging removes one scalar; the remaining modes need an explicit fate. The reduction relation \(M_P^2=M_6^4\mathcal A\) does not derive either input scale.

### Rejected alternatives remain part of the result

The one-parent/index-three six-sector alternative encounters a primitive tensor-charge integrality witness (norm −8/3; a displayed charged repair still gives −4/3). A simple constant-Higgs/shared-wavefunction limit also gives degenerate singular values. The U(8) Dirac-adjoint alternative has mirror content, a flux centralizer U(4)×U(2)×U(2), and an unstable split-bundle Yang–Mills Hessian in the stated fixed-metric test (complex Morse index 312, including 180 colored directions). These results do not by themselves establish the coupled AdS stability problem.

**Source:** [AXG-04](studies/06_AXG_04/MTFT_AXG04_Report.md), including its Higgs and larger-parent derivation notes.

## 8. What the implementation should make possible

The immediate result is a controlled research workspace: propose a spectrum and global gauge group, compute its anomaly polynomial, check the integral charge and tensor lattices, specify its mode operator, count actual zero modes, reduce its action, and test the background and mass normalization. A failed gate must return its witness and assumptions as clearly as a passed identity.

The most consequential next questions are a global GS construction for C3X; a justified mechanism changing the AdS/scale-separation result; the fate of remaining light fields and moduli; and a family/flavor mechanism that does not conceal the input multiplicity. These are listed with acceptance criteria in [R_AND_D_NEXT_GATES.md](R_AND_D_NEXT_GATES.md).
