# MTFT R2C-01: Six-dimensional chirality, family indices, and the Higgs interaction

**An exact finite investigation anchored to the supplied v0.32.0 archive**  
18 September 2026

## Result

With M1's original stack fluxes and one complex six-dimensional Weyl field in each of its ten oriented bifundamentals, **no chirality assignment both preserves the six named three-family sectors and permits their four ordinary elementary-scalar Yukawa bilinears**.

The exhaustive scan has 1,024 assignments. Sixteen preserve the six family indices; four preserve the complete original sector-index ledger. Sixty-four permit all four scalar Lorentz contractions. The intersection with the family-preserving assignments is empty. In fact, every family-preserving assignment forbids all four of those scalar contractions.

This is a conditional no-go result for a specified field content, flux background and interaction type. It does not rule out M1's arithmetic spaces, its section-product tensors, or every possible higher-dimensional parent. It identifies why assigning previously unspecified chiralities cannot by itself complete an elementary-scalar parent.

The investigation also independently reconstructs the **complete fermionic anomaly polynomial of the declared spectrum**, including gravitational terms. It does not include unspecified gauginos, gravitinos, chiral tensors or additional fermions.

An internal vector or form origin for the Higgs remains open at the level of Lorentz representation theory: an odd Clifford insertion has the opposite selection rule. That possibility needs an action, a physical mode operator, a full anomaly cancellation mechanism and a stable background.

## 1. Inputs and conventions

The frozen model is

\[
G=U(3)_c\times U(2)_L\times U(1)_a\times U(1)_b\times U(1)_d,
\]
\[
N=(3,2,1,1,1),\qquad m=(0,-3,3,3,0),
\]
\[
Y=\left(\tfrac16,0,-\tfrac12,\tfrac12,-\tfrac12\right),\qquad
B-L=\left(\tfrac13,0,0,0,-1\right).
\]

The global group is this direct product of unitary groups. The local SU(3)/SU(2) plus central-curvature decomposition used below is a splitting-principle calculation; it does not replace U(n) by an independently chosen SU(n) × U(1) global group.

For stack order \(c,L,a,b,d\), the oriented sector \(ij\), with \(i<j\), means

\[
R_{ij}=\mathbf N_i\otimes\overline{\mathbf N_j},\qquad
q_{ij}=e_i-e_j,\qquad d_{ij}=m_i-m_j.
\]

Assign one six-dimensional Weyl chirality \(\epsilon_{ij}\in\{+1,-1\}\) to each sector. We choose the overall chirality convention so that

\[
\boxed{n_L(R_{ij})-n_R(R_{ij})=\epsilon_{ij}\,d_{ij}.}
\]

A negative signed index can be described using left-handed fields in the conjugate gauge representation. This conjugation does not provide a freely adjustable six-dimensional chirality sign.

The ten sectors are fixed in the order

\[
(cL,ca,cb,cd,La,Lb,Ld,ab,ad,bd).
\]

| Sector | Gauge-representation dimension | Internal degree | \((h^0,h^1)\) for the specified bundles | Desired signed family index |
|---|---:|---:|---:|---:|
| cL | 6 | 3 | (3,0) | 3, giving \(Q\) |
| ca | 3 | −3 | (0,3) | −3, giving \(u^c\) after conjugation |
| cb | 3 | −3 | (0,3) | −3, giving \(d^c\) after conjugation |
| cd | 3 | 0 | (2,2) | — |
| La | 2 | −6 | (0,6) | — |
| Lb | 2 | −6 | (0,6) | — |
| Ld | 2 | −3 | (0,3) | −3, giving \(L\) after conjugation |
| ab | 1 | 0 | (2,2) | — |
| ad | 1 | 3 | (3,0) | 3, giving \(\nu^c\) |
| bd | 1 | 3 | (3,0) | 3, giving \(e^c\) |

The six-family condition is deliberately weaker than requiring a complete realistic low-energy spectrum. The La/Lb sectors and the index-zero sectors remain in the calculation. Vector-like under the Standard Model subgroup does not automatically mean a mass term is invariant under every stack U(1).

### Why these are pure three-mode sectors

Let \(\Delta=P_1+P_2+P_3\), with distinct chosen W13-fixed CM points containing both signs of \(u\). The input from the existing arithmetic construction is

\[
S_0^2=K,\quad \deg S_0=12,\quad
H^0(S_0)=\operatorname{span}\{1,u\},\quad
u(P_i)=\pm i/\sqrt{13}.
\]

The stack line bundles are \((\mathcal O,\mathcal O(-\Delta),\mathcal O(\Delta),\mathcal O(\Delta),\mathcal O)\). A representative evaluation matrix is

\[
E=\begin{pmatrix}
1&i/\sqrt{13}\\
1&i/\sqrt{13}\\
1&-i/\sqrt{13}
\end{pmatrix}.
\]

The minor using rows 1 and 3 is \(-2i/\sqrt{13}\ne0\). Thus \(H^0(S_0(-\Delta))=0\). Serre duality and Riemann–Roch give

\[
h^1(S_0(\Delta))=0,\qquad
h^0(S_0(\Delta))-h^1(S_0(\Delta))=\deg\Delta=3.
\]

The dual bundle has \((h^0,h^1)=(0,3)\). Replacing \(\Delta\) by \(2\Delta\) gives (6,0) and (0,6). These are exact deductions from the frozen arithmetic input; this study does not independently reconstruct the theta characteristic or CM points.

For the trivial degree-zero difference, \(h^0(S_0)=h^1(S_0)=2\). Therefore **index zero does not imply an empty sector**. A generic flat twist can change these kernels, but it is an additional bundle choice whose effects on the other sectors must also be checked.

## 2. The scalar obstruction, reproducible with matrices

Let \(\Gamma_*\) be the six-dimensional chirality operator and

\[
P_\epsilon=\frac{1+\epsilon\Gamma_*}{2}.
\]

The Dirac conjugate of a Weyl field carries the complementary projector. Hence the local scalar contraction contains

\[
\overline{\Psi_\epsilon}\,\Chi_{\epsilon'}
=\overline\Psi P_{-\epsilon}P_{\epsilon'}\Chi,
\]

which vanishes for \(\epsilon'=\epsilon\). A nonzero ordinary scalar bilinear requires opposite six-dimensional chiralities.

The transpose/charge-conjugated contraction has the same scalar selection rule. Here is a complete finite matrix recipe to check it. Define

\[
\sigma_1=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
\sigma_2=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
\sigma_3=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\]

With \(I=I_2\), take

\[
\begin{aligned}
\gamma_1&=\sigma_1\otimes I\otimes I,&
\gamma_2&=\sigma_2\otimes I\otimes I,\\
\gamma_3&=\sigma_3\otimes\sigma_1\otimes I,&
\gamma_4&=\sigma_3\otimes\sigma_2\otimes I,\\
\gamma_5&=\sigma_3\otimes\sigma_3\otimes\sigma_1,&
\gamma_6&=\sigma_3\otimes\sigma_3\otimes\sigma_2.
\end{aligned}
\]

Then

\[
\{\gamma_a,\gamma_b\}=2\delta_{ab}I_8,\qquad
\Gamma_*=i\gamma_1\gamma_2\gamma_3\gamma_4\gamma_5\gamma_6
=\sigma_3\otimes\sigma_3\otimes\sigma_3,
\]
\[
\mathcal C=\gamma_2\gamma_4\gamma_6,\qquad
\mathcal C\gamma_a\mathcal C^{-1}=-\gamma_a^T,\qquad
P_\epsilon^T\mathcal C P_\epsilon=0.
\]

These are matrices for the complexified Clifford algebra. They test algebraic Lorentz invariant tensors; no Euclidean dynamical parent is being substituted for the Lorentzian theory. Equivalently, for \(\operatorname{Spin}(6,\mathbb C)\cong SL_4(\mathbb C)\), \(4\otimes4=6\oplus10\) has no scalar, whereas \(4\otimes4^*=1\oplus15\) does.

The exact matrix ranks are:

| \(\epsilon,\epsilon'\) | Scalar \(P_{-\epsilon}P_{\epsilon'}\) | Transpose scalar \(P_\epsilon^T\mathcal C P_{\epsilon'}\) | Vector \(P_{-\epsilon}\gamma_5P_{\epsilon'}\) |
|---|---:|---:|---:|
| −,− | 0 | 0 | 4 |
| −,+ | 4 | 4 | 0 |
| +,− | 4 | 4 | 0 |
| +,+ | 0 | 0 | 4 |

All six expanded 8 × 8 gamma matrices, \(\Gamma_*\), and \(\mathcal C\) are included in results/clifford_matrices.json.

### Applying the rule to the four Higgs triangles

The four-dimensional left-handed charges are

\[
\begin{aligned}
Q&=e_c-e_L,&u^c&=e_a-e_c,&d^c&=e_b-e_c,\\
L&=e_d-e_L,&\nu^c&=e_a-e_d,&e^c&=e_b-e_d,\\
H_u&=e_L-e_a,&H_d&=e_L-e_b.
\end{aligned}
\]

Consequently \(Qu^cH_u\), \(Qd^cH_d\), \(L\nu^cH_u\) and \(Le^cH_d\) have zero total stack charge. Gauge invariance is satisfied. Their elementary-scalar six-dimensional Lorentz contractions separately require

\[
\epsilon_{cL}=-\epsilon_{ca}=-\epsilon_{cb},\qquad
\epsilon_{Ld}=-\epsilon_{ad}=-\epsilon_{bd}.
\]

Preserving the six original signed indices instead forces

\[
\epsilon_{cL}=\epsilon_{ca}=\epsilon_{cb}
=\epsilon_{Ld}=\epsilon_{ad}=\epsilon_{bd}=+1.
\]

These conditions cannot hold together. The up-quark pair alone proves the contradiction. Regrouping the extra lepton-like doublets does not change that quark-sector witness.

## 3. Exhaustive assignment results

| Requirement | Number of assignments |
|---|---:|
| Arbitrary signs for ten Weyl sectors | 1,024 |
| Preserve the six named family indices | 16 |
| Preserve all original M1 sector indices | 4 |
| Permit all four elementary-scalar Lorentz contractions | 64 |
| Preserve six family indices and permit all four | **0** |
| Preserve six family indices and retain vanishing 4D anomaly on span{phase, Y, B−L} | 8 |

The counts have short independent combinatorial proofs:

- Six fixed family signs leave four free signs: \(2^4=16\).
- Preserving the two additional nonzero-index La/Lb sectors fixes two more signs: \(2^2=4\).
- Scalar compatibility leaves one sign for the quark triple, one for the lepton triple, and four isolated sector signs: \(2^6=64\).

The eight anomaly-free members of the 16-family set have \(\epsilon_{La}=\epsilon_{Lb}\). The other eight alter the additional doublet balance and fail this particular 4D anomaly check. This does not affect the scalar no-go: it already holds on the larger set of 16.

Every row, including all ten signs and signed indices, is recorded in results/assignments.csv. Allowing a contraction is only a necessary Lorentz condition; it does not establish a nonzero overlap integral, Higgs vacuum expectation value or physical mass.

## 4. Complete fermionic anomaly reconstruction

Write

\[
C=c_2(SU(3)),\quad D_3=c_3(SU(3)),\quad W=c_2(SU(2)),
\]

and let \(f_i\) denote the normalized central curvature in the fundamental of stack \(i\). Define

\[
A_8=\frac{7p_1^2-4p_2}{5760},\qquad
H(q)=\frac{q^4}{24}-\frac{p_1q^2}{48}+A_8,\qquad
T(q)=\frac{q^2}{2}-\frac{p_1}{24}.
\]

For a positive complex Weyl fermion,

\[
I_8(R)=[\widehat A(T)\operatorname{ch}(R)]_8.
\]

Four formulas reproduce all ten sector polynomials by hand:

\[
\begin{aligned}
F_1(q)&=H(q),\\
F_3(q)&=3H(q)-C\,T(q)+\frac{D_3q}{2}+\frac{C^2}{12},\\
F_2(q)&=2H(q)-W\,T(q)+\frac{W^2}{12},\\
F_{32}(q)&=6H(q)-(2C+3W)T(q)+D_3q
 +\frac{C^2}{6}+\frac{W^2}{4}+CW.
\end{aligned}
\]

Use \(q=f_i-f_j\) and the following map:

| Sectors | Polynomial |
|---|---|
| cL | \(F_{32}(q)\) |
| ca, cb, cd | \(F_3(q)\) |
| La, Lb, Ld | \(F_2(q)\) |
| ab, ad, bd | \(F_1(q)\) |

Then \(I_8(\epsilon)=\sum_{i<j}\epsilon_{ij}I_8(R_{ij})\). The bundle contains every expanded sector expression and the four complete-ledger survivor polynomials.

### An independent route

For color roots \((\alpha,\beta,-\alpha-\beta)\), weak roots \((\gamma,-\gamma)\), and zero roots for the three rank-one stacks, enumerate every bifundamental weight

\[
w=f_i-f_j+r_i-r_j.
\]

Sum \(w^4/24-p_1w^2/48+A_8\), and substitute

\[
C=-(\alpha^2+\alpha\beta+\beta^2),\quad
D_3=-\alpha\beta(\alpha+\beta),\quad W=-\gamma^2.
\]

This agrees with the Chern-character construction in each of the ten sectors. Linearity extends that identity to every chirality assignment.

The four-dimensional polynomial is independently built from

\[
I_6(R)=\sum_{\text{weights }w}\left(\frac{w^3}{6}-\frac{p_1w}{24}\right).
\]

Every sector satisfies the symbolic pushforward identity

\[
\boxed{\sum_km_k\frac{\partial I_8(\epsilon)}{\partial f_k}
=\sum_{i<j}\epsilon_{ij}(m_i-m_j)I_6(R_{ij}).}
\]

This is why **three internal modes multiply the reduced four-dimensional anomaly, not the uncompactified six-dimensional fermion count**. The two quantities must be tracked separately.

### Irreducible gravitational term

For the 16 family-preserving assignments, put

\[
r=\epsilon_{cd},\quad s=\epsilon_{ab},\quad
\ell_a=\epsilon_{La},\quad\ell_b=\epsilon_{Lb}.
\]

Their signed six-dimensional representation dimension is

\[
n_{\rm grav}=16+3r+2\ell_a+2\ell_b+s\ge8,
\qquad [p_2]I_8=-\frac{n_{\rm grav}}{1440}\ne0.
\]

For the four assignments preserving the whole original ledger, \(\ell_a=\ell_b=1\):

| \(\epsilon_{cd}\) | \(\epsilon_{ab}\) | \(n_{\rm grav}\) | \([p_2]I_8\) |
|---:|---:|---:|---:|
| −1 | −1 | 16 | −1/90 |
| −1 | +1 | 18 | −1/80 |
| +1 | −1 | 22 | −11/720 |
| +1 | +1 | 24 | −1/60 |

Ordinary Green–Schwarz counterterms built from products of characteristic 4-forms cannot remove an irreducible \(p_2\) term. Additional chiral fields, including possible tensor or gravity-sector contributions, can change the total polynomial and must be specified explicitly. This is an obstruction to the declared spectrum as a standalone anomaly-free gravitational parent, not a calculation of a graviton or Newton's constant.

### Color-cubic term

For any assignment,

\[
\frac{\partial I_8}{\partial D_3}
=\frac12\left[
2\epsilon_{cL}(f_c-f_L)+
\epsilon_{ca}(f_c-f_a)+
\epsilon_{cb}(f_c-f_b)+
\epsilon_{cd}(f_c-f_d)\right].
\]

It cannot vanish as a polynomial in independent stack curvatures: the coefficient of \(f_a\), for example, is \(-\epsilon_{ca}/2\).

For the family-preserving assignments this becomes

\[
\frac12\left[(4+r)f_c-2f_L-f_a-f_b-rf_d\right].
\]

A characteristic 4-form factor cannot contain \(D_3\), which has degree six. Thus ordinary 2-form Green–Schwarz factorization alone cannot cancel this term. A separately declared shifting scalar with a 6-form coupling, additional charged matter, or a changed gauge-group construction is a different mechanism. Neutral spectator fermions do not change this color term.

## 5. Controlled change: allow the stack fluxes to move

To see what the scalar rule demands, keep the same six named family sectors and charges but allow a new flux vector. Set

\[
\epsilon_{cL}=a,\quad \epsilon_{ca}=\epsilon_{cb}=-a,\qquad
\epsilon_{Ld}=b,\quad \epsilon_{ad}=\epsilon_{bd}=-b,
\quad a,b\in\{\pm1\}.
\]

Fix the irrelevant common flux shift by \(m_c=0\). Solving the six target-index equations gives exactly

\[
\boxed{m(a,b)=(0,-3a,-3a,-3a,\,3(b-a)).}
\]

| \(a\) | \(b\) | New flux vector | \(d_{cd}\) | Six-sector required tensor norm below |
|---:|---:|---|---:|---:|
| +1 | +1 | (0,−3,−3,−3,0) | 0 | −8/3 |
| +1 | −1 | (0,−3,−3,−3,−6) | 6 | −4/3 |
| −1 | +1 | (0,3,3,3,6) | −6 | 4/3 |
| −1 | −1 | (0,3,3,3,0) | 0 | 8/3 |

None is the original M1 flux. If the cd sector is retained, opposite choices \(a=-b\) introduce six net colored modes. For all four branches, \(d_{La}=d_{Lb}=d_{ab}=0\).

The family sectors can still use pure degree-\(\pm3\) bundles from the same divisor construction. However, both Higgs line degrees are now zero. With the explicit untwisted bundles, the previous space

\[
H^1(\mathcal O(-2\Delta)),\quad \dim=18
\]

has been replaced by \(H^1(\mathcal O)\), of dimension 13. An elementary scalar on the trivial bundle would instead have the constant covariant-Laplacian zero mode. Neither statement reproduces the original 18-dimensional Higgs tensor construction. Flat twists can change this further and must be recorded.

### Primitive tensor-charge witness

First retain only the six named family sectors; this removal of the other four sectors is an explicit control. Their pure gravitational anomaly cancels, and their color-cubic direction is proportional to

\[
k_1=(0,-2,1,1,0).
\]

The primitive cocharacter

\[
t\longmapsto(f_c,f_L,f_a,f_b,f_d)=(0,0,t,-t,0)
\]

lies in \(\ker k_1\). Setting the non-Abelian and gravitational characteristic classes to zero gives

\[
I_8\big|_t=-\frac{3a+b}{12}\,t^4.
\]

Under the explicitly assumed ordinary integral tensor-charge lattice convention

\[
I_8=\tfrac12\Omega(X_4,X_4),\qquad
X_4=\tfrac12 b_{tt}t^2+\cdots,
\]

one necessarily needs

\[
\Omega(b_{tt},b_{tt})=8[t^4]I_8
=-\frac{2(3a+b)}3\in\mathbb Z.
\]

Every branch fails this necessary condition. This generalizes the previously studied \((a,b)=(1,1)\) AXG-04 witness to all four scalar-compatible core branches. The role of the global gauge group and integral string-charge data is explained in [Monnier, Moore and Park, §§2–3](https://arxiv.org/pdf/1711.04777). Their supergravity framework is a reference for the stated quantization assumptions, not a claim that M1 has acquired supersymmetry.

The net number of six-dimensional complex SU(2) doublets in a core is \(3a+b=\pm4,\pm2\). This also fails the usual mod-six condition for that restricted fermion content. These are bulk field multiplicities. The six-dimensional assignments and the familiar generation-count argument appear in [Dobrescu and Poppitz, equations (1)–(4)](https://arxiv.org/pdf/hep-ph/0102010); they must not be reinterpreted as multiplying one bulk field by its flux zero-mode degeneracy. The modern interpretation of those congruences requires a properly quantized Green–Schwarz construction; they are not an additional torsion bordism claim after full anomaly cancellation. See [Lee and Tachikawa, §§1 and 2.7](https://arxiv.org/pdf/2012.11622).

### Retain all ten sectors instead

We also checked all 16 sign choices of cd, La, Lb and ab for each changed-flux core: 64 completions in total.

With no other chiral fields, cancelling the fermionic pure gravitational term forces

\[
\epsilon_{cd}=\epsilon_{ab}=r,\qquad
\epsilon_{La}=\epsilon_{Lb}=-r,\qquad r=\pm1.
\]

Thus there are two such completions per branch. Their required norm along the same primitive cocharacter is

\[
\Omega(b_{tt},b_{tt})=-2a-\frac{2b}{3}+4r,
\]

which is never an integer. None of these eight gravitationally balanced completions passes this norm test.

This conclusion assumes that the tested U(1) background remains allowed and its anomaly must be cancelled by the stated tensor mechanism. Extra shifting scalars acting on that direction, a different allowed background category, new matter, inflow or a different global group require a new calculation. The result is not a no-go for arbitrary extensions.

## 6. What is left open at the interaction level

The Clifford table shows

\[
P_{-\epsilon}\Gamma^MP_\epsilon\ne0
\]

for an odd gamma insertion. If the four-dimensional Higgs originates from the internal components of a six-dimensional Lorentz vector, the scalar-bilinear obstruction does not apply to the corresponding vector interaction.

For example, let \(\mathcal H_{uM}\) and \(\mathcal H_{dM}\) have stack charges \(e_L-e_a\) and \(e_L-e_b\). The following are formal gauge- and Lorentz-invariant interaction structures for equal chirality within each family triple:

\[
\begin{aligned}
\mathcal L_{\rm int}={}&
\lambda_u\overline{\Psi}_{cL}\Gamma^M\Psi_{ca}\mathcal H^\dagger_{uM}
+\lambda_d\overline{\Psi}_{cL}\Gamma^M\Psi_{cb}\mathcal H^\dagger_{dM}\\
&+\lambda_\nu\overline{\Psi}_{Ld}\Gamma^M\Psi_{ad}\mathcal H_{uM}
+\lambda_e\overline{\Psi}_{Ld}\Gamma^M\Psi_{bd}\mathcal H_{dM}
+\mathrm{h.c.}
\end{aligned}
\]

Gauge-index contractions are implicit; every stack charge sum is checked explicitly. Selecting internal components occurs after compactification, not by inserting a preferred internal vector into an otherwise unspecified Lorentz-invariant scalar action.

These expressions establish allowed tensor structures only. A charged vector field needs a consistent kinetic sector and parent origin. Its physical internal Hessian need not equal the scalar covariant Laplacian or the Dolbeault Laplacian. Magnetized gauge examples explicitly distinguish these operators; see [Cremades, Ibáñez and Marchesano, equations (3.65)–(3.70)](https://arxiv.org/pdf/hep-th/0404229). Their toroidal calculation is not a spectrum calculation on \(X_0(143)\).

The next useful experiment is therefore to specify a parent representation and action for these form-valued Higgs modes while preserving the original family bundles, then derive the quadratic mode operator and overlap map. It must also supply the missing anomaly-cancelling content. The previously rejected U(8) adjoint construction is not revived by this Clifford observation.

## 7. C3X check and the meaning of three families

As a positive algebraic control, take the six-sector \((a,b)=(1,1)\) charge sum, restrict

\[
f=6Yh+(0,-1,-1,-1,0)x,
\]

and take three copies. The independently reconstructed polynomial is

\[
I_8^{\rm C3X}
=(W+9h^2)\left(3C+W+\frac{p_1}{2}-27h^2-6x^2\right).
\]

This calculation starts from fermion weights and reaches the factorization. It supplements the packaged identity check that starts from already defined factors.

It verifies this local polynomial only. It does not establish a global Green–Schwarz functional, a vacuum or an absolute scale. C3X continues to have **three input copies × one internal mode**. M1's intended target continues to be **one parent × three internal modes**.

## 8. Reproduction and integration record

The supplied archive's SHA-256 is

~~~text
46ffc563a0cd81a7a82d975497bf64d82f8471eb424409e3d154532b5d8d6069
~~~

It differs from the release hash quoted in the conversation. This report is pinned to the supplied archive; PROVENANCE.json records both hashes and the hashes of every included source snapshot.

Run from the extracted study directory:

~~~bash
python -m pip install -r requirements.txt
python chirality_audit.py --output results
~~~

The tested environment is Python 3.12.14 and SymPy 1.14.0. The audit does not import the MTFT package, use numerical tolerances, or require NumPy/SciPy. **293 exact checks pass.** This is this study's check count, not a rerun of the release-wide test suite.

| File | Purpose |
|---|---|
| chirality_audit.py | Full executable finite scan and anomaly/Clifford/control calculations |
| results/assignments.csv | All 1,024 sign assignments and their outcomes |
| results/audit_results.json | Counts, witnesses, source anchor and changed-flux controls |
| results/sector_anomaly_basis.json | Complete I8 and I6 building blocks for every sector |
| results/clifford_matrices.json | Expanded matrices and contraction ranks |
| results/verification.json | Named exact checks and independent calculation routes |
| HANDOFF.md | Proposed module contracts and regression witnesses for integration |
| ANALYSIS_SPECIFICATION.md | Fixed scope recorded before the enumeration |
| PROVENANCE.json | Source hashes and archive identity |

## Glossary

| Term | Meaning in this study |
|---|---|
| Six-dimensional chirality | The sign of the Weyl representation of the parent Lorentz group. |
| Internal chirality | The two-dimensional spinor grading represented here by H0 versus H1. |
| Signed index | Left minus right zero modes in one specified gauge representation; it is not the total number of modes. |
| Pure family sector | A kernel with (3,0) or (0,3), so that the index counts all modes in that sector. |
| Bifundamental | A field fundamental under one stack and antifundamental under another. |
| Elementary scalar | A six-dimensional Lorentz scalar, as distinct from an internal component of a vector. |
| Anomaly polynomial | A characteristic-class polynomial encoding the perturbative anomaly of specified chiral fields. |
| Irreducible term | A term outside the products of 4-forms available to ordinary 2-form Green–Schwarz cancellation. |
| Cocharacter | A homomorphism U(1) → G, fixing an allowed integral gauge background. |
| Tensor-charge lattice | The integral charge lattice and bilinear pairing used in a quantized tensor cancellation mechanism. |
| Necessary test | A failure excludes the stated class; passing it alone does not establish a consistent parent. |

The durable outcome is a sharper boundary between the arithmetic three-mode construction and a physical parent action: the flux already fixes the family chiralities, and an elementary-scalar interaction cannot be added just by choosing signs.
