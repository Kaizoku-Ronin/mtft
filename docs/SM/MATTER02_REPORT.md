# Matter-02: exact modes and a rank-three holomorphic coupling

1 October 2026 · MTFT 0.34.0 + unreleased Matter-01/02 additions

The specified collinear blow-up background admits a nonzero bulk lepton–Higgs coupling with **holomorphic rank at least three**. This is a stronger result than a charge-selection rule. It uses actual cohomology classes and a nonzero exact residue minor. The physical electron mass, the selection of three light families, and the massless photon remain open.

Matter-01 supplied a measured-input QED reference and showed that net family indices alone do not determine the mode spectrum. Matter-02 specifies a holomorphic representative, computes its mode dimensions, and tests the next interaction gate. It preserves the released source and the earlier archived calculations.

## Background and assumptions

Let C = X₀(143), a closed genus-13 modular curve, and E = 143a1, with elliptic origin O. Choose the rational cusp P = ∞ on C. Let π : S → C × E be the blow-up at (P,O), with exceptional divisor D. This is an explicit new choice of representative; it is not the earlier speculative CM point. Set

\[
M=\operatorname{pr}_C^*\mathcal O_C(P),\qquad
\mathcal L=\pi^*M\otimes\mathcal O_S(-3D).
\]

There is no nontrivial flat twist on E. In the archived lattice convention, c₁(ℒ) = (1,0,−3). The five stack exponents, ordered (c,L,A₁,A₂,B), are (0,1,2,−3,−1), with ranks (3,2,1,1,1). Their rank-weighted exponent sum vanishes. A Hom block from stack j to stack i has exponent q = nᵢ−nⱼ.

The common slope is zero when ε = A_E/3. The explicit ample example (A_X,A_E,ε) = (3,3,1) satisfies this relation. This specifies a Kähler class, not an explicit Hermitian metric or numerical harmonic wavefunctions. The split background has Φ_background = 0. Its canonical class is K_S = π* K_(C×E) + D; K_S² = −1 and χ(𝒪_S) = 0.

The interaction calculation assumes the partially twisted eight-dimensional Yang–Mills bulk field dictionary and superpotential of Beasley–Heckman–Vafa (BHV), not a derived global F-theory embedding. A spinᶜ/R-twist and ultraviolet completion still require an explicit parent. This internal surface has real dimension four; the earlier C3X breathing calculation had a different internal dimension. Its numerical mass and coupling coefficients are not transferred here.

## Exact cusp certificate and modes

The archived 13-dimensional weight-two cusp-form basis gives differential f(z) dz/z at the width-one cusp z=q. Vanishing to order m means that a₁,…,a_m vanish. Both the original and Atkin–Lehner adapted bases have exact jet rank m for m = 1,…,13. The full first-13-row determinant is

\[
-18\,196\,056\,815\,164\,784\,640.
\]

Consequently P is ordinary: the canonical vanishing orders are 0,…,12, and h⁰(K_C(−mP)) = max(13−m,0). The certificate checks the frozen input hashes and the basis transformation across all 141 archived coefficient rows. It relies on the archive being the stated complete modular-form basis; it does not regenerate that basis from PARI/GP in this milestone.

For nonzero integral q,

\[
h^0(\mathcal L^q)=0,\qquad
h^2(\mathcal L^q)=
\begin{cases}
\max(13-q,0),&q>0,\\
\max(14+2q,0),&q<0,
\end{cases}
\qquad
h^1=h^2+\frac{9q^2+3q}{2}.
\]

For q = 0 the dimensions are (1,14,13). These formulas apply to this specified cusp representative and trivial elliptic flat factor. They are not dimensions for arbitrary representatives of the same Néron–Severi class.

| q | h⁰(ℒᑫ) | h¹(ℒᑫ) | h²(ℒᑫ) | h⁰(K_Sℒᑫ) |
|---:|---:|---:|---:|---:|
| −5 | 0 | 109 | 4 | 8 |
| −4 | 0 | 72 | 6 | 9 |
| −3 | 0 | 44 | 8 | 10 |
| −2 | 0 | 25 | 10 | 11 |
| −1 | 0 | 15 | 12 | 12 |
| 0 | 1 | 14 | 13 | 13 |
| 1 | 0 | 18 | 12 | 12 |
| 2 | 0 | 32 | 11 | 10 |
| 3 | 0 | 55 | 10 | 8 |
| 4 | 0 | 87 | 9 | 6 |
| 5 | 0 | 128 | 8 | 4 |

The last column follows from Serre duality: h⁰(K_Sℒᑫ) = h²(ℒ⁻ᑫ). In particular, charged Higgs-type Φ modes survive. Inferring their absence from h⁰(ℒᑫ) = 0 would be incorrect.

Summing h¹ over the twenty off-diagonal stack Hom orientations gives 864; the corresponding h² sum is 198. These are unweighted vector-space dimensions, not counts of particles or independent observed families. They refine the earlier Riemann–Roch lower bounds rather than contradicting them.

For q>0, the pushforward π_*ℒᑫ is Mᑫ I_p^(3q). Base sections depend only on the curve coordinate and cannot have a zero of order 3q when their total divisor degree is q. Serre duality gives h² from K_C(−qP). For q=−r<0, the dual canonical bundle pushes to (K_C(rP) ⊠ 𝒪_E) I_p^(3r−1); its sections are independent of the elliptic coordinate, giving K_C(−(2r−1)P). Surface Riemann–Roch then gives h¹. The full Leray and local-jet derivation is in COHOMOLOGY_DERIVATION.md.

## A convention correction before assigning physical chirality

The archived function `net` uses χ(ℒᑫ)−χ(ℒ⁻ᑫ) = −3q. In the BHV chiral convention, for nonzero q the chiral space has dimension

\[
N(q)=h^1(\mathcal L^q)+h^0(K_S\mathcal L^q)
     =h^1(q)+h^2(-q),\qquad N(q)-N(-q)=+3q.
\]

These are opposite named orientations. The released index remains unchanged. The new API exposes both explicitly; assigning a physical representation versus its conjugate must fix the convention. A net index of magnitude three does not identify exactly three light Standard Model families. For q=0, the H⁰ gauge/vector mode is reported separately from adjoint chiral modes.

## The global coupling and its exact minor

For the selected bulk parent,

\[
W\propto\int_S\operatorname{Tr}(F^{0,2}\wedge\Phi^{2,0}),\qquad
H^1(\mathcal L)\otimes H^1(\mathcal L^{-3})
\otimes H^0(K_S\mathcal L^2)\longrightarrow\mathbb C.
\]

| Electron-type block | q | Hypercharge Y | Cohomology dimension |
|---|---:|---:|---:|
| A₁ ← L, lepton-type a | 1 | −1/2 | h¹ = 18 |
| B ← A₁, antilepton-type a | −3 | +1 | h¹ = 44 |
| L ← B, down-Higgs-type Φ | 2 | −1/2 | h⁰(K_Sℒ²) = 10 |

Both charges close. The U(8) root bracket is nonzero: Tr([E_(A₁,L), E_(B,A₁)] E_(L,B)) = −1. Gauge selection therefore does not kill this term.

In local coordinates z on C and w on E, the positive q=1 jet quotient has monomials of total degree below three. Its single global-section image is z. Choose the independent classes w, w², zw. On the negative q=−3 side, the principal parts z⁻⁶w⁻², z⁻⁷w⁻², z⁻⁶w⁻³ are allowed exceptional classes and have zero global obstruction because the global dual sections have no w dependence.

A canonical differential α vanishing exactly to order three at P exists by the exact jet ranks. Relative to the K_Cℒ² frame, its coefficient is

\[
c(z)=v_5z^5+v_6z^6+\cdots,  v_5\ne0.
\]

Taking the coefficient of z⁻¹w⁻¹ in each product produces

\[
Y_{\rm minor}=
\begin{pmatrix}
v_5&v_6&0\\
0&0&v_5\\
0&v_5&0
\end{pmatrix},\qquad \det Y_{\rm minor}=-v_5^3\ne0.
\]

One exact rational combination of the archived basis sets a₁=a₂=a₃=a₅=0 and a₄=1. It gives v₅=1, v₆=0 and determinant −1. An independent calculation chooses a different differential with v₆=1186/3083 and obtains the same determinant. These are basis-dependent holomorphic units, not a physical Yukawa value of one.

The argument is global: multiplication by Φ sends the three positive fat jets to independent classes in the target fat-jet cokernel, and global Serre duality supplies the negative partners. It is not merely a visually plausible local residue. This proves a rank lower bound for a slice of the full 18 × 44 × 10 holomorphic interaction tensor. It does not compute the full tensor or select the observed three-family subspace.

No Higgs expectation value is chosen. A nonzero off-diagonal vev can change the BPS equations and the zero-mode complex, so the split-spectrum dimensions cannot simply be carried unchanged into a recombined vacuum. Positive kinetic normalization preserves the rank of the full coupling under invertible changes of field basis, but its physical singular values and masses require the Hermitian Gram matrices and a vacuum solution.

## Photon gate

The traceless hypercharge lift is Y = (1/6,0,−1/2,−1/2,1/2). Exact fundamental-trace projections give

\[
\operatorname{Tr}(Yc_1)=0,\qquad
\int_S\operatorname{Tr}(Yc_1^2)=54.
\]

A common-center lift Y+sI keeps adjoint root charges fixed. Its quadratic projection is 54−144s, which vanishes at s=3/8, while Tr(Y+sI)=8s. An independent trace axion would instead force s=0. Thus the vanishing linear flux projection alone cannot establish photon masslessness. Neither does the nonzero quadratic projection alone establish a photon mass: the parent must say which gauged axion terms exist.

The photon audit supplies explicitly labelled INPUT matrices to demonstrate the target condition: after electroweak breaking the joint axion/Higgs mass kernel should be the electromagnetic direction Q. It also checks that the rank-two mixed-anomaly kernel contains a cubic-anomalous direction; mixed-anomaly cancellation alone cannot define a Green–Schwarz charge matrix. An ordinary real-adjoint local I₁₀=0 calculation does not settle spinᶜ/R-twist or global anomalies. Full kinetic, shift, counterterm and duality data remain missing. PHOTON_KERNEL_REPORT.md separates the exact traces from these hypothetical action probes.

## Implementation and reproduction

The new `mtft.matter.collinear` module provides exact dimensions, cusp-jet certificates and the rank-three witness. The standalone cohomology and Yukawa certificates use Python's standard library, including rational arithmetic; no new PARI calculation is required. The photon audit and integration tests use the repository's existing Python dependencies. The lab is a self-contained local HTML file.

From the candidate repository:

```bash
python -m pip install -e .
python -m pip install pytest matplotlib
python studies/matter02/run_all.py
```

The archived focused run passes **220 exact arithmetic gates** (148 cohomology, 34 Yukawa, 38 photon), **133 Python tests**, and **138 numerical/mock-DOM lab checks**. These verify finite inputs and algebra accompanying the geometric proofs; they are not a formal sheaf-cohomology proof assistant certificate. The standalone scientific overview was rendered and visually checked.

The runner prints verdicts and retains full logs/JSON under `studies/matter02/results/`. The focused gates cover the new exact calculations and the earlier Matter-01 QED/reference tests. They do not run the entire historical release suite or its hours-long investigations. The archived source basis, derived certificates, scripts, tests, reports and integration patches are included in the handoff bundle. Node checks exercise lab logic with a mock DOM; browser layout remains unverified.

## Next gate

The next physical calculation is to solve a Higgs/recombined vacuum and determine which modes remain light, with a consistent chirality dictionary. Then construct the Hermitian mode metrics and evaluate the normalized coupling. In parallel, derive the neutral gauge/axion mass operator from one complete parent action. These are separate requirements: neither a nonzero holomorphic cubic nor a visual analogy with geometry replaces them.

The older Feigenbaum/oloid questions can later supply geometric or dynamical diagnostics once an evolution operator and invariant are specified. They are not used as inputs to this exact coupling certificate. Matter-02 does not establish a Lorentzian spacetime, gravitational backreaction, or a unification theorem.

## Primary references

- Beasley, Heckman and Vafa, *GUTs and Exceptional Branes in F-theory—I*, arXiv:0802.3391v2, §§3.3.2–3.3.3, equations (3.44)–(3.55), and Appendix C: https://arxiv.org/pdf/0802.3391 . Parent field dictionary and bulk superpotential.
- Patrick Brosnan, *UMD 808J: Algebraic Surfaces*, §4.2, propositions 4.7–4.11: https://terpconnect.umd.edu/~pbrosnan/notes/Surfaces/sect0028.html . Blow-up intersection, canonical-class and direct-image facts.
- Frozen MTFT modular-form archive: `src/mtft/canonical/_data/PROVENANCE.txt` and the X₀(143) basis files; content hashes are in the certificates. New mode/rank results are derived here, not attributed to BHV.

