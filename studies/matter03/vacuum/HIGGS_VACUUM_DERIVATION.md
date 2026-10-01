# Matter-03: a classical Higgs branch exists, but its scale is unselected

The single directed Matter-02 Higgs cannot be a nonzero supersymmetric vacuum of the stated compact bulk parent. This remains true after allowing an arbitrary smooth positive Hermitian bundle metric. An explicit reverse Higgs, already present in the exact mode census, repairs that obstruction. The resulting paired background is polystable and therefore admits a solution of the surface D-term equation. This is an existence theorem for a classical branch, not a computation of a Higgs expectation value, kinetic metric, particle mass, or Standard Model metastability.

## Fixed geometry and parent

Let C=X0(143), T=143a1, P=infinity the width1 cusp, Q the elliptic identity, Y=C×T, and pi:S→Y the blow-up at (P,Q), with exceptional divisor D. Put

\[
 M=pr_C^*O_C(P),\qquad L=\pi^*M\otimes O_S(-3D),
\]

with no additional flat twists. The bundle of the U(8) bulk parent is

\[
 V=O_S^{\oplus3}\oplus L^{\oplus2}\oplus L^2\oplus L^{-3}\oplus L^{-1}.
\]

Its determinant is trivial. The chosen ample wall polarization has

\[
 \mu(L)=A_E-3\epsilon=0,
 \qquad 0<\epsilon<\min(A_E,A_X).
\]

An explicit integral example is (A_E,A_X,epsilon)=(3,3,1). All summands then have slope0. The surface is smooth, compact and projective. These hypotheses matter for the correspondence used below.

The defect-free, classically twisted bulk theory uses a holomorphic canonical-valued endomorphism

\[
 \Phi\in H^0(S,End(V)\otimes K_S)
\]

and the equations F_A^(0,2)=0, bar-partial_A Phi=0, together with a curvature-plus-Higgs-commutator moment-map equation. These are the bulk equations3.32–3.34 of [Beasley–Heckman–Vafa](https://arxiv.org/pdf/0802.3391), with superpotential3.28 proportional to integral Tr(F^(0,2) wedge Phi). Normalizations and anti-Hermitian conventions change coefficients, not the following stability obstruction. No defect sources, soft terms, axion potential or quantum correction has been added.

This is a K_S-twisted surface pair. It is not an ordinary Simpson Higgs bundle with an Omega_S^1-valued field. The applicable primary theorem is [Tanaka, arXiv:1312.2673](https://arxiv.org/pdf/1312.2673), Definitions1.1–1.2 and Theorem1.3: polystability of a locally free canonical-twisted pair is equivalent to existence of a Hermitian solution of the relevant Vafa–Witten/Higgs moment-map equation. The theorem supplies an existence deduction; the certificate does not numerically solve that equation.

## Why a single directed Higgs is excluded

Write W=L^(oplus2), B=L^-1 and take only

\[
 H_d:B\longrightarrow W\otimes K_S,
 \qquad H_d\ne0.
\]

All other Phi blocks vanish. This is holomorphic and hence F-flat, but W is Phi-invariant because Phi(W)=0. Its slope0 equals that of V. For a solution of the moment-map equation, the integrated projection identity for a Phi-invariant holomorphic subbundle F is, with consistent positive norm conventions,

\[
 rank(F)\mu(V)-deg(F)
 =c_A\|\beta_F\|^2+c_\Phi\|\pi_F\Phi(1-\pi_F)\|^2,
 \qquad c_A,c_\Phi>0.
\]

Here pi_F is orthogonal projection in the proposed metric and beta_F its second fundamental form. The identity follows by projecting the curvature equation and applying the Chern–Weil subbundle formula. It is equally the equality case of the slope inequality in the correspondence. At the wall its left side vanishes for F=W, so both nonnegative terms must vanish. The orthogonal complement is then holomorphic and Phi-invariant. But Phi(V) is contained in W, and a nonzero map into W cannot also take an invariant complement into that complement. Thus Phi must vanish, a contradiction.

This excludes nonzero single-Hd solutions for every smooth positive Hermitian metric in the given holomorphic bundle and polarization. It is stronger than observing a nonzero commutator in a frozen diagonal metric. It does not assert that an arbitrary nonzero Higgs pair on this surface is excluded.

Moving the Kähler class alone does not rescue the full unchanged bundle. Let x=mu(L)=A_E−3epsilon. The inactive line subbundles L² and L^-3 are still Phi-invariant, with slopes2x and−3x. Since mu(V)=0, semistability requires

\[
 2x\le0,\qquad -3x\le0,
\]

which forces x=0. There the preceding equality obstruction applies. An isolated active L⊕L^-1 pair with a single arrow would instead require x<0, but discarding the spectator constraints is not a vacuum solution of the full rank8 construction. Independent FI shifts or additional charged condensates would change the parent assumptions and require separate analysis.

## The available reverse field

The Matter-02 exact census gives

\[
 h^0(K_S L^2)=10,\qquad h^0(K_S L^{-2})=11.
\]

Choose one weak copy L_a within W and set

\[
 V_2=L_a\oplus L^{-1},\qquad
 \Phi_2=\begin{pmatrix}0&H_d\\H_u&0\end{pmatrix}.
\]

Here H_d lies in H0(K_S L²) and H_u lies in H0(K_S L^-2), with H_u:L_a→L^-1 tensor K_S. The remaining weak copy and the other stacks have Phi0.

Both fields can be specified exactly from the frozen rational cusp-form basis. Retain the Matter-02 H_d canonical differential with coefficients a1=a2=a3=a5=0, a4=1. It has canonical order3 at P. Choose H_u with a1=a2=a4=0, a3=1, hence canonical order2. The existing cusp jet ranks prove both linear systems are solvable. The independent script supplies the rational basis vectors and checks them in exact Fraction arithmetic. It also checks the full original and adapted jet ranks1 through13 against provenance hashes.

The rational vectors in the original13-element basis are

\[
 H_d:\ (8235/34426,0,-603/17213,171/9836,81/9836,
 593/68852,0,0,0,0,0,0,0),
\]

\[
 H_u:\ (27/3083,0,-37/6166,25/6166,139/6166,
 0,0,0,0,0,0,0,0).
\]

These are algebraic coordinate normalizations. They are not unit norm in a physical L² metric.

## An odd exceptional valuation proves stability

Let alpha_d and alpha_u denote the two underlying curve canonical differentials, and eta a nonvanishing elliptic canonical differential. The relevant bundles are

\[
 K_S L^2=\pi^*(K_Y M^2)\otimes O_S(-5D),
 \qquad
 K_S L^{-2}=\pi^*(K_Y M^{-2})\otimes O_S(7D).
\]

The prescribed cusp orders are exactly sufficient to produce regular global sections of these bundles. To verify their product, use local coordinates z,w at the blow-up point and the chart z=u,w=uv. Then

\[
 dz\wedge dw=u\,du\wedge dv,
 \qquad e_L=u^3/z=u^2.
\]

The local K_S L² frame is (du wedge dv)u4 and the K_S L^-2 frame is (du wedge dv)u^-4. Relative to them, H_d has leading scalar order0 and H_u has leading scalar order7. Therefore

\[
 s=H_dH_u\in H^0(K_S^2),\qquad ord_D(s)=7.
\]

Equivalently, the underlying product alpha_d alpha_u has order3+2=5 at P, and the two pullbacks of top forms each contribute one exceptional Jacobian zero, giving5+2=7. The exact coefficient convolution independently verifies the leading coefficient and the same order. Neither the line-frame choice nor a rational trivialization of K_S changes this parity: changing the latter subtracts an even valuation from s.

Over the function field C(S), a square has even valuation at every prime divisor. The odd valuation7 therefore proves s is not a square. Consequently

\[
 \det(\lambda I-\Phi_2)=\lambda^2-s
\]

is irreducible over C(S). A rank1 Phi-invariant coherent subsheaf would supply an invariant line in the generic2-dimensional fiber and an eigenvalue in C(S), contradicting irreducibility. Thus no proper positive-rank invariant subsheaf exists, and (V2,Phi2) is stable. This argument works for the canonical-twisted pair on a surface and does not require borrowing a curve-only spectral correspondence.

At the original wall the spectators

\[
 O_S^{\oplus3},\quad L,\quad L^2,\quad L^{-3}
\]

are direct sums of stable slope0 line pairs with Phi0. The active pair is stable with determinant O and slope0. Hence the full rank8 pair is polystable. The cited surface correspondence now proves existence of a smooth Hermitian metric and unitary connection solving the classical F/D equations. No physical metric has yet been computed. Away from the wall the inactive spectators again obstruct full polystability, even though the active rank2 pair remains stable.

## What the branch does and does not select

Replace both arrows by tH_d,tH_u for any nonzero complex t. The spectral product is t²s, still of odd exceptional valuation7. The stable active pair and the polystable full bundle remain valid. At t0 the original split Phi0 wall background is also polystable. Therefore this construction supplies a continuous classical family of supersymmetric vacua. Its gauge invariant spectral coefficient changes with t²; it is not all a single holomorphic gauge orbit.

The classical equations do not choose the magnitude of t. They establish a Higgs branch rather than an isolated vacuum with a selected electroweak scale. In a4D rigid supersymmetric description, a genuine F/D-flat branch has no restoring potential along these directions. A soft or quantum potential may lift it, but must be specified and computed. This result says nothing about whether the observed Standard Model vacuum is stable or metastable at very high energies.

Giving these fields a vacuum value also changes the fluctuation complex. The original H1 and H0(K) mode spaces are not automatically the massless modes of the deformed background. A rank3 Yukawa minor is useful input into that new problem, but it does not by itself select three observed families or furnish three physical mass eigenvalues.

## Classical gauge symmetry

The initial constant U(8) gauge centralizer of the flux is U(3)×U(2)×U(1)³, of real dimension16. H_d and H_u in one common weak direction identify that direction's diagonal generator with the B generator and remove the two weak off-diagonal generators. The resulting compact classical centralizer is

\[
 U(3)\times U(1)^4,
\]

of dimension13. For an SU(8) parent impose the overall trace relation, obtaining S(U(3)×U(1)^4), of dimension12. This agrees with the Higgs automorphism description: the irreducible rank2 stable pair has only scalar automorphisms, the three trivial color copies supply GL3, and the three distinct spectator lines supply separate scalars. There are no additional global maps between distinct line powers because H0(Lq)=0 for q≠0, as certified in Matter-02.

In the ordered indices(c1,c2,c3,w1,w2,A1,A2,B), use the earlier hypercharge lift

\[
 Y=(1/6,1/6,1/6,0,0,-1/2,-1/2,1/2),
 \quad T_3=(0,0,0,1/2,-1/2,0,0,0).
\]

Then Q=T3+Y has equal value1/2 on w1 and B, so both arrows commute with Q. This proves the paired neutral branch preserves the proposed electromagnetic generator in the classical bulk Higgs mechanism. Extra Abelian generators also survive; this is not a unique-photon theorem. Stückelberg/axion terms, kinetic mixing and global completion remain separate gates.

## Reproduction and status

Run `python vacuum_certificate.py --output vacuum_exact.json`. This uses only the standard library and passes57 exact arithmetic gates. It constructs both forms, checks their archived provenance and cusp orders, verifies the exceptional valuation by two coefficient/frame routes, and computes the constant gauge-constraint rank and Q commutators. The slope no-go, generic-fiber stability proof and surface correspondence are supplied above; they are not substitutes for a numerical metric solve, nor independently formalized proofs.

The exactness boundary is: correct archived full cusp basis; specified holomorphic representatives; smooth compact projective surface; ample slope-zero wall; classical defect-free BHV canonical twist; and the cited canonical-twisted surface correspondence. Within these hypotheses, the single-arrow no-go and paired existence branch follow. The remaining physical gates include the deformed zero-mode/hypercohomology census, family and extra-mode lifting, Higgs-modulus stabilization, numerical D-term metric, positive kinetic normalization, and axionic photon kernel matching.
