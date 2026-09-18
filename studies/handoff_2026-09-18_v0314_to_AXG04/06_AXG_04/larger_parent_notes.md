# AXG-04 control: a U(8) parent with an adjoint Dirac fermion

This is a deliberately specified, nonsupersymmetric effective field theory. It supplies a single gauge action and charged internal gauge components, but it does **not** supply a stable chiral Standard Model vacuum. All numerical statements below concern the original M1 flux, not a replacement flux used in other AXG-04 branches.

## 1. An explicit action and representation

Take a smooth compact spin curve X of genus 13 and a six-dimensional Lorentzian manifold locally of product form M4 × X. A control action is

\[
S=\int d^6x\sqrt{-g}\left[
\frac{M_6^4}{2}(\mathcal R_6-2\Lambda_6)
-\frac{1}{2g_6^2}\operatorname{Tr}_{\mathbf8}(F_{MN}F^{MN})
+i\operatorname{Tr}(\bar\Psi\Gamma^M D_M\Psi)
-M_\Psi\operatorname{Tr}(\bar\Psi\Psi)\right],
\]

\[
F=dA-iA\wedge A,\qquad D_M\Psi=\nabla_M\Psi-i[A_M,\Psi].
\]

Here A is a U(8) connection and Ψ is one **complex Dirac** spinor in its adjoint. Tr is the fundamental matrix trace, with generators normalized by Tr(Ta Tb)=δab/2. The identity generator acts trivially on adjoint matter. The metric, M6, g6, Λ6, the optional Dirac mass, and the compactification data are inputs. In mass units, [M6]=[MΨ]=1, [Λ6]=2, [g6]=−1, [A]=1, and [Ψ]=5/2. The zero-mass control sets MΨ=0; it does not derive this choice.

On a fixed unwarped product with area A_X, the coefficients obey M_Pl,4²=M6⁴ A_X and 1/g4²=A_X/g6². These follow by integration and predict neither coefficient without the parent parameters and a radius. Writing the Einstein term makes the metric dynamical by assumption; it does not prove that this product solves the Einstein equations.

Decompose the rank-eight bundle as

\[
E=(\mathbb C^3\otimes L_c)\oplus(\mathbb C^2\otimes L_L)
\oplus L_a\oplus L_b\oplus L_d,
\quad m=(0,-3,3,3,0).
\]

For a matrix block A_ij, the representation is (N_i, conjugate N_j), the line bundle is L_i L_j^{-1}, and the Abelian charge is e_i−e_j. A useful block drawing is

\[
A_M=\begin{pmatrix}
A_c&X_{cL}&X_{ca}&X_{cb}&X_{cd}\\
X_{Lc}&A_L&H_u&H_d&X_{Ld}\\
X_{ac}&X_{aL}&A_a&X_{ab}&X_{ad}\\
X_{bc}&X_{bL}&X_{ba}&A_b&X_{bd}\\
X_{dc}&X_{dL}&X_{da}&X_{db}&A_d
\end{pmatrix}_M.
\]

The symbols Hu and Hd denote potential Higgs components only when M is internal. Hermiticity relates opposite gauge-field blocks. Complexified Lie-algebra dimensions are

\[
(N_iN_j)=\begin{pmatrix}9&6&3&3&3\\6&4&2&2&2\\3&2&1&1&1\\3&2&1&1&1\\3&2&1&1&1\end{pmatrix},
\]

whose entries sum to 64=16 diagonal +48 off-diagonal. The flux matrix is

\[
\mathcal M=\operatorname{diag}(0,0,0,-3,-3,3,3,0),
\quad\operatorname{Tr}\mathcal M=0,\quad\operatorname{Tr}\mathcal M^2=36.
\]

**Enhancement caveat.** This flux alone has centralizer U(4)×U(2)×U(2), not the requested five-block group. The c,d blocks have equal curvature, as do a,b. If their flat holonomies also agree, eight extra gauge generators remain. Inequivalent flat twists can distinguish equal-degree lines, but are additional choices. Such twists do not change the following index and Morse-index conclusions.

## 2. Charged Higgs interactions genuinely exist in the gauge action

The internal part of the gauge vertex is Tr(barΨ Γ^m [A_m,Ψ]). Matrix units obey

\[
E_{ij}E_{kl}=\delta_{jk}E_{il},\quad
[E_{jk},E_{ki}]=E_{ji},\quad
\operatorname{Tr}\{E_{ij}[E_{jk},E_{ki}]\}=1
\]

for three distinct indices. Thus the four block cycles c→L→a→c, c→L→b→c, d→L→a→d, and d→L→b→d have nonzero group coefficients. Their charges sum to zero. With y=(1/6,0,−1/2,1/2,−1/2), the internal La and Lb components have hypercharges +1/2 and −1/2.

The gamma matrix matters: an internal component of a six-dimensional vector is not an elementary six-dimensional Lorentz scalar. Consequently, the AXG-03 same-chirality scalar-bilinear obstruction does not forbid this vertex. The nonzero group coefficient does not establish a nonzero wavefunction overlap, a selected light Higgs, or the observed Yukawa matrices. The general gauge-vertex reduction is standard; see [Cremades, Ibáñez and Marchesano, equation (1.1) and Appendix A](https://arxiv.org/pdf/hep-th/0404229). Their explicit torus wavefunctions are not being used on X0(143).

## 3. A Dirac adjoint cancels anomalies by retaining mirrors

The degree matrix is

\[
\Delta_{ij}=m_i-m_j=
\begin{pmatrix}
0&3&-3&-3&0\\
-3&0&-6&-6&-3\\
3&6&0&0&3\\
3&6&0&0&3\\
0&3&-3&-3&0
\end{pmatrix}.
\]

Let S be a theta characteristic, deg(S)=g−1=12. Riemann–Roch gives

\[
\operatorname{ind}\bar\partial_{S\otimes L_iL_j^{-1}}
=(12+\Delta_{ij})+1-13=\Delta_{ij}.
\]

This is the number of chiral representation multiplets in the index; multiplying by N_iN_j counts individual gauge components. Actual kernel dimensions can exceed the index and depend on the spin structure and holomorphic bundles.

Write Ψ=Ψ+⊕Ψ− with opposite six-dimensional chirality in the **same** gauge bundle. Product chirality satisfies χ6=χ4 χX. Every fixed internal zero mode therefore supplies opposite four-dimensional chiralities from Ψ+ and Ψ−. Their index matrices are Δ and −Δ; their sum vanishes block by block. Smooth changes of the common gauge flux cannot project out just one partner.

The same paired representation cancels the opposite local six-dimensional anomaly polynomials. This is the ordinary vectorlike determinant control, not a UV-completion theorem. Giving a gauge-invariant Dirac mass pairs the modes; it does not leave a chiral Standard Model. Boundary projections, defects, asymmetric matter, or another mechanism would change the specified model.

## 4. The eighteen H1 classes acquire a physical operator—but the background is a saddle

For a constant-curvature split unitary connection, standard Yang–Mills theory on a closed curve has a negative normal space H1(ad_negative). Its real Morse index is twice that complex dimension. This is [Atiyah–Bott, Proposition 5.4 and equation (5.10), printed pages 556–559](https://www.uvm.edu/~cvincen1/files/teaching/spring2019-math382/atiyahbott.pdf).

For a map from a slope-μ_high block to a slope-μ_low block, δ=μ_high−μ_low>0. Its Hom line has degree −δ and hence h0=0. Riemann–Roch gives h1=g−1+δ per matrix entry. Grouping the blocks by slopes gives

| Source slope and rank | Target slope and rank | Complex negative dimension |
|---|---|---:|
| +3, rank 2 | 0, rank 4 | 2·4·(12+3)=120 |
| +3, rank 2 | −3, rank 2 | 2·2·(12+6)=72 |
| 0, rank 4 | −3, rank 2 | 4·2·(12+3)=120 |

Thus the complex negative dimension is **312** and the real Morse index is **624**. A second arithmetic calculation sums (12+m_r−m_s) over the positive roots of the explicit eight-by-eight diagonal flux and gives the same result.

In particular, Hom(L_a,L_L) and Hom(L_b,L_L) have degree −6, giving eighteen unstable doublet multiplets each: 36 complex components per sector. This is an actual Yang–Mills Hessian interpretation of the H1 dimension; it is not a minimally coupled scalar zero mode. The gauge one-form Hessian has a background-curvature term absent from the scalar Bochner Laplacian. There are also **180 complex colored negative directions**. Their existence prevents interpreting this saddle as a selected electroweak vacuum.

These are directions in the gauge-fixed transverse Hessian, not merely infinitesimal gauge transformations. The count assumes the standard YM action and holds for the fixed compact curve. Gravity backreaction, additional interactions, or a different compactification require a new stability analysis. A negative gauge subspace at a genuine full critical point would remain an obstruction to a positive full Hessian, but no such full critical solution is asserted here.

The total U(8) bundle has degree zero. Its split line degrees are not separately protected U(8) topological charges once all off-diagonal connections are allowed. A flat connection has lower YM energy in this total topological class. That comparison reinforces the local Hessian result without specifying an endpoint or decay dynamics.

## Outcome

This control makes two desired ingredients explicit: a dynamical metric term and gauge-origin charged Higgs vertices. It simultaneously exposes mirror fermions, gauge enhancement unless extra holonomies are supplied, and a 624-dimensional real unstable gauge subspace. It therefore demonstrates why enlarging the gauge group alone is insufficient. Any next version must specify how it preserves chirality and stabilizes the unwanted charged directions while retaining the useful Higgs interactions.

The accompanying script verifies **31 exact finite checks**, using only the Python standard library. The index and Hessian theorems are cited mathematical inputs; the script checks their explicit substitutions, not a new proof of those theorems.
