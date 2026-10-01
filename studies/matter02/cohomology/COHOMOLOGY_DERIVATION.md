# Matter-02: exact cohomology of the collinear control

This calculation fixes an actual holomorphic representative, rather than inferring cohomology dimensions from its first Chern class. It is an exact algebraic result conditional on the archived modular-form basis being the stated basis of the full cusp-form space. No metric, period matrix, mass, or physical field interpretation enters the calculation.

## Geometry and conventions

Work in characteristic zero, initially over Q and then after base change to C. Let

- C = X0(143), the smooth compact modular curve of genus13;
- T = 143a1, with rational identity Q;
- P = infinity, the rational width1 cusp of C;
- Y = C×T and pi:S→Y the blow-up at p=(P,Q);
- E_exc the exceptional divisor;
- M=pr_C*O_C(P), L=pi*M⊗O_S(-3E_exc).

There are no additional degree-zero twists. The elliptic identity is enough to specify the point; no affine Weierstrass coordinates for Q are needed. In the repository's (a,b,e) convention, c1(L)=(1,0,-3). The five stack powers are (c,L,A1,A2,B)=(0,1,2,-3,-1). Every off-diagonal Hom line is a power L^q with |q|≤5.

The intersection and canonical classes are

\[
 c_1(L)^2=-9,\qquad K_S=\pi^*K_Y+E_{\rm exc},\qquad c_1(L)K_S=3.
\]

Since chi(O_Y)=(1−13)(1−1)=0 and blow-up preserves chi(O), surface Riemann–Roch gives

\[
 \boxed{\chi(L^q)=\frac{-9q^2-3q}{2}.}
\]

This sign also gives chi(L^q)−chi(L^(-q))=−3q, the repository's recorded net index. A later physical field dictionary must specify its representation orientation separately.

## The cusp is ordinary: an exact finite certificate

The frozen input contains the 141×13 matrix of a basis of S2(Gamma0(143)), indexed by qexpansion degree0 through140. Its first13 nonzero-degree rows form a square matrix with determinant

\[
 -18196056815164784640\ne0.
\]

The adapted basis is related to the original by the archived integral change-of-basis matrix of determinant−1078272. The verifier checks that transformation in every one of the141 coefficient rows and matches the provenance SHA1 of each input file.

At the width1 cusp, a weight2 cusp form f(z)=sum a_n z^n determines a regular canonical differential f(z) dz/z. Vanishing to order at least m therefore means a1=...=am=0. The rank of these first m coefficient rows is exactly m, for m=1,...,13. Both archived bases give that rank, with nonzero exact rational minors and explicit rational kernels included in the result file.

Consequently the canonical vanishing sequence at P is (0,1,...,12), and the Weierstrass gaps are (1,2,...,13). Thus

\[
 h^0_C(K_C(-mP))=\max(13-m,0),\qquad m\ge0,
\]

and

\[
 h^0_C(O_C(nP))=\max(1,n-12),\qquad n\ge0.
\]

Only coefficient rows up to9 are necessary for |q|≤5. The extra full-rank13 certificate allows an exact formula for every integer q. This is a statement at this chosen cusp; it does not assume every point on C has the same gap sequence.

## Blow-up pushforwards and local jets

For a point blow-up of a smooth surface,

\[
 \pi_*O_S(-kE_{\rm exc})=I_p^k,\quad R^1\pi_*O_S(-kE_{\rm exc})=0,\qquad k\ge0.
\]

One obtains these statements by induction from the ideal of E_exc and the restriction O_S(-E_exc)|E_exc=O_P1(1). The successive ideal quotients downstairs are the degree k homogeneous polynomials in two local coordinates. The required higher cohomology on P1 vanishes.

For positive exceptional powers, pi_*O_S(kE_exc)=O_Y. The exact sequences

\[
 0\to O_S((k-1)E_{\rm exc})\to O_S(kE_{\rm exc})
 \to O_{\mathbb P^1}(-k)\to0
\]

show that R1pi_*O_S(kE_exc) is supported at p and has length

\[
 \sum_{j=1}^k(j-1)=k(k-1)/2.
\]

Locally its principal-part basis may be represented by z^(-i)w^(-j), with i,j≥1 and i+j≤k. The fat-point quotient O_Y/I_p^k has length k(k+1)/2 and monomial basis z^i w^j with i,j≥0, i+j<k.

## Zeroth cohomology

For q>0,

\[
 H^0_S(L^q)=H^0_Y(M^q\otimes I_p^{3q}).
\]

All base sections are pulled from C because H0(T,O_T)=C. A nonzero section of O_C(qP) has an effective zero divisor of total degree q, and hence cannot vanish to order3q at P. Therefore h0(L^q)=0.

For q=−r<0, positive exceptional pushforward gives

\[
 H^0_S(L^{-r})=H^0_Y(M^{-r})=0
\]

because O_C(-rP) has negative degree. At q0, h0(O_S)=1.

## Second cohomology by Serre duality

For q>0,

\[
 K_S\otimes L^{-q}
 =\pi^*\bigl(pr_C^*(K_C(-qP))\bigr)
 \otimes O_S((3q+1)E_{\rm exc}).
\]

Positive exceptional pushforward yields

\[
 \boxed{h^2_S(L^q)=h^0_C(K_C(-qP))=\max(13-q,0).}
\]

For q=−r<0,

\[
 K_S\otimes L^r
 =\pi^*\bigl(pr_C^*(K_C(rP))\bigr)
 \otimes O_S(-(3r-1)E_{\rm exc}).
\]

A section of K_C(rP) pulled to Y depends only on z. Requiring its coefficient to lie in I_p^(3r−1) forces canonical vanishing of order(3r−1)−r=2r−1 at P. Thus

\[
 \boxed{h^2_S(L^{-r})=h^0_C(K_C(-(2r-1)P))
 =\max(14-2r,0).}
\]

There is no fiber-negativity argument that eliminates these spaces. The restriction to a generic elliptic fiber is trivial, and the genus13 curve supplies surviving canonical sections.

At q0, Kunneth and birational invariance give (h0,h1,h2)=(1,14,13).

## First cohomology and a second exact route

With h0=0 for q≠0, Riemann–Roch determines

\[
 \boxed{h^1_S(L^q)=h^2_S(L^q)+(9q^2+3q)/2.}
\]

For positive q, the fat-point exact sequence

\[
 0\to M^q I_p^{3q}\to M^q\to M^q|_{3q\,p}\to0
\]

provides a separate dimension check. Write a=h0_C(O(qP)), b=h1_C(O(qP)). Kunneth gives base dimensions(a,a+b,b). The global-section-to-fat-jet map is injective by the same degree argument. If N=3q(3q+1)/2, the long exact sequence gives h1=b+N and h2=b, in agreement with the displayed formula.

For q=−r, the positive exceptional Leray sequence has the terms

\[
 0\to H^1_Y(M^{-r})\to H^1_S(L^{-r})
 \to H^0_Y(R^1\pi_*L^{-r})
 \xrightarrow{\delta}H^2_Y(M^{-r})\to H^2_S(L^{-r})\to0.
\]

The two base spaces have dimension r+12. The torsion space has length3r(3r−1)/2. Serre duality identifies the rank of delta with (r+12)−h2_S(L^(-r)), determined by the canonical cusp jets above. This gives the same h1 independently of the surface intersection calculation. The verifier checks both routes for all Hom powers and checks the closed arithmetic formula on q=−25,...,25.

## Exact dimension table

| q | h0 | h1 | h2 | chi |
|---:|---:|---:|---:|---:|
| −5 | 0 |109|4|−105|
| −4 | 0 |72|6|−66|
| −3 | 0 |44|8|−36|
| −2 | 0 |25|10|−15|
| −1 | 0 |15|12|−3|
| 0 | 1 |14|13|0|
| 1 | 0 |18|12|−6|
| 2 | 0 |32|11|−21|
| 3 | 0 |55|10|−45|
| 4 | 0 |87|9|−78|
| 5 | 0 |128|8|−120|

The unweighted sum of H1 dimensions over all20 oriented off-diagonal Hom lines is864, compared with the Matter-01 RR lower bound666. The corresponding H2 sum is198. These are vector-space bookkeeping totals, without representation dimensions or a physical field dictionary. They are not particle counts.

In particular the color/weak Hom powers±1 have H1 dimensions15 and18, while the weak/B powers±2 have dimensions25 and32. These are exact spaces for this holomorphic representative before any couplings or symmetry breaking; they cannot be reduced to the net chiral index alone.

## Independent check of the local Yukawa witness

For q=(1,−3,2), the proposed pairing is

\[
 H^1(L)\otimes H^1(L^{-3})\otimes H^0(K_S L^2)\to H^2(K_S)\cong\mathbb C.
\]

Canonical order3 at P exists because the cusp jet ranks grow from3 to4. A section phi2 of K_S L² is represented downstairs by such a canonical differential, whose coefficient relative to K_C⊗O(2P) starts with z5 after algebraic leading-coefficient normalization. The H1(L) class from the fat jet w² modulo I³ survives because global M sections have no w dependence.

Multiplication by phi2 maps it to z5 w² modulo I⁸ in the fat-point model for H1(K_S L³). This total-degree7 monomial survives the quotient by I⁸ and by all global base sections, which again have no w dependence. The product is therefore nonzero. Perfect Serre duality guarantees a class in H1(L^(-3)) that pairs nontrivially with it.

The explicit exceptional principal part z^(-6)w^(-3) lies in R1pi_*O(9E_exc), since6+3=9. Its global obstruction pairs only with w-independent base canonical coefficients and hence vanishes. Its local residue against z5w² dz dw is nonzero. Coordinate leading-coefficient normalization here is algebraic; it is not normalization in an L² metric and does not determine a physical Yukawa value.

The stronger rank witness uses positive jets (w,w²,zw) and negative principal parts (z^(-6)w^(-2),z^(-7)w^(-2),z^(-6)w^(-3)). All three positive classes survive the global-section quotient, and the three negative classes lie independently in the kernel of the Leray obstruction because their w poles have order at least2. For phi2 coefficient v5z5+v6z6+... the residue submatrix is

\[
 \begin{pmatrix}v_5&v_6&0\\0&0&v_5\\0&v_5&0\end{pmatrix},
 \qquad\det=-v_5^3\ne0.
\]

The same conclusion follows directly because multiplication sends the positive jets to three independent w-dependent target fat jets. Hence the holomorphic pairing has rank at least3 in this chosen sector. It is not a complete mass matrix and does not identify a physically selected set of three families. The independent certificate constructs an actual rational canonical form of exact order3 from the archived basis, scales its a4 to1, and verifies the residue matrix and determinant−1 exactly.

## Verification and provenance

Run `python cohomology_certificate.py`. The script uses only the Python standard library and exact Fraction arithmetic. It writes a JSON result with148 passing gates, nonzero minors, explicit rational cusp-kernel vectors, source hashes, all Hom orientations and the dimension table. This count records finite arithmetic checks of the provided proof; it does not replace the geometric assumptions or independently rederive the full modular-form input from the definition of Gamma0(143).

The archived input files and their provenance are under `data/`. The claim assumes they are the correct full13-dimensional weight2 cusp-form basis as documented by the existing repository. No numerical metric or period data are introduced.

Supporting mathematical sources:

- Patrick Brosnan, [Algebraic Surfaces, §4.2](https://terpconnect.umd.edu/~pbrosnan/notes/Surfaces/sect0028.html): exceptional curve, canonical divisor and higher direct images of the structure sheaf. The ideal-power and positive-exceptional formulas above are derived from the displayed exact sequences.
- The Stacks Project, [Blowing up and de Rham cohomology](https://stacks.math.columbia.edu/tag/0FUC): birational/Hodge behavior under a smooth blow-up.
- William Fulton, [Algebraic Curves](https://www.math.lsa.umich.edu/~wfulton/CurveBook.pdf): Riemann–Roch, canonical differentials and the gap sequence.

The remaining physics questions include the twisted field dictionary, normalized bases, actual Yukawa matrices and their ranks, Higgs expectation values and paired-mode lifting, and the gauge kinetic and massless-photon matching problem.
