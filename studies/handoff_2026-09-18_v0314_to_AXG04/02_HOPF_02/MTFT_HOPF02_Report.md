# MTFT HOPF-02 — Arithmetic Hopf geometry and what level 143 factors

Investigation date: 17 September 2026. Based on the supplied MTFT 0.31.4
candidate and the preceding spin-lift audit. This is an additive research
note, not a modification of the package or the study paper.

**Execution result:** 85 exact checks passed; zero failed. The full output is
in `results.json`. The initial generic symbolic rank routine was interrupted
for excessive runtime; the final run uses SymPy's exact rational-domain rank
routine, evaluates the same commutator system, and completes successfully.

## Findings

The existing half-form pencil supplies an explicit map from X₀(143) to a
two-sphere inside the quaternionic Hopf base S⁴. Pulling back the Hopf
connection along this map yields a reducible SU(2) connection: its two complex
line summands have degrees −12 and +12. This is an arithmetic construction of
a nontrivial U(1) reduction, conditional on the already selected theta
characteristic and its section basis. It does not establish a new SU(2)
physical gauge sector.

The original quartic Hecke block cannot carry a quaternionic structure
commuting with T₂. Doubling its multiplicity does permit an explicit such
structure. That positive control enlarges the state space and introduces
choices; it is not a reinterpretation of the original eight real modes.

The factorization 143 = 11·13 controls congruence and moduli data. It gives
P¹(Z/143) ≅ P¹(F₁₁) × P¹(F₁₃), a set of 12·14 = 168 elements, and maps from
X₀(143) to X₀(11) and X₀(13). It does not give a sphere factorization.

There is, however, a precise appearance of the two proposed spheres in a
different object: the Hecke split of real H¹ into the sextic block and its
complement has dimensions 12+14. Its positive Hodge unit sphere is the join
S¹¹ * S¹³ = S²⁵. The script verifies the two complementary projectors.

| Statement | Status and scope |
|---|---|
| The spin-compatible matrices generate Q₈ | Exact matrix identities; the theta-square correction is from the previous audit |
| The section pencil has degree 12 | Exact divisor calculation |
| Its Hopf connection reduces to U(1), with line degrees ±12 | Exact differential and topological calculation for the specified connection |
| The full SU(2) bundle over the curve is topologically trivial | Standard bundle classification over a surface; the reduction is nontrivial |
| Quaternionic structure commuting with T₂ on the original quartic block | Impossible by real eigenspace dimensions |
| A doubled algebraic block admits it | Exact matrices, with a constructed control metric |
| Level 143 decomposes by CRT | Exact enumeration and proof |
| S¹¹ and S¹³ arise as unit spheres of the 12+14 Hecke split | Exact complementary projectors; their join is the Hodge unit S²⁵ |
| S¹⁴³ is a Shimura variety associated to X₀(143) | No such identification; the sphere has odd real dimension |
| A dynamical spacetime or parent gauge theory has been derived | No |

## 1. Level, dimension, and genus are different labels

With standard sphere notation,

\[
S^{143}=\{x\in\mathbb R^{144}:\|x\|=1\}
\]

is a real 143-dimensional manifold. In contrast,

\[
\Gamma_0(143)=\left\{
\begin{pmatrix}a&b\\c&d\end{pmatrix}\in SL_2(\mathbb Z):143\mid c
\right\},
\quad
Y_0(143)=\Gamma_0(143)\backslash\mathbb H,
\]

and X₀(143) is the compactification obtained by adjoining the cusps. Here
\(\mathbb H\) means the complex upper half-plane. Below, \(\mathbb H_q\)
will mean Hamilton's quaternions, to keep these two uses separate.

The open modular curve is a one-complex-dimensional Shimura variety (more
precisely the relevant connected component in the Shimura description).
X₀(143) is its compactification. A Shimura variety is not a sphere assigned to
the integer in a curve's name. Its smooth complex locus has even real
dimension, so the odd-dimensional sphere S¹⁴³ cannot itself be such a complex
variety. See Milne [1].

The compact curve X₀(143) has complex dimension 1, real dimension 2, and genus
13. Its first real cohomology has dimension 26. These are different objects
and different dimension counts.

One further notation distinction: the compact quaternionic Shimura curve
sometimes denoted X¹⁴³ uses 143 as a quaternion-algebra discriminant. It is
still a complex curve, not S¹⁴³. An indefinite quaternion algebra ramified at
11 and 13 is split over the reals; it is not the definite Hamilton algebra
used in the compact Hopf construction. Sharing the word “quaternionic” does
not identify these geometries; see Voight [7, §1].

## 2. What really factors at 143

The Chinese remainder theorem gives

\[
\mathbb Z/143\mathbb Z\simeq\mathbb F_{11}\times\mathbb F_{13},
\qquad
(a,b)\longmapsto78a+66b\pmod{143}.
\]

Check the coefficients by hand:

\[
78\equiv(1,0),\qquad66\equiv(0,1)
\quad\text{modulo }(11,13).
\]

They also satisfy 78² ≡ 78, 66² ≡ 66, 78·66 ≡ 0, and 78+66 ≡ 1 modulo
143. This is an actual ring product decomposition.

A point of P¹(Z/143) is a primitive pair (c,d), up to multiplication by a
unit. Reduction at the two primes gives the bijection

\[
\boxed{\mathbb P^1(\mathbb Z/143)\simeq
\mathbb P^1(\mathbb F_{11})\times\mathbb P^1(\mathbb F_{13}).}
\]

For a field Fₚ, use representatives (0,1) and (1,t), t ∈ Fₚ. There are p+1
points. Hence the two finite factors have 12 and 14 elements, and their
product has 168. This is the same index set used by the package's Manin
model; compare the modular-symbol construction in Stein [3].

For row vectors the modular generators act by

\[
S=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
T=\begin{pmatrix}1&1\\0&1\end{pmatrix},
\qquad
(c,d)S=(d,-c),\quad(c,d)T=(c,c+d).
\]

The script records all 168 images of each generator and verifies that these
permutations agree with the package after reordering. It checks S² = 1 and
(ST)³ = 1 on the projective set.

There is a useful hand picture. Arrange the set as a 12-by-14 grid. At each
prime T fixes (0,1), while cycling all p points (1,t). The grid separates into
four T-orbits:

| Coordinate types | Number of points | T-cycle length |
|---|---:|---:|
| Fixed at both primes | 1 | 1 |
| Finite at 11, fixed at 13 | 11 | 11 |
| Fixed at 11, finite at 13 | 13 | 13 |
| Finite at both | 143 | 143 |

The last row is one orbit because 11 and 13 are coprime. These are precisely
the four cusp widths. Their sum is 168. The S and ST orbit counts are 84 and
56, so the compact surface has

\[
\chi=4-84+56=-24=2-2\cdot13.
\]

The factorization survives into moduli. A cyclic subgroup C₁₄₃ of an
elliptic curve has unique subgroups of orders 11 and 13, and conversely two
such subgroups on the *same elliptic curve* sum to a cyclic order-143
subgroup. This uses the elliptic-curve moduli interpretation in Milne [2].
It produces forgetful maps

\[
X_0(143)\longrightarrow X_0(11),\qquad
X_0(143)\longrightarrow X_0(13),
\]

of generic degrees 14 and 12. Both maps retain the underlying elliptic curve.
Thus the relevant product is over the common j-line, rather than an
unrestricted product of two independent elliptic-curve moduli. On the generic
locus this identifies the level-143 moduli with the fiber product over
X(1); compact coarse curves require normalization of the dominating
component, with ramification at special points and cusps taken into account.

\[
\begin{array}{ccc}
X_0(143)&\longrightarrow&X_0(11)\\
\downarrow&&\downarrow\\
X_0(13)&\longrightarrow&X(1).
\end{array}
\]

The genus calculation checks the distinction:

| Curve | Modular index | Genus |
|---|---:|---:|
| X₀(11) | 12 | 1 |
| X₀(13) | 14 | 0 |
| X₀(143) | 168 | 13 |

In particular, X₀(13)(C) is topologically S²; its level does not make it S¹³.

## 3. What operations on spheres actually do

| Operation | Result or ambient dimension |
|---|---|
| S¹¹ × S¹³ | A 24-dimensional product, not a sphere |
| S¹¹ * S¹³, the join | S²⁵ |
| R¹² ⊕ R¹⁴ | R²⁶, with unit sphere S²⁵ |
| R¹² ⊗ R¹⁴ | R¹⁶⁸, with unit sphere S¹⁶⁷ |
| R¹¹ ⊗ R¹³ | R¹⁴³, with unit sphere S¹⁴² |

The join identity follows from

\[
(v,w,\theta)\longmapsto(\cos\theta\,v,\sin\theta\,w),
\quad 0\leq\theta\leq\frac\pi2.
\]

At the endpoints one factor is collapsed. This has dimension 11+13+1 = 25.
The tensor map (v,w) ↦ v⊗w reaches only the rank-one unit tensors, not the
whole S¹⁶⁷; simultaneously changing both signs leaves its image unchanged.

MTFT's positive Hodge metric on H¹(X;R) really does define a unit S²⁵.
Here the **actual Hecke decomposition** supplies a distinguished 12+14 split:

\[
H^1(X;\mathbb R)=V_{12}\oplus V_{14},\quad
V_{12}=V_{\rm sextic},\quad
V_{14}=V_{\rm elliptic}\oplus V_{\rm old}\oplus V_{\rm quartic}.
\]

The dimensions are 12 and 2+4+8=14. Since the good-prime T₂ is self-adjoint
for the standard positive Hodge/Petersson inner product, its spectral blocks
with disjoint eigenvalues are orthogonal. Thus

\[
\boxed{S(H^1(X;\mathbb R))\cong S(V_{12})*S(V_{14})
\cong S^{11}*S^{13}\cong S^{25}.}
\]

The displayed join map realizes this identification with the Hodge norms on
the factors. The sphere factors are not individually assigned to the level
primes 11 and 13. They arise from the sextic block versus its complement.
The modular curve itself remains two-dimensional.

This split is reproducible without choosing eigenvectors. Write

\[
h=x^6-10x^4+2x^3+24x^2-7x-12,\qquad
r=x(x+2)(x^4-3x^3-x^2+5x+1).
\]

The extended Euclidean algorithm gives a(x)h(x)+b(x)r(x)=1, with

\[
b(x)=\frac{-2641x^5+984x^4+23014x^3-16562x^2-37008x+40927}{74424}.
\]

Then P₁₂=b(T₂)r(T₂) and P₁₄=I−P₁₂. The script verifies their idempotence,
complementarity, dimensions, and annihilator identities, and records both
26-by-26 rational matrices. The self-adjointness theorem identifies these
algebraic projectors with orthogonal Hodge projectors; no measured Gram
matrix was replaced by the control metric of §7. This post-preregistration
deduction is a consequence of the already recorded block decomposition.
The JSON matrices use the package's **homology** basis. Their transposes give
the induced operators in the dual **cohomology** basis; the dimensions,
annihilator polynomials, and sphere decomposition are unchanged.

Likewise the 168-dimensional permutation module on the CRT
index set has a genuine tensor decomposition

\[
\mathbb R[\mathbb P^1(\mathbb Z/143)]\simeq
\mathbb R[\mathbb P^1(\mathbb F_{11})]\otimes
\mathbb R[\mathbb P^1(\mathbb F_{13})].
\]

The modular generators act on the two factors simultaneously. Passing to
Manin relations, the relative homology quotient, and then its cuspidal part
is additional work; the finite tensor factorization is not a product
decomposition of H¹ or of X₀(143).

S¹⁴³ nevertheless has a legitimate generalized Hopf construction:

\[
\boxed{S^3\hookrightarrow S^{143}\longrightarrow\mathbb{HP}^{35}.}
\]

It is the unit sphere in quaternionic dimension 36: 4·36−1 = 143. The base
has real dimension 140, not 4. The quaternionic projective line alone gives
the special spherical base HP¹ = S⁴. The same S¹⁴³ is the unit sphere in C⁷²
and gives S¹ → S¹⁴³ → CP⁷¹. See the general bundle construction in Hatcher
[4, §4.2]. Neither construction selects level 143.

## 4. The arithmetic half-form pencil selects a Hopf reduction

Keep the earlier selected theta characteristic

\[
S_0=\mathcal O_X(D),\quad D=6c_0+6c_{11},\quad
S_0^{\otimes2}\simeq K_X,
\]

with sections 1 and u, where

\[
u=\frac{\eta(13\tau)\eta(143\tau)}{\eta(\tau)\eta(11\tau)}.
\]

For eta exponents r_d = (−1,−1,1,1), the Ligozat order calculation used in
the supplied arithmetic module is

\[
\operatorname{ord}_{c}(u)=\frac{143}{24\gcd(c,143/c)c}
\sum_{d\mid143}\frac{\gcd(c,d)^2r_d}{d}.
\]

In the cusp order (0,1/11,1/13,∞), it gives

\[
\operatorname{div}(u)=(-6,-6,6,6).
\]

The two sections have no common zero: the section 1 vanishes at the two
poles of u as a section of O(D), where u is nonzero in a local bundle frame;
u vanishes at the other two cusps, where 1 is nonzero. Therefore

\[
f:X\longrightarrow\mathbb{CP}^1,\qquad
f(p)=[1:z(p)],\qquad z=\sqrt{13}\,u
\]

is a globally defined holomorphic map of degree 12, and

\[
f^*\mathcal O_{\mathbb{CP}^1}(1)\simeq S_0.
\]

This map is specified once S₀ and the normalized section basis have been
selected. It is not a claim that the whole MTFT program uniquely forces that
selection.

The previously established unit identities imply

\[
z\circ W_{11}=-z,\qquad z\circ W_{13}=-1/z.
\]

In the column convention [1:z], their projective matrices are

\[
A=\begin{pmatrix}i&0\\0&-i\end{pmatrix},\qquad
B=\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]

Each squares to −I, while AB = −BA; their projective actions commute. Thus
the construction retains the Q₈ lift and V₄ base action from the spin audit.
The script rechecks the matrices and resulting rational transformations;
the eta transformation identities themselves are inherited exact inputs.

Embed CP¹ in HP¹ using complex coordinate pairs. Under

\[
h(q_1,q_2)=\left(2q_1\bar q_2,|q_1|^2-|q_2|^2\right),
\]

the normalized pair (1,z)/√(1+|z|²) maps into an S² ⊂ S⁴: the j and k
components of the quaternion coordinate vanish. This is an explicit
arithmetic route into the Hopf geometry, but its image is two-dimensional.

## 5. A distinction between two Q₈ actions

The Hopf structure group acts on the **right**:

\[
(q_1,q_2)\mapsto(q_1v,q_2v),\quad v\in Sp(1).
\]

This leaves h unchanged. Restricting to the right Q₈ action gives the valid
associated bundle S³/Q₈ → S⁷/Q₈ → S⁴ discussed previously. Its fiber is not
a quotient group: Q₈ is not normal in SU(2).

The arithmetic projective matrices above instead have a natural realization
acting on the **left**, mixing the two quaternionic coordinates. They are
elements of Sp(2), and generally move the base point. For example, A sends
(1,1)/√2 to (i,−i)/√2, changing the first Hopf coordinate from +1 to −1.
Right multiplication by i sends (1,1)/√2 to (i,i)/√2 and leaves that coordinate
at +1.

Thus the abstract inclusion Q₈ ⊂ SU(2) in the previous discussion remains
correct. The arithmetic realization does not automatically identify its
geometric action with the fiber action. This explicit distinction is needed
before interpreting either SU(2) as a physical gauge group.

## 6. Connection, curvature, and flux: formulas to reproduce

On the unit sphere in \(\mathbb H_q^2\), the Hopf connection is

\[
\mathcal A=\operatorname{Im}(\bar q_1\,dq_1+\bar q_2\,dq_2).
\]

Choose the local section (q,1)/√(1+|q|²), with
q = x₀+x₁i+x₂j+x₃k and D = 1+Σxᵤ². Its gauge potential is

\[
A=\frac{\operatorname{Im}(\bar q\,dq)}{D}.
\]

Use imaginary quaternion generators i,j,k with [i,j]=2k, and the orientation
dx₀∧dx₁∧dx₂∧dx₃. Computing F = dA + A∧A gives

\[
\begin{array}{lll}
F_{01}=2i/D^2,&F_{02}=2j/D^2,&F_{03}=2k/D^2,\\
F_{23}=-2i/D^2,&F_{13}=2j/D^2,&F_{12}=-2k/D^2.
\end{array}
\]

Here F = Σ_{μ<ν} F_{μν} dxμ∧dxν. These formulas give *F = −F. In the
complex two-dimensional matrix representation tr(i²)=tr(j²)=tr(k²)=−2,

\[
\operatorname{tr}(F\wedge F)=\frac{48}{(1+r^2)^4}\,d^4x.
\]

With anti-Hermitian curvature and c₂ = tr(F∧F)/(8π²),

\[
\int_{S^4}c_2=
\frac{48\cdot2\pi^2}{8\pi^2}
\int_0^\infty\frac{r^3\,dr}{(1+r^2)^4}
=12\cdot\frac1{12}=1.
\]

This checks the basic instanton described by Lévay [5]. Orientation changes
or opposite curvature conventions alter signed charges, not their magnitude.

On the complex plane q=x+iy, the potential and curvature reduce to

\[
A=i\frac{x\,dy-y\,dx}{1+x^2+y^2},\qquad
F=\frac{2i\,dx\wedge dy}{(1+x^2+y^2)^2}.
\]

Only one Lie-algebra direction remains, and A∧A = 0. This is the U(1)
monopole reduction. For c₁=iF/(2π), the two complex line summands have Chern
numbers −1 and +1 on CP¹. Pullback by the degree-12 map gives

\[
\boxed{E\simeq S_0^{-1}\oplus S_0,\qquad
\deg(S_0^{-1})=-12,\quad\deg(S_0)=12.}
\]

This is a smooth complex-bundle identification with the indicated natural
line reduction. The induced connection is reducible. The underlying SU(2)
bundle over the surface is topologically trivial: principal SU(2) bundles on
a two-dimensional CW complex have no nontrivial obstruction, since SU(2) is
connected and simply connected. Its chosen U(1) reduction can nonetheless
have nonzero degree. H⁴(X;Z)=0 also prevents retaining an instanton number on
the curve. Vanishing c₂ alone is not being used as the general classification
argument for arbitrary base spaces.

These ±12 degrees are not the selected M1 ±3 family fluxes. Relating those
bundles requires additional data. No family count or mass scale is inferred
from the instanton number.

There is a further geometric limitation, derived by Riemann–Hurwitz:

\[
\deg R_f=(2g_X-2)-12(2g_{\mathbb{CP}^1}-2)=24+24=48.
\]

The four order-six zeros and poles account for 4·5=20 ramification units;
28 remain away from those cusps, counted with multiplicity. Thus pulling
back the round metric of the image sphere produces a degenerate metric at
ramification points. It cannot simply replace the smooth Hodge or internal
hyperbolic metric. The pulled-back bundle connection remains smooth.

## 7. The Hecke obstruction and a positive control

The script reconstructs the package's T₂ on the 26-dimensional rational
cuspidal homology, rather than accepting a stored eigenvalue list. It checks

\[
\chi_{T_2}=x^2(x+2)^4g_4(x)^2h_6(x)^2,
\quad g_4=x^4-3x^3-x^2+5x+1,
\]

and restricts to ker g₄(T₂). On that real eight-dimensional block, g₄
annihilates the operator, is squarefree, and has four distinct real roots.
Consequently there are four real two-dimensional eigenspaces. The real
centralizer is M₂(R)⁴, of dimension 16; the script also computes that dimension
by solving the commutator equations.

Quaternionic generators commuting with T₂ would preserve every eigenspace.
On a real plane, set

\[
I=\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]

Solving IJ = −JI forces J = [[a,b],[b,−a]], so J²=(a²+b²)I₂. It cannot square
to −I₂ over the reals. This proves the obstruction without a nonlinear
numerical search. It concerns quaternionic actions commuting with T₂, not
all nonabelian symmetries, all quaternionic identifications, or all possible
physical transformations.

For a positive control, let

\[
C=\begin{pmatrix}0&0&0&-1\\1&0&0&-5\\0&1&0&1\\0&0&1&3\end{pmatrix},
\quad T_{16}=I_4\otimes C.
\]

This doubles the original multiplicity from two to four at each real
eigenvalue. On the four-dimensional multiplicity factor take

\[
L_i=\begin{pmatrix}0&-1&0&0\\1&0&0&0\\0&0&0&-1\\0&0&1&0\end{pmatrix},
\quad
L_j=\begin{pmatrix}0&0&-1&0\\0&0&0&1\\1&0&0&0\\0&-1&0&0\end{pmatrix}.
\]

Then I₁₆=Lᵢ⊗I₄ and J₁₆=Lⱼ⊗I₄ square to −I, anticommute, and commute with
T₁₆. They commute with every polynomial in that operator as well.

An exact positive control metric is I₄⊗G, where

\[
G_{ab}=\operatorname{tr}(C^{a+b}),\qquad 0\leq a,b\leq3.
\]

Explicitly,

\[
G=\begin{pmatrix}
4&3&11&21\\3&11&21&55\\11&21&55&128\\21&55&128&323
\end{pmatrix}.
\]

Its leading principal minors are 4, 35, 216, and 1957, all positive.

Writing the four real roots as λᵣ shows
Gₐᵦ=Σᵣλᵣᵃλᵣᵇ: it is a positive Vandermonde Gram matrix. The code verifies
Sylvester's criterion and CᵀG=GC. Both quaternion generators preserve this
metric. This constructs an algebraic extension, not a canonical arithmetic
selection, an integral lattice identification, or equality to the package's
Petersson metric.

## 8. What the continuation establishes

There is now an explicit connection between the selected arithmetic spin
pencil and Hopf gauge geometry. The resulting connection preserves a U(1)
reduction and has degree fixed by the pencil. The calculation also specifies
two obstacles to a stronger interpretation: the original quartic Hecke
block forbids a commuting quaternionic structure, and the natural map into
S⁴ has image only S² and is ramified.

The next stronger model would have to specify which extra quaternionic
directions or mode multiplicities it introduces, how their operators act,
and which arithmetic symmetries they preserve. Its field equations and
metric dynamics would then be independent tests. No new physical dimension
is inferred solely from a level, a representation dimension, or a sphere
bundle.

## Reproduction and provenance

Run `python investigate.py` with Python ≥3.10 and SymPy installed. The script
uses the unmodified `input/hecke_v0314.py` snapshot from the user-supplied
0.31.4 archive. Its MIT license is included. No network or installed MTFT
package is required. Results include the full rational operator, quartic
restriction basis, control matrices, curvature components, and modular
permutations. The JSON records the source hash and all assertion outcomes.

The reported equalities use rational/symbolic arithmetic. Real-root isolation
uses exact rational intervals. Symbolic integration is over the displayed
rational radial integrands. Topological conclusions are proved in the text;
they are not disguised as numerical tests.

**Inherited inputs:** the theta-square audit, h⁰(S₀)=2, and the two exact
Atkin–Lehner unit identities. The cusp orders, matrix relations, induced
projective actions, operator restriction, and curvature are checked here.
The Ligozat route is shared with the package and is not claimed as an
independent second proof of its divisor theorem.

## References

1. J. S. Milne, [Introduction to Shimura Varieties](https://www.jmilne.org/math/xnotes/svi.pdf), introduction and §§4–5. Modular curves and Shimura data.
2. J. S. Milne, [Modular Functions and Modular Forms](https://www.jmilne.org/math/CourseNotes/MF.pdf), §8, especially Lemma 8.5. The cyclic-subgroup moduli interpretation.
3. W. Stein, [Modular Forms: A Computational Approach](https://wstein.org/books/modform/modform/modular_symbols.html), modular-symbol and congruence-group constructions.
4. A. Hatcher, [Algebraic Topology](https://pi.math.cornell.edu/~hatcher/AT/AT.pdf), §4.2 and the quaternionic projective-space bundles.
5. P. Lévay, [The Geometry of Entanglement: Metrics, Connections and the Geometric Phase](https://arxiv.org/abs/quant-ph/0306115), §§III–IV. Quaternionic Hopf connection and instanton geometry.
6. Supplied MTFT 0.31.4 candidate, `hecke.py` and `surface/arithspin.py`; preceding `MTFT_v0314_Spin_and_Vacuum_Audit.md` and the September 16 study paper. Package-specific arithmetic inputs are distinguished from new deductions above.
7. J. Voight, [Shimura Curves of Genus at Most Two](https://arxiv.org/abs/0802.0911), §1. Quaternion-algebra discriminants, real splitting, and Shimura-curve notation.
