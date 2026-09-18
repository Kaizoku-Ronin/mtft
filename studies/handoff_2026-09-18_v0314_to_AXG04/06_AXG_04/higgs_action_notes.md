# AXG-04: an elementary Higgs action and its exact tree-level flavor limit

This is an explicitly changed, nonsupersymmetric six-dimensional EFT candidate. It does not preserve the old all-bifundamental matter assignment or its flux. Anomaly completion, a scalar potential and a gravitational solution must be tested separately.

## 1. Field orientation and index convention

Keep stack order `(c,L,a,b,d)`, ranks `(3,2,1,1,1)`, and the original hypercharge vector `(1/6,0,-1/2,1/2,-1/2)`. For an unordered sector `ij`, use the representation `(N_i,N_j-bar)` and central charge `e_i-e_j`. Write the intrinsic 6D chirality as `s_ij=+1` or `-1`, in the index convention

\[
I_{ij}=s_{ij}(m_i-m_j).
\]

A negative index means an excess of left-handed four-dimensional fields in the conjugate representation.

| 6D sector | 6D chirality | New bundle degree | Signed 4D index | Left-handed field |
|---|---:|---:|---:|---|
| `cL` | + | +3 | +3 | Q |
| `ca` | − | +3 | −3 | u^c |
| `cb` | − | +3 | −3 | d^c |
| `Ld` | + | −3 | −3 | lepton doublet |
| `ad` | − | −3 | +3 | nu^c |
| `bd` | − | −3 | +3 | e^c |

These values use

\[
m'=(0,-3,-3,-3,0).
\]

There are no `La`, `Lb`, `cd` or `ab` Weyls in this six-field test. Extra spectator fields, if added later, change the action and must be counted separately.

At the original flux `(0,-3,+3,+3,0)`, the same chirality assignment instead gives indices `(3,3,3,-3,-3,-3)`. Four sectors have the wrong SM orientation. Solving all six desired signed-index equations forces

\[
m_L=m_a=m_b=m_c-3,\qquad m_d=m_c.
\]

The common shift can be set to `m_c=0` for this bifundamental calculation. In particular, the flux realignment is required by the combined chirality and family-sign choices; it is not merely a numerical convenience.

## 2. Explicit elementary scalar vertices

Add weak scalar doublets

\[
\Phi_u:\quad q_u=e_L-e_a,\qquad
\Phi_d:\quad q_d=e_L-e_b.
\]

They carry hypercharges `+1/2` and `-1/2`, respectively. A conventional gauge kinetic term, positive scalar kinetic terms and Weyl kinetic terms can be augmented by

\[
\begin{aligned}
\mathcal L_{Y,6}=-[&y_{6u}\,\overline{\Psi_{ca}^{-}}\,\Psi_{cL}^{+}\,\Phi_u
+y_{6d}\,\overline{\Psi_{cb}^{-}}\,\Psi_{cL}^{+}\,\Phi_d\\
&+y_{6\nu}\,\overline{\Psi_{Ld}^{+}}\,\Psi_{ad}^{-}\,\Phi_u
+y_{6e}\,\overline{\Psi_{Ld}^{+}}\,\Psi_{bd}^{-}\,\Phi_d
+\mathrm{h.c.}].
\end{aligned}
\]

Color indices contract between a fundamental and its barred partner. A weak anti-fundamental contracts with the fundamental Higgs; the usual SU(2) epsilon identification rewrites this in the familiar left-handed SM notation if desired. For example,

\[
-(e_c-e_a)+(e_c-e_L)+(e_L-e_a)=0
\]

checks the up vertex, and

\[
-(e_L-e_d)+(e_a-e_d)+(e_L-e_a)=0
\]

checks the neutrino vertex. All four gauge charge sums vanish.

The Lorentz contractions also exist. Since `bar(P_+ Psi)=bar(Psi)P_-`, opposite-chirality bilinears contain `P_+P_+` or `P_-P_-`, rather than the zero product `P_-P_+` of AXG-03. There is one 6D field per sector; the three families arise as internal modes, so these four `y_6f` are numbers, not arbitrary three-by-three family matrices. With canonical 6D fields, `[Psi]=5/2`, `[Phi]=2`, hence `[y_6f]=-1` in mass units. This is a local nonrenormalizable EFT action, not a UV completion.

The scalar potential `V_6(Phi_u,Phi_d,...)` is additional data. Setting quadratic masses and nonminimal curvature couplings to zero is the minimal kinetic test below, not a symmetry-protected prediction. The chosen mixed chiralities are not asserted to fit ordinary 6D `(1,0)` hypermultiplets; their fixed chirality is visible in [Park and Taylor, Table 1](https://arxiv.org/pdf/1110.5916).

## 3. Exact Higgs kinetic zero modes need matching connections

The new flux gives

\[
\deg(L_L L_a^{-1})=\deg(L_L L_b^{-1})=0.
\]

Degree zero alone is insufficient for a scalar zero mode. For a closed connected curve,

\[
\ker(D^\dagger D)=\{\text{parallel sections}\}.
\]

A nonzero parallel section exists exactly when the relevant unitary line connection has trivial holonomy. A convenient sufficient choice is to identify `L_L`, `L_a` and `L_b` as Hermitian line bundles **with their connections**. Then the ratio bundles have a unit parallel section `t_0`; their normalized scalar profiles are

\[
h_0=\frac{t_0}{\sqrt{\mathcal A}},\qquad
\int_X |h_0|^2\,d\mathrm{vol}=1.
\]

There is one complex scalar profile per elementary Higgs field, hence one 4D doublet per field. These are zero modes of the minimal kinetic operator; the full mass also depends on the unspecified scalar potential.

## 4. The arithmetic spin structure can supply exactly three families

Let `D=P_1+P_2+P_3` be any three of the four W13-fixed CM points, and choose

\[
L_c=L_d=\mathcal O,\qquad
L_L=L_a=L_b=\mathcal O(-D).
\]

Use the audited arithmetic theta characteristic `S_0` with `S_0^2=K_X`, `deg S_0=12` and `H^0(S_0)=span{1,u}`. At the four CM points, `u` takes `+i/sqrt(13)` twice and `-i/sqrt(13)` twice. Every triple therefore contains both signs. Its evaluation matrix has rank two; a representative is

\[
\begin{pmatrix}
1&i/\sqrt{13}\\
1&i/\sqrt{13}\\
1&-i/\sqrt{13}
\end{pmatrix}.
\]

Thus no nonzero section of `S_0` vanishes at all three points, so

\[
h^0(S_0(-D))=0.
\]

Serre duality and Riemann–Roch then give

\[
(h^0,h^1)(S_0(D))=(3,0),\qquad
(h^0,h^1)(S_0(-D))=(0,3).
\]

These dimensions establish exactly three chiral zero modes and no paired zero modes for each displayed sector under this specific realization. An arbitrary degree-three flux establishes only `h^0-h^1=3`; `(3+k,k)` is possible. For example, if the degree-three divisor lies in the zero divisor of a section of `S_0`, the negative twist can have a section. The purity assertion is consequently arithmetic data, not a consequence of the integer index alone.

## 5. Constant Higgs profiles make the tree-level flavor matrices degenerate

The line-connection identification has a cost. The quark sectors `cL`, `ca` and `cb` now share the same internal Dirac operator and zero-mode Hilbert space. Opposite 6D chirality supplies opposite 4D chirality while leaving the same internal chirality and mode functions available. The scalar bilinear therefore contracts the ordinary Hermitian inner product of those internal modes.

Choose one common orthonormal basis `f_k`, and let the displayed sector bases differ by arbitrary unitary matrices:

\[
\chi_{Q,i}=\sum_k f_k(U_Q)_{ki},\qquad
\chi_{u,j}=\sum_k f_k(U_u)_{kj}.
\]

The normalized overlap for the Lagrangian convention `bar(u_R) Y_u Q_L` is

\[
(Y_u)_{ji}
=y_{6u}\int_X \chi_{u,j}^{\dagger}\chi_{Q,i}h_0\,d\mathrm{vol}
=\frac{y_{6u}}{\sqrt{\mathcal A}}(U_u^\dagger U_Q)_{ji}.
\]

Consequently,

\[
\boxed{Y_u^\dagger Y_u=\frac{|y_{6u}|^2}{\mathcal A}I_3,\qquad
Y_d^\dagger Y_d=\frac{|y_{6d}|^2}{\mathcal A}I_3.}
\]

The analogous inner-product argument applies to the `Ld`, `ad` and `bd` spaces and gives equal singular values within each lepton sector as well. Identifying `L_c` and `L_d` is not necessary for this within-sector conclusion, although it is a simple sufficient choice for the exact three-family realization above. Independent `y_6f` coefficients can change relative scales between species but cannot generate a hierarchy inside any one species.

CKM can be chosen to be the identity after independent right-handed basis rotations. With exactly degenerate quark masses, however, a displayed CKM matrix is basis-undetermined and unobservable; this is **not** a unique prediction that the measured CKM matrix equals the identity. The invariant statement is the absence of physical tree-level flavor mixing and

\[
[Y_u^\dagger Y_u,Y_d^\dagger Y_d]=0,
\qquad \det[Y_u^\dagger Y_u,Y_d^\dagger Y_d]=0.
\]

The overlap method and canonical normalization are standard in dimensional reduction; see [Cremades, Ibáñez and Marchesano, sections 2 and 3.1.3](https://arxiv.org/pdf/hep-th/0404229). The constant-profile identity here is an exact Hilbert-space argument, not a torus spectrum imported onto the genus-13 curve.

## 6. Radius, scope and next missing interaction

The relation `Y_4=y_6/sqrt(Area)` is dimensionally consistent: `y_6` has units of length and `sqrt(Area)` has units of length. In the Einstein ansatz

\[
ds_6^2=\rho^{-2}ds_{4,E}^2+\rho^2ds_{X,0}^2,
\]

canonical fermion normalization gives the same result,

\[
Y_4=\frac{y_6}{\rho\sqrt{\mathcal A_0}}.
\]

This fixes neither the radius nor a Higgs expectation value. At tree level and within this minimal action, the exact connection matching that allows a massless constant elementary Higgs also erases within-sector mass hierarchy.

This is a conditional tree-level obstruction. Space-dependent interactions, nonconstant Higgs profiles, differing fermion operators, localized terms, higher-derivative terms or radiative corrections can change it. None is computed or excluded by the identity above. Simply changing orthonormal bases cannot change it. The full KK interactions need not enjoy a U(3) flavor symmetry, so the argument must not be promoted to an all-orders degeneracy theorem.

## 7. Separate control: three parent copies and an explicit unit-flux line

Another proposed parent replaces the five-stack quiver by `SU(3) x SU(2) x U(1)_h x U(1)_X`, with `h=6Y`. It has three copies of the physical SM six-dimensional representations `Q+, L+, U-, D-, E-, N-`. Every one of these 6D fields has `X=+1`. An elementary `H` carries `(h,X)=(3,0)`, and its conjugate weak doublet `Htilde=i sigma_2 H*` has `(-3,0)`. The ordinary scalar vertices

\[
\overline Q\,\widetilde H\,\mathbf y_{6u}U,
\quad \overline Q\,H\,\mathbf y_{6d}D,
\quad \overline L\,H\,\mathbf y_{6e}E,
\quad \overline L\,\widetilde H\,\mathbf y_{6\nu}N
\]

are gauge invariant and connect opposite 6D chiralities. The bold coefficients now have independent three-by-three parent-generation indices. This is a different theory, not a reinterpretation of the six-field quiver.

Unit `X` flux supplies Dirac index one per parent copy. The arithmetic spin structure requires care: `L_X=O(P)` does **not** produce a pure one-mode spectrum. Because `S_0` is basepoint-free and has two sections,

\[
h^0(S_0(-P))=1,
\qquad (h^0,h^1)(S_0(P))=(2,1).
\]

Three parent copies would consequently yield six desired modes and three opposite-chirality modes, with net index three.

There is, however, an explicit pure degree-one choice. Take three distinct CM points `P,Q,R` with `u(P)=+i/sqrt(13)` and `u(Q)=-i/sqrt(13)`, and set

\[
\boxed{L_X=\mathcal O(P+Q-R).}
\]

The proof proceeds in two small evaluation matrices. Evaluation at `R` is a nonzero map from the two-dimensional `H^0(S_0)` to a one-dimensional fiber. Thus `h^0(S_0(-R))=1`. Riemann–Roch and Serre duality give

\[
h^0(S_0(R))-h^0(S_0(-R))=1,
\qquad h^0(S_0(R))=2.
\]

The natural inclusion `H^0(S_0) -> H^0(S_0(R))` is therefore an isomorphism. Since `P,Q` differ from `R`, imposing vanishing at `P,Q` now gives the matrix

\[
\begin{pmatrix}1&i/\sqrt{13}\\1&-i/\sqrt{13}\end{pmatrix},
\qquad \det=-2i/\sqrt{13}\ne0.
\]

Consequently,

\[
h^0(S_0(R-P-Q))=0,
\qquad h^1(S_0L_X)=0,
\qquad \boxed{h^0(S_0L_X)=1.}
\]

This supplies exactly one zero mode per parent copy and therefore exactly three in each representation sector. The negative coefficient in the divisor `P+Q-R` is not a physical flux singularity: it presents a smooth holomorphic degree-one line, which admits a smooth constant-curvature HYM connection. No cusp singularity or localized source has been introduced by this notation.

With no internal hypercharge flux, a covariantly constant elementary Higgs, and the same `X` line for all fermions, the single internal mode still supplies only its normalized Hermitian overlap. Therefore

\[
\boxed{\mathbf Y_{4f}=\frac{\mathbf y_{6f}}{\sqrt{\mathcal A}}.}
\]

Unlike the six-field quiver, this parent can have unequal singular values and observable mixing: they can be inserted through its arbitrary parent-generation matrices. The curve has not predicted those matrices, and this is not a way to derive the observed hierarchy from the one internal wavefunction. Anomaly, tensor, global and gravitational consistency are separate tests.

The reproduction script checks exact input, gauge, chirality, index, rank, normalization and matrix identities; its current count is recorded in the generated JSON. These assertions do not independently construct the anomaly or gravity sectors.
