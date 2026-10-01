# Matter-02 photon gate: what the collinear background does and does not fix

1 October 2026 · MTFT v0.34.0 baseline · 38 exact checks PASS · physical photon OPEN

The simpler collinear background cancels the **linear** hypercharge flux projection but has a nonzero **quadratic** projection. This does not determine whether its photon is massive. It does determine exactly which candidate axion couplings would preserve or remove a hypercharge lift. The audit also finds that the mixed-anomaly kernel is too large to be called a fully anomaly-free gauge sector.

`check_photon.py` independently constructs the trace and anomaly ledgers, compares them with released arithmetic functions, and checks several explicitly supplied example mass matrices. The examples are INPUT actions for testing the missing gate; none is a derived B3 completion.

## Frozen data and the missing action

Stack order is `(c,L,A1,A2,B)`, ranks `(3,2,1,1,1)`, and line powers `(0,1,2,-3,-1)` of `v=(1,0,-3)`. The traceless hypercharge representative is

\[
Y=(1/6,0,-1/2,-1/2,1/2),\qquad v^2=-9.
\]

`surface/blowup.py` supplies the intersection lattice, canonical class, twist index, ampleness, slopes and split-background search. It does not supply an eight-dimensional form-field action with kinetic terms, gauge transformations and counterterms, or its four-dimensional reduction. `parent.py` catalogs product-surface records, not a complete B3 action. The existing `research/anomalies.py` and `charge_lattices.py` can analyze a **given** charge matrix; their M1 shift conventions cannot be transported to B3 as physical input without a new reduction.

The released ordinary adjoint calculation has `I10(adjoint U(8))=0`, so that particular local polynomial demands zero GS fields. This does **not** decide anomalies of the non-spin/spin^c twist, global anomalies, the full compactified sector, or extra fields required by a ultraviolet completion. In particular, it supplies no license to import the two GS fields required by the distinct B1 complex-representation construction.

## Exact trace projections

Writing `Tr` for the fundamental U(8) trace, and using normalized Chern-class variables,

\[
\operatorname{Tr}(Yc_1)=0,
\qquad \int_{\widetilde S}\operatorname{Tr}(Yc_1^2)=54.
\]

The underlying matrix moments `Tr(Y n^j)` for `j=0,1,2,3,4` are `(0,0,-6,9,-48)`. The last two numbers are **not** integrals of internal `c1³` and `c1⁴`: those differential forms vanish on a real four-dimensional surface. Mixed external/internal terms must be expanded before reducing. For example,

\[
\operatorname{Tr}(\bar F+Yf)^3\supset3f\operatorname{Tr}(Y\bar F^2),
\qquad
\operatorname{Tr}(\bar F+Yf)^4\supset6f^2\operatorname{Tr}(Y^2\bar F^2).
\]

Their internally integrated coefficients are `162` and `−189` in these normalized class conventions; factorials and a selected parent term can change the convention. A topological `f∧f` coefficient is not a photon mass or a positive energy density. The old score-31 background instead has quadratic projection `233`; both backgrounds have zero linear projection.

A common-center lift `Y(s)=Y+s I` leaves all adjoint root charges unchanged but changes other data:

\[
\operatorname{Tr}Y(s)=8s,\quad
\operatorname{Tr}(Y(s)c_1)=0,\quad
\int\operatorname{Tr}(Y(s)c_1^2)=54-144s,\quad
\operatorname{Tr}Y(s)^2=5/6+8s^2.
\]

Thus `s=3/8` cancels the quadratic projection while retaining the linear one. It changes the group normalization and is excluded by any independent trace shift requiring `TrY(s)=0`. This is why the common center must be included when selecting the physical lift.

The following are **conditional coupling probes**, not assumptions about MTFT's parent. Form degrees identify the relevant projection:

| Supplied eight-dimensional term and zero mode | Possible four-dimensional BF charge row, up to nonzero normalization |
|---|---|
| `B6 ∧ TrF`, with B6 having four internal indices and two external indices | `r=(3,2,1,1,1)` |
| `B4 ∧ TrF²`, with B4 having two internal and two external indices | `r n=(0,2,2,−3,−1)`, paired with each independent internal two-class |
| `B2 ∧ TrF³`, with an external B2 | `r n²=(0,2,4,9,1)` times `v²` |

Such forms might be absent, dual to one another, constrained or coupled differently. A parent and its duality constraints must decide which rows really coexist. A neutral axion multiplying `TrF⁴` alone is a theta interaction and does not imply a Stückelberg shift; a gauged axion kinetic term is separate data.

The linear row alone has rank one and contains Y. The linear and quadratic rows together have rank two and contain `Y+(3/8)I` rather than the displayed traceless Y. Adding the trace row gives rank three and excludes **every** `Y+sI` lift. These statements are exact for those INPUT rows, not photon-mass conclusions for B3.

## A new anomaly warning

For the conditional oriented adjoint net ledger `phi_i=−3 n_i`, define

\[
P(q)=\sum_{i<j}(\phi_i-\phi_j)r_i r_j(q_i-q_j)^3.
\]

The rows for SU(3)²U(1), SU(2)²U(1) and gravitational-U(1) anomalies are

\[
A=\begin{pmatrix}
0&-3&-3&9/2&3/2\\
9/2&-12&-3/2&6&3\\
0&-48&-48&72&24
\end{pmatrix}.
\]

Its rank is two and it annihilates Y and the center. However, writing its three-dimensional rational nullspace as `q=T t`, with T recorded in the JSON, gives

\[
P(Tt)=-72(t_0-t_1)^2(4t_0-t_1-3t_2).
\]

Therefore the whole mixed-anomaly kernel is **not** a fully cubic-anomaly-free gauge subspace. The generic repository description equating mixed-anomaly row kernels with all anomaly-free massless U(1)s needs this qualification in B3. The center/Y plane does have identically zero restricted cubic polynomial, including its polarized self/mixed cubic traces. Vector-like cohomological pairs do not change this net anomaly ledger.

Adding the illustrative row `(0,0,1,−1,0)` to two independent mixed-anomaly rows gives a rank-three matrix K3 whose kernel is precisely `span{Y,I}`. Adding the trace row gives the following rank-four INPUT matrix with kernel `span{Y}`:

\[
K_4=\begin{pmatrix}
0&2&2&-3&-1\\
3&-8&-1&4&2\\
0&0&1&-1&0\\
3&2&1&1&1
\end{pmatrix}.
\]

This shows an algebraic target a completion could achieve. Vanishing anomalies on its kernel is a necessary compatibility check; it is not the construction of quantized GS counterterms, periodic axions, an anomaly-free twist, or a full parent action. K4 was supplied for the demonstration.

## Electroweak breaking and positive kinetic normalization

In the neutral **connection** basis `(T3,c,L,A1,A2,B)`, the desired electromagnetic generator is `Q=(1,Y)`. The Q-neutral component of the `L/B` Hd doublet has charge covector

\[
q_{H_d}=(1/2,0,1,0,0,-1),\qquad q_{H_d}Q=0.
\]

Extend INPUT K4 by a zero T3 column and stack this Higgs row beneath it. The resulting constraint matrix has rank five and exactly one kernel direction, Q. A neutral `L/A1` Hu VEV also annihilates Q and does not increase its rank. This conditional example makes the missing photon gate concrete.

For positive gauge kinetic H and positive axion kinetic G, the mass form is `Kᵀ G K`, plus positive Higgs charge outer products. Its kernel is the intersection of the axion and Higgs kernels. An invertible canonical normalization changes coordinates and couplings, **never its rank or the existence of that kernel**. The script checks this exactly with a non-diagonal positive H. A supplied kinetic matrix cannot repair a charge direction already removed by an independent axion shift.

Matter-mode Gram matrices are a different issue. Within one unbroken electromagnetic charge sector, a constant generator `q I` remains `q I` after canonical normalization: normalizing the kinetic term also normalizes its minimal photon vertex. Such overlap matrices can change normalized Yukawa masses and field bases. They cannot turn an open photon-mass gate into a masslessness theorem.

## Next data required

1. Select an explicit B3 parent form-field and twist action, including gauge transformations and periodicities.
2. Reduce its independent fields and couplings on the blow-up, preserving form degree and duality constraints.
3. Check the full four-dimensional anomaly polynomial, including cubic anomalies and allowed counterterms.
4. Derive the shift matrix in the same neutral connection basis as the actual selected Higgs VEVs.
5. Verify the positive kinetic mass operator has the desired electromagnetic kernel, then compute its normalized coupling and threshold matching.

The conclusion remains **physical photon OPEN**. This investigation replaces a vague gate with exact projections, a warning about the anomaly kernel, and reproducible examples of both preserved and removed hypercharge lifts.

## Reproduce and primary references

```bash
python photon/check_photon.py --repo /path/to/mtft-repository --output photon/photon_checks.json
```

Only SymPy and the released MTFT source are required. Source hashes and all 38 exact checks are recorded.

- [Ruegg and Ruiz-Altaba, The Stückelberg Field](https://arxiv.org/abs/hep-th/0304245): gauge invariance can coexist with an Abelian vector mass, including a massive photon after electroweak mixing.
- [Beasley, Heckman and Vafa, GUTs and Exceptional Branes in F-theory I](https://arxiv.org/abs/0802.3391): the standard twisted surface gauge framework referenced by the B3 module; it does not provide MTFT's missing global completion.
- [Anderson et al., Hypercharge Flux in Heterotic Compactifications](https://arxiv.org/abs/1411.0034): a primary example of model-dependent massless-hypercharge conditions; its heterotic action is not adopted here.
