# Matter-01: one audited blow-up baseline and the electromagnetic bridge

1 October 2026 · MTFT v0.34.0 · Exact index/charge audit, physical completion OPEN

The B3 blow-up construction passes useful arithmetic gates, but it does not yet supply a physical electron–photon sector. This audit fixes one concrete background from the archived finite search, checks its slopes and net Standard Model charges, and records the missing steps. It also exposes a new selection issue: the assignment with the smallest **net-index extras score** can force very large cohomology spaces. That score should not select a physical vacuum by itself.

The accompanying `studies/vacuum_audit.py` runs **45 exact checks**, all passing, in roughly one second in the recorded environment. It writes `results/vacuum_audit.json`, with source hashes, witnesses, and eleven explicitly OPEN completion gates. The result does not certify a stable spacetime vacuum, full particle spectrum, measured mass, or derived electromagnetic coupling.

## 1. Sources and selection

The first lookup was the live Legend card `blowup_susy_twist`. The subsequent sources were `surface/blowup.py`, `research/anomaly8d.py`, PV02, the parent-action ledger, the previous independent parent audit, CC26, and ADS03. The generic `parent.gates()` interface handles product bidegrees; it is not a B3 blow-up validator. Its generic hypercharge-normalization verdict cannot be promoted to a result for this five-stack background.

PV02 and the earlier independent audit record **116 assignments** in the box `n_box=8, box_a=4, box_b=1, box_e=30`, with minimum extras score 31. These are exact finite-box records, not a global minimum theorem. Matter-01 did not rerun the full box. A small selection check at the recorded normal `(8,3,1)` returned four assignments and reproduced its first score-31 witness. All displayed witness arithmetic was then checked independently against the released functions.

The audit retains that witness as an archived baseline, **not an endorsed physical vacuum**. A simpler one-line-bundle construction is retained as a separate cohomology control.

## 2. Frozen baseline

The conditional parent is the topologically twisted eight-dimensional adjoint U(8) gauge construction on the blow-up of `X0(143) × 143a1`. Classes here are `(a,b,e)` with zero graph component; they mean `a f1 + b f2 + e E`. The canonical class is `24 f1 + E`, with `E²=-1`.

| Stack | Rank | Class `(a,b,e)` | Chirality potential `e−24b` | Slope at `(A_E,A_X,ε)=(8,3,1)` |
|---|---:|---:|---:|---:|
| c | 3 | (0,0,0) | 0 | 0 |
| L | 2 | (−3,1,21) | −3 | 0 |
| A1 | 1 | (−1,0,8) | 8 | 0 |
| A2 | 1 | (4,−1,−29) | −5 | 0 |
| B | 1 | (3,−1,−21) | 3 | 0 |

The Kähler point is strictly inside the stated product-polarization cone: `0<1<min(8,3)`. All slopes `a A_E+b A_X+e ε` vanish. The rank-weighted first Chern class also vanishes, so this is traceless flux data. These facts establish conditional split HYM/D-flatness; they do not compute the full fluctuation Hessian, select electroweak breaking, fix moduli, or prove the global parent consistent.

The net indices of Q, uᶜ, dᶜ, L, and eᶜ are each three. The index-forced lower bounds are five up-type pairs, five charged-lepton-singlet pairs, eight weak-doublet pairs, and thirteen neutral singlets. The published score 31 adds those different counting units. It is neither 31 particles nor the full spectrum. More cohomological pairs can occur.

## 3. New exact cohomology lower bounds

Riemann–Roch provides information even before constructing the individual cohomology spaces:

\[
\chi(L)=\frac{L^2-L\cdot\widetilde K}{2},\qquad
h^1(L)=h^0(L)+h^2(L)-\chi(L)\geq\max(0,-\chi(L)).
\]

For `(a,b,e)`, `L²=2ab−e²` and `L·K̃=24b−e`. This bound uses only nonnegativity of cohomology dimensions; it does not assume vanishing of h⁰ or h². It applies to any holomorphic representative with the specified class. Full dimensions, normalized bases, and their physical field dictionary remain uncomputed.

| Hom line | Class | χ(L) | χ(L⁻¹) | h¹(L) lower bound | h¹(L⁻¹) lower bound |
|---|---:|---:|---:|---:|---:|
| c/L | (3,−1,−21) | −222 | −225 | 222 | 225 |
| c/A1 | (1,0,−8) | −36 | −28 | 36 | 28 |
| c/A2 | (−4,1,29) | −422 | −427 | 422 | 427 |
| c/B | (−3,1,21) | −225 | −222 | 225 | 222 |
| L/A1 | (−2,1,13) | −92 | −81 | 92 | 81 |
| L/A2 | (−7,2,50) | −1263 | −1265 | 1263 | 1265 |
| L/B | (−6,2,42) | −897 | −891 | 897 | 891 |
| A1/A2 | (−5,1,37) | −683 | −696 | 683 | 696 |
| A1/B | (−4,1,29) | −422 | −427 | 422 | 427 |
| A2/B | (1,0,−8) | −36 | −28 | 36 | 28 |

For example, the weak-doublet line A1/L has class `(2,−1,−13)` and h¹≥81; it is **not** the charged-lepton singlet. The c/L net difference is only three while both orientations already have h¹ bounds in the hundreds. Thus three net family indices must not be relabeled as exactly three normalized family modes. Whether these cohomological modes correspond to particular massless fields, acquire couplings, or are lifted requires the complete twisted reduction and mass matrices.

The published collinear control uses `v=(1,0,−3)` and stack powers `(c,L,A1,A2,B)=(0,1,2,−3,−1)`. It is D-flat for `ε/A_E=1/3`, `A_X/A_E>1/3`. Its net-index score is 36 rather than 31, yet its c/L h¹ bounds are only 3 and 6, and its L/B bounds are 21 and 15. The unweighted sums of both-orientation h¹ lower bounds over the ten Hom lines are 666 for this control and 8588 for the score-31 baseline. These sums are bookkeeping comparisons, not particle counts or an optimized vacuum-selection functional. The simpler collinear family is a sensible next cohomology control; it is not established viable either.

## 4. Which U(1) becomes electromagnetism?

The starting reference QED engine uses **U(1) electromagnetism**. The proposed higher-energy MTFT construction instead assigns **U(1) hypercharge**, which combines with weak isospin after a neutral Higgs VEV:

\[
Q=T_3+Y.
\]

In this convention a traceless U(8) representative is

\[
(y_c,y_L,y_{A1},y_{A2},y_B)
=\left(\tfrac16,0,-\tfrac12,-\tfrac12,\tfrac12\right).
\]

Adjoint root differences yield Q:`1/6`, uᶜ:`−2/3`, dᶜ:`1/3`, L:`−1/2`, and eᶜ:`+1`. The lepton doublet becomes electric charges `(0,−1)`, while eᶜ is a left-handed anti-electron with charge +1. A charged Dirac electron requires pairing the electron component with the conjugate of eᶜ and supplying a nonzero normalized mass term. Neither this pairing nor its observed mass follows from a net index alone.

Hu-type and Hd-type doublets have Y=+1/2 and −1/2, respectively, so each has a Q-neutral component. The charged-lepton bundle triangles

\[
(A_i-L)+(B-A_i)+(L-B)=0
\]

also close exactly in class and hypercharge. This is a **necessary compatibility check**. It does not establish the required Dolbeault form degrees, zero-mode dimensions, a nonzero trilinear overlap, selected Higgs direction, kinetic Gram matrices, or a physical Yukawa matrix. Those are the next calculations.

For the displayed traceless representative `Tr₈(Y²)=5/6`, whereas `Tr₈(T₃²)=1/2`. A single universal parent Yang–Mills kinetic term, without thresholds or mixing corrections, would therefore give `gY²/g2²=3/5` and `sin²θW=3/8` **at its matching scale**. These are conditional group-normalization identities, not a prediction of the measured low-energy coupling.

There is a further freedom: adding the same number s to all five stack hypercharges leaves every adjoint root charge unchanged, while `Tr(Y(s)²)=5/6+8s²`. The actual extra-U(1), axion, and kinetic sectors must decide the physical lift. A preferred traceless representative is not a substitute for that calculation.

## 5. Photon masslessness is an additional gate

An Abelian gauge symmetry alone does not forbid a Stückelberg mass. An axion sector schematically contains `(∂a+K A)ᵀG(∂a+K A)`, giving a gauge mass matrix `KᵀG K`. The hypercharge direction must lie in its kernel in the same connection basis, and the selected electroweak VEV must leave Q unbroken. Canonical normalization must include the positive gauge kinetic matrix and any extra-U(1) mixing.

For the baseline the exact linear flux projection `Tr(Y c1)=0`, but `∫Tr(Y c1²)=233`. These are topological projections, **not** a supplied GS/axion mass matrix. The latter nonzero value is not by itself a photon-mass theorem, while the former zero value is not a proof of masslessness. The required higher-dimensional terms have not been selected and reduced.

With only the usual two electroweak neutral gauge fields and no extra mixing, the Higgs-generated mass matrix is

\[
M^2_{\rm EW}=\frac{v^2}{4}
\begin{pmatrix}g_2^2&-g_2g_Y\\-g_2g_Y&g_Y^2\end{pmatrix}.
\]

It has a null direction and gives `1/e²=1/g2²+1/gY²`. Adding an unspecified axion mass matrix can change that result. Matter-01 therefore leaves the physical photon kernel and gauge matching OPEN. The separate QED benchmark explicitly assumes an ordinary massless electromagnetic photon and inputs its coupling.

## 6. Anomalies and the remaining ledger

The net three-family Standard Model content has vanishing SU(3)³, SU(3)²Y, SU(2)²Y, Y³, and mixed gravitational–Y coefficients. Its twelve net SU(2) doublets satisfy the usual parity test. Vector-like pairs do not change those necessary anomaly checks. The parent module also returns local I₁₀=0 for an 8D adjoint. These exact checks do not establish global parent consistency or determine the other U(1) sectors.

| Completion gate | Matter-01 verdict |
|---|---|
| A specified arithmetic blow-up point and actual holomorphic bundles | OPEN |
| Full h⁰,h¹,h² and normalized mode spaces | OPEN; new exact RR lower bounds supplied |
| Nonzero normalized electron Yukawa and selected Higgs VEV | OPEN |
| Gauge-invariant masses lifting unwanted vector-like matter | OPEN |
| Hypercharge in every axion/Stückelberg kernel | OPEN |
| Physical gauge kinetic matrix, mixing, thresholds and RG matching | OPEN |
| Electroweak vacuum and full electromagnetic null direction | OPEN |
| Global anomalies and globally consistent twisted parent | OPEN |
| Complete stability and moduli stabilization | OPEN |
| Compactification, Higgs, gauge and electron scales | OPEN |
| A common parent linking B3 to C3X/ADS03 | OPEN |

The last row matters: ADS03 is a conditional classical Higgs–breathing calculation descended from the distinct six-dimensional C3X benchmark. Its radius and field parameters cannot be inserted into the eight-dimensional blow-up construction as derived matching data. The oloid/triangle identity still describes that scalar potential exactly; it supplies no missing fermion zero modes or electromagnetic coupling here. Feigenbaum scaling was not established in ADS03.

## 7. Reproduction and next use

From the candidate repository root, with dependencies installed:

```bash
python studies/matter01/studies/vacuum_audit.py
```

The script writes its source hashes alongside all exact gates. The distributed record was evaluated against the unchanged v0.34.0 attachment; a rerun against a separately patched candidate may update hashes of files not used in the arithmetic. No long enumeration is needed for these witness checks.

The next MTFT-specific calculation should choose an actual blow-up point and holomorphic bundle representatives, then compute cohomology and normalized lepton/Higgs Yukawas for the simpler collinear control. The immediate computational QED benchmark may proceed independently with explicitly supplied parameters. Success of that benchmark tests the matter engine; it does not close this vacuum ledger.

Primary physics references: [PDG Standard Model review index](https://pdgweb.lbl.gov/2025/reviews/standard_model_and_related.html); [A. Pich, Standard Model of Electroweak Interactions](https://cds.cern.ch/record/819632/files/0502010.pdf); [Ruegg and Ruiz-Altaba, The Stückelberg Field](https://arxiv.org/abs/hep-th/0304245). [Hypercharge Flux in Heterotic Compactifications](https://arxiv.org/abs/1411.0034) supplies an example of model-dependent axion constraints; its parent theory is not assumed to be MTFT's.
