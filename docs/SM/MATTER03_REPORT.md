# Matter-03: Higgs-branch existence, excess light modes, and normalization

The specified B3 construction admits a paired classical Higgs branch. The
single directed Higgs used for the Matter-02 residue witness cannot itself
solve the full compact D-term problem. The explicit pair leaves **237 charged
degree-one holomorphic deformation classes per electric-charge sign**. Its
common scaling remains unselected, and its physical kinetic metrics have not
been computed. These findings advance the vacuum gate and identify a concrete
family-selection obstruction; they do not predict observed particle masses.

This is an unreleased extension of the v0.34.0 source candidate. Released
functions, earlier archives, source data and version pins are preserved.

## What has changed since Matter-02

| Question | Matter-03 result | Scope |
|---|---|---|
| Can the single H_d be a nonzero vacuum? | No, on the unchanged full rank-eight background | Exact conditional stability obstruction |
| Can a reverse H_u repair that obstruction? | Yes, for the certified pair | Classical solution exists by the canonical-twisted surface correspondence |
| Is an expectation-value scale selected? | No | Nonzero common complex scaling preserves stability |
| Does this select three light charged families? | No: the explicit pair leaves 237 classes of each charge sign | Exact conditional holomorphic count |
| Could another aligned pair in these spaces do so? | At least 160 classes per sign remain | Dimension bound independent of the particular section profiles |
| Do normalized mass calculations have usable software? | Yes, with explicitly supplied inputs | Finite numerical linear algebra, no B3 metric or mass prediction |
| Is the physical photon established? | Q commutes with both Higgs arrows | Additional Abelian and axion/Stueckelberg gates remain open |

## The classical vacuum gate

Keep the Matter-02 holomorphic representative

\[
 S=\mathrm{Bl}_{(P,O)}(X_0(143)\times143a1),\quad
 \mathcal L=\pi^*\mathrm{pr}_C^*\mathcal O_C(P)(-3D),
\]

with P the ordinary width-one cusp at infinity and no elliptic flat twists.
The rank-eight bundle is

\[
 V=\mathcal O^{\oplus3}\oplus\mathcal L^{\oplus2}
   \oplus\mathcal L^2\oplus\mathcal L^{-3}\oplus\mathcal L^{-1}.
\]

At the ample wall mu(L)=A_E-3 epsilon=0 every summand has slope zero. The
defect-free classical bulk parent is the BHV canonical two-form twist,
Phi in H^0(End(V) tensor K_S). It is not an ordinary one-form Higgs bundle.

For a single map H_d:L^-1 to L^{oplus2} tensor K_S, the weak
subbundle is invariant and has the full bundle's slope. The equality case
of the integrated projection identity forces the incoming Higgs map to
vanish. This excludes every smooth positive Hermitian solution for nonzero
H_d, not just a frozen diagonal metric. Moving the polarization alone also
fails: the untouched L^2 and L^-3 spectators require 2 mu(L)<=0 and
-3 mu(L)<=0, so the wall is forced.

Choose one weak copy L+ and add the reverse arrow:

\[
 V_2=\mathcal L\oplus\mathcal L^{-1},\qquad
 \Phi_2=\begin{pmatrix}0&H_d\\H_u&0\end{pmatrix}.
\]

The Matter-02 H_d is a canonical differential of order three at P. The
certified H_u has order two. The product H_d H_u is a section of K_S^2 with
exceptional-divisor valuation 3+2+2=7. An odd divisor valuation cannot be a
meromorphic square, so the spectral polynomial lambda^2-H_d H_u is
irreducible over the surface function field. An invariant rank-one subsheaf
would give a generic eigenline and contradict that irreducibility. Thus the
active rank-two pair is stable; adjoining the zero-slope spectator lines
makes the complete pair polystable. Tanaka's surface-specific canonical-
twisted correspondence then supplies a Hermitian F/D solution.

This is an existence deduction under the specified classical assumptions.
No nonlinear metric equation has been numerically solved. Multiplying both
arrows by any nonzero t keeps the valuation odd and the pair stable, so no
physical Higgs scale is selected. The preserved compact classical gauge
centralizer is U(3) times U(1)^4; an SU(8) parent imposes the overall trace
constraint. The proposed Q=T3+Y commutes with both arrows. Extra Abelian
generators remain, and Higgs neutrality does not determine the axion kernel.

Detailed proof and rational representatives:
[Higgs vacuum derivation](HIGGS_VACUUM_DERIVATION.md).

## Recompute the light modes after choosing a Higgs

The split counts from Matter-02 cannot simply be reused in a nonzero Higgs
background. The new charged zero-mode gate is the two-term complex

\[
 \operatorname{End}V\xrightarrow{\operatorname{ad}\Phi}
 K_S\otimes\operatorname{End}V.
\]

Its degree-one hypercohomology is an extension of coker(d0) by ker(d1).
The complex has only two columns, so no later spectral-sequence differential
changes this count. Here d1 is the full holomorphic aa mass pairing, and d0
accounts for the charged gauge quotient. This includes charged Phi modes
that the isolated 18-by-44 residue block omits.

| Electric charge | Split a modes | Split Phi modes | d0 rank | d1 rank | Paired degree-one classes |
|---|---:|---:|---:|---:|---:|
| -1 | 216 | 63 | 1 | 41 | **237** |
| +1 | 217 | 62 | 1 | 41 | **237** |

The two charges have equal total split count 279, although their a and Phi
parts differ. The six independent rectangular aa blocks are:

| Spectator | Higgs | Block size | Rank of the certified pair | Maximum rank from size |
|---|---|---:|---:|---:|
| L- | down | 14 x 25 | 2 | 14 |
| L- | up | 32 x 14 | 2 | 14 |
| A1 | down | 18 x 44 | 5 | 18 |
| A1 | up | 55 x 15 | 3 | 15 |
| A2 | down | 72 x 32 | 17 | 32 |
| A2 | up | 25 x 87 | 12 | 25 |
| Total | | | **41** | **118** |

The exact rank proof uses projection to the curve, fat-point torsion and
the curve/elliptic contributions. It controls extension cross terms rather
than adding ranks of associated-graded maps without justification. Exact
finite rational matrices certify the local arithmetic; global duality and
Leray arguments supply the rest of the proof. The derivation was independently
reviewed before these ranks were frozen.

Even without the profile-specific rank calculation, the maximum 118 and
the one-dimensional charged gauge quotient give 279-118-1=160 remaining
classes per sign. This bound applies to aligned neutral two-arrow backgrounds
in the same fixed bundle and parent, with no a-field condensates or new
interactions. It is not a no-go theorem for every possible MTFT background.

The 237 classes are geometric infinitesimal deformations in the declared
bulk theory, not a catalogue of 237 experimentally observed particles.
Their physical interpretation needs the parent dictionary and the Hermitian
vacuum. They do obstruct identifying this restricted zero-mode spectrum
with only three charged families. Positive kinetic normalization preserves
the kernels; it cannot repair that obstruction.

Detailed proof: [charged deformation complex](CHARGED_COMPLEX_DERIVATION.md).

## Physical mass gate

The new field-coordinate normalization API implements

\[
 M_{\rm can}=\frac{\kappa v_H}{\sqrt2}
 C_L^T Y(C_H n) C_E,\qquad C_X^\dagger G_XC_X=I,
 \qquad n^\dagger n=1.
\]

It requires explicit positive Hermitian Grams, a tensor, a canonical Higgs
direction, a nonnegative expectation-value scale, and the overall parent
coefficient kappa. Complex changes of raw field coordinates preserve the
singular values when the same physical Higgs direction is transported.
The ordinary transpose in the holomorphic slot is essential. The historical
section-basis normalization API is unchanged.

Missing data include a selected scale, an acceptable light spectrum, actual
B3 surface mode metrics and parent kinetic factors, and the complete coupling
operator. A positive input metric example changes the rank-three witness's
singular values from (1,1,1) to (1,1/8,1/15). That is a diagnostic of dependence
on missing metrics, not a prediction.

See [physical inputs and Gemini vacuum clarification](MATTER03_PHYSICAL_INPUTS.md)
and [kinetic derivation](KINETIC_NORMALIZATION_REPORT.md). Standard Model
high-scale metastability is a different question from existence and scale
selection on this classical Higgs branch.

## Reproduce and inspect

```bash
python studies/matter03/run_all.py
python -m mtft.matter.higgs_vacuum --output higgs_vacuum.json
python -m mtft.matter.charged_modes --output charged_modes.json
```

The study includes independent certificates, public-API regressions,
normalization diagnostics, Node logic checks and a scientific figure. See
[the study guide](../../studies/matter03/README.md) and the retained
`studies/matter03/results/run_summary.json` for actual checks and timings.
This runner does not run the whole historical release suite or hours-long
studies. Node/mock-DOM tests do not verify browser layout.

Open [the offline Higgs lab](../../viz/matter03_higgs_lab.html) and
[scientific overview](../../viz/matter03_overview.png).

## Next investigation

Family selection is now a structural gate. A proposed change to the bundle,
additional condensates or interactions must recompute the full charged
complex and explicitly lift the excess modes while preserving an admissible
vacuum and photon. An arbitrary fitted soft mass term would change the model
without deriving that bridge. Metric computation and physical singular values
should follow an identified low-energy sector, with all supplied parameters
distinguished from arithmetic predictions.

Primary parent references:
[BHV, arXiv:0802.3391](https://arxiv.org/pdf/0802.3391), bulk equations and
couplings; [Tanaka, arXiv:1312.2673](https://arxiv.org/pdf/1312.2673),
Definitions 1.1-1.2 and Theorem 1.3 for canonical-twisted surface pairs.
The new obstruction, representatives, valuations and counts are deductions
for the frozen Matter-02 geometry, not results asserted by those papers.
