# MTFT: how the local active sector crosses the Hecke decomposition

**8 September 2026 · release 0.26.2 · continuation of the active-module study**

We now have a structural explanation for two previous observations. The
local active sector spans portions of every arithmetic component, while the
tested arithmetic operators preserve each component exactly. Its mixing
with the five local singlets therefore comes from the mismatch between
these two decompositions. W143's unusually small coupling has an exact
rank bound and an exact sign cancellation behind it.

The main findings are:

1. Exact rational projectors onto the four Hecke components agree by two
   constructions. Every registered operator preserves all four components.
2. The local active sector has substantial overlap with all four components,
   but contains none of them in full at the frozen diagnostic tolerance.
3. The largest W11 and W13 commutator contribution comes from the q4/q6
   component of the local projector. For W143 that contribution vanishes
   exactly because both arithmetic components have the same sign.
4. W143 has at most four real active-to-fixed coupling directions for any
   split. Our selected split reaches four numerically; the other tested
   non-STAR operators reach ten. STAR has zero resolved coupling.

## 1. Exact arithmetic components

The unchanged four-component decomposition is:

| Component | Meaning | Real dimension | Complex dimension |
|---|---|---:|---:|
| ell | elliptic factor 143a1 | 2 | 1 |
| old | level-11 oldspace | 4 | 2 |
| q4 | quartic newform orbit | 8 | 4 |
| q6 | sextic newform orbit | 12 | 6 |

In particular, the real eight-dimensional q4 component is different from
the complex eight-dimensional local active space.

We constructed rational idempotents E_b by the Chinese remainder theorem
applied to the squarefree T2 minimal factors x, x+2, g4 and h6. An independent
construction uses the exact block bases returned by the package. Exact
arithmetic verifies equality of these projectors, E_b E_c=delta_bc E_b,
and sum_b E_b=I, with ranks 2,4,8,12.

T2,T3,T5,T7,W11,W13,W143,STAR and replication operators T17,T19 commute
with every E_b exactly. Thus E_b A E_c=0 for b different from c. These
tested operators do not move vectors between distinct Hecke components.
Their algebraic properties are certified for the packaged matrices; the
original geometric identification of the integer-recognized W matrices is
inherited from MTFT.

## 2. Where the active eight lies

Let P be the local active orthogonal projector from the previous study,
in the same Hodge-orthonormal real coordinates. We measure
tr(P E_b)/2 in complex-dimension units and tr(P E_b)/dim_R(b) as a fraction
of each arithmetic component.

| Component | Active trace overlap, complex units | Active fraction of the component |
|---|---:|---:|
| ell | 0.9904595054 | 99.04595% |
| old | 1.4143810994 | 70.71905% |
| q4 | 2.4765873599 | 61.91468% |
| q6 | 3.1185720353 | 51.97620% |
| Total | 8.0000000000 | — |

These noninteger overlaps quantify alignment. They are not intersection
dimensions, probability assignments, particle counts, or fitted mass fractions.
The complementary overlaps sum to five complex units.

For finer geometry we computed the eigenvalues of U_b^T P U_b, where U_b
is an orthonormal basis of the arithmetic block. These are squared cosines
of principal angles, using the standard singular-value characterization.
[Zhu and Knyazev, *Angles between subspaces and their tangents*, Theorem 2.1](https://arxiv.org/pdf/1209.0523)

One value per complex direction is listed below; all values occur in real
pairs, as expected from complex-structure compatibility:

| Component | Squared principal cosines, ascending |
|---|---|
| ell | 0.990460 |
| old | 0.432015, 0.982367 |
| q4 | 0.217122, 0.486836, 0.873589, 0.899040 |
| q6 | 0.003925, 0.111393, 0.517798, 0.680737, 0.804719, 1.000000 |

At the frozen 1e-8 endpoint tolerance, only q6 has an active intersection:
one complex dimension. That minimum is already forced by 8+6-13=1 and is
not evidence of special arithmetic selection. No block has an intersection
with the fixed sector at this tolerance. No endpoint classification was
ambiguous. In particular, the approximately 99% elliptic overlap is still
not an exact elliptic component contained in the active sector.

Deleting the off-diagonal Hecke components of P does not produce another
projector. The resulting positive contraction sum_b E_b P E_b has
idempotence defect 0.79551758 in Frobenius norm. The discarded squared norm
is 3.24042693. Thus simply treating the local split separately inside each
Hecke component would lose part of its geometry.

## 3. What the coupling map means

The measured mismatch is C=[P,A]. Because A preserves the arithmetic
components, its commutator satisfies

\[
E_b C E_c=E_b P E_c A_c-A_b E_b P E_c,
\qquad A_b=E_b A E_b.
\]

The sixteen squared norms ||E_b C E_c||_F^2 are additive. They attribute the
failure to preserve the local active space to different pieces of P in
arithmetic coordinates. They must not be read as A transporting particles
or vectors between invariant Hecke components.

| Operator | Coupling rank, complex | Share of commutator squared norm in off-diagonal Hecke components |
|---|---:|---:|
| T2 | 5 | 78.27% |
| T3 | 5 | 62.75% |
| T5 | 5 | 55.73% |
| T7 | 5 | 48.01% |
| W11 | 5 | 100.00% |
| W13 | 5 | 99.07% |
| W143 | 2 | 82.54% |
| STAR | 0 resolved | undefined for vanishing coupling |
| T17, replication | 5 | 55.94% |
| T19, replication | 5 | 62.67% |

These are numerical coupling ranks under the registered singular-value bands.
The complete spectra and ordered 4x4 attribution matrices are in the JSON
ledger. The figure displays W11 and W143 on one common scale.

For T2, the q4/q6 pair accounts for 39.21% of the commutator squared norm,
old/q6 for 34.11%, and within-q6 for 16.10%. For W11, q4/q6 accounts for
74.57% and old/q6 for 25.14%. W13 has a similar 73.97% q4/q6 share. W143's
q4/q6 contribution is exactly absent for the algebraic reason below.

## 4. Why W143 couples less

The exact signs on the newform components are:

| Operator | ell | old | q4 | q6 |
|---|---|---|---|---|
| W11 | + | - | - | + |
| W13 | + | two real +, two real - | + | - |
| W143 | + | two real +, two real - | - | - |

If W acts by signs s_b and s_c on two components, then for any P,

\[
E_b[P,W]E_c=(s_c-s_b)E_bPE_c.
\]

The q4/q6 component contributes to W11 and W13 because their signs are
opposite. It contributes zero to W143 because both signs are -1. This is
an exact cancellation independent of the numerical local projector.

There is also an exact bound on the number of coupling directions. Let
R_+=(I+W143)/2. Exact arithmetic gives rank_R(R_+)=4: two real directions
in ell and two in old. For any idempotent P and Q=I-P,

\[
QW_{143}P=2QR_+P,\qquad \operatorname{rank}(QW_{143}P)\le4.
\]

The rank-four numerical result saturates this bound. The four retained
singular values, normalized by ||W143||_F, are 0.0677251 twice and
0.0370214 twice; all other singular values are below 7e-16. This accounts
for two complex coupling directions. A proof and the general involution
bound appear in `W143_COUPLING_LEMMA.md`.

The factorization through R_+ does not confine the coupling's endpoints to
ell and old: the local projectors P and Q themselves cross the arithmetic
components. The signs here are eigenvalues of an involution, not a spacetime
signature. The explanation gives mathematical structure without assigning
new particle or force labels.

## 5. Component amplitudes can cancel

A different decomposition writes QAP=sum_b K_b, with K_b=Q A_b P.
These maps need not be orthogonal. We therefore recorded their signed Gram
matrix G_bc=tr(K_b^T K_c)/||A||_F^2, whose full sum reconstructs the
active-to-fixed squared norm fraction.

| Operator | Diagonal Gram sum | Off-diagonal Gram sum | Net |
|---|---:|---:|---:|
| W11 | 0.1246318050 | +0.0974581822 | 0.2220899872 |
| W13 | 0.1267126294 | +0.0971688711 | 0.2238815005 |
| W143 | 0.1267126294 | -0.1147980849 | 0.0119145446 |
| STAR | 0.1246318050 | -0.1246318050 | consistent with zero |

Thus adding the individual component norms would substantially overstate
W143's coupling and would incorrectly assign coupling to STAR. The additive
commutator map in Section 3 and the signed amplitude decomposition answer
different questions; both reconstructions were checked.

For STAR, adding order-0.1 Gram entries leaves a tiny negative value of order
1e-17 in floating-point arithmetic, while direct squared coupling is of order
1e-28 on the independent route. That discrepancy is recorded using absolute
error. It is cancellation roundoff, not negative physical energy.

## 6. Verification and limits

Exact arithmetic checked every block projector, its agreement with the
independent basis construction, all ten operators' block preservation, and
the W143 rank/sign explanation. The geometric route used exact block bases
with independent QR and the earlier localization projector, rather than
CRT projectors and the direct channel kernel used by the main calculation.

The independent route passed 14 consistency checks. Across three coordinate
seeds, all overlaps, full principal spectra, coupling ranks, signed Gram
tables and commutator attributions agree within 1.79e-13 absolute difference.
Every discrete diagnostic count agrees. These routes share the release's
period/arithmetic inputs; they are computationally distinct, not independent
experimental data.

All period-dependent geometry remains float64 DIAGNOSTIC. Exact input bases
do not turn numerical intersections or rank saturation into exact theorems.
The previously recorded Lie-closure ambiguities were neither rerun nor
repaired; this study does not depend on them. No release source was modified.

## 7. Consequence for the next construction

The arithmetic decomposition and the local operator decomposition are both
useful, but they organize the same space differently. A future effective
operator must retain the measured coupling rather than identify the active
eight with an exact arithmetic component.

W143 now supplies a particularly small test case: at most two complex
coupling directions connect the local sectors. The next bounded construction
can reduce this involution to its coupled two-subspace blocks, retain the
complement through an exact block-resolvent identity, and check what an
eight-dimensional effective description preserves. This would be finite
operator mathematics; interpreting it as physical dynamics would require
additional structure. No such reduction or physical identification is
claimed in the present study.

The bundle contains the protocol, exact rational projectors, complete
ledgers, independent verification, figures, proofs, source archive, necessary
frozen inputs from the prior study, and a reproduction runner.
