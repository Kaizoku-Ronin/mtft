# T2 memory approximation: frozen selection and held-out assessment

Frozen locally on 8 September 2026 before evaluating any reduced-memory
candidate. The preceding study already measured full coupling rank five,
smallest B singular value approximately 0.08997, and exact five-state memory.
Those facts motivate the experiment; they are not new blinded observations.

## Inputs and scope

Use fixed extracts of the completed MTFT_T2_Effective_Operator_v0.1.0 study.
Primary input T2_blocks.npz SHA256:
af01726b204310d959241020af12594b3687814eeb2cf29f85072fbabfd97f8b.
Independent full-real-frame input T2_independent_geometry.npz SHA256:
2990bd3f4ec917f8ae1a264da3688bc235388a716124d9e371ce1026a232237f.
EXTRACTION_PROVENANCE.json records the parent files and retained keys.
The parent supplied MTFT 0.26.2 source archive is unchanged, SHA256
35cfc00c32877bd6c12b16464ddc8190bcca9e2f054474f4ab0a6c605fe9da38.

The primary matrix is H=[[A,B],[B*,D]] with 8 active and 5 complementary
complex dimensions, Hermitian in the frozen frame. Preserve A exactly.
No reweighting, renormalization, fitted counterterm, operator search, or
physical time/energy identification is permitted. Recheck block reconstruction,
selfadjointness, and basis orthogonality at tolerance 1e-9. All measured
geometric quantities and error bounds remain floating-point DIAGNOSTIC.
The replay consumes frozen upstream geometry rather than reconstructing it.

## Three prescribed candidate families

For an orthonormal complement basis Qr (5xr), r=0,...,5, define
Br=B Qr, Dr=Qr*D Qr and Hr=[[A,Br],[Br*,Dr]]. This is a Hermitian Galerkin
reduction. For r=0, Hr=A. For r=5, the response must recover the full model.

1. B-SVD: Qr consists of the top r right singular vectors of B, ordered by
   descending singular value. No training-score search within this family.
2. Pole subsets: diagonalize D in increasing eigenvalue order. At each r,
   evaluate every subset of r eigenvectors on the training spectral grid and
   select the subset minimizing the maximum relative Frobenius response error
   across that grid. There are 32 subsets total across all six ranks. Scores
   within 1e-12 of the best are ties; break ties by lexicographic index tuple.
3. Response snapshots: for each training z, form
   X(z)=(zI-D)^-1 B*. Define M=mean_z X(z)X(z)*. Qr uses the top r eigenvectors
   of this positive semidefinite Gram matrix, in descending eigenvalue order.
   Uniform weight is used per training point; no held-out information enters M.

Record every training subset score and the chosen bases before evaluating
held-out errors. Record eigenvalues/singular values and their gaps. A cutoff
inside a degeneracy at relative gap <=1e-9 is ambiguous and must be flagged;
do not present that particular truncated subspace as uniquely selected.
The r=0 and r=5 repetitions across families are controls, not independent models.

## Training and held-out points

Training spectral grid: x=linspace(-3,3,61), eta in {0.05,0.2,0.7},
z=x+i eta: 183 points. Only this grid enters pole-subset selection and the
snapshot Gram matrix. Model dimension is not chosen using held-out data
and then described as independently validated.

After freezing all per-rank choices, evaluate the prescribed 18 family/rank
models on the disjoint held-out grid x=-2.995+0.01*j, j=0,...,599,
with the same three eta values: 1800 points. This is deterministic numerical
holdout, not independent physical data. The prior full model may already have
been studied nearby; none of these reduced-model errors has been used before.

At each z, compare G_r(z)=[(zI-Hr)^-1]_active with the full projected
G(z)=U*(zI-H)^-1U. The primary code uses full/reduced spectral decompositions;
the independent route uses Schur solves in a separately reconstructed frame.
Record relative Frobenius error, absolute operator-norm error, and per-eta maxima.

Held-out dynamics are never used to choose a basis. Use t=(j+1/2)*2pi/256,
j=0,...,255, plus {0,0.25,1,pi,2pi}, sorted and deduplicated: 261 points.
Compare full and reduced active propagators for exp(-itH) and exp(-itHr).
Record relative Frobenius amplitude error, absolute operator-norm error, and
absolute difference of average complementary squared-norm transfer.
Compute the average from 1-||U_PP||_F^2/8. Retain any roundoff excursions in
raw data rather than silently clipping them to [0,1].

## Fixed budgets and decision rule

Primary accuracy budget: BOTH maximum spectral response error <=0.01 and
maximum active-propagator error <=0.01, using relative Frobenius norms.
For each family, select the smallest r whose TRAINING spectral error meets
the 1% budget; freeze this r before opening held-out outputs. Declare that
selection validated only if BOTH held-out gates pass. No replacement r may
be selected using a failed holdout and called independently validated.

Also freeze analogous training-based choices for budgets 5% and 10%, which
are secondary descriptive tiers, not substitutes for the primary result.
If a training choice fails holdout, report the failure and retain all ranks
as a descriptive error curve. Any smallest rank that passes after examining
all held-out ranks is only a scan result, not a fresh held-out validation.

Report error by eta and over time prefixes t<=0.25, t<=1, and the full 2pi
horizon as prescribed diagnostics. A frequency-only or short-time pass does
not count as passing the joint primary gate. No claim of a continuous-domain
or all-time error bound follows from maxima on these finite grids.

## Why a small B singular value can be misleading

For discarded complement Qd, the full-to-decoupled embedding error contains
C=column_stack vertically(B Qd, Qr*D Qd). Record direct discarded coupling,
retained/discarded D mixing, and ||C||_2. The Hermitian embedding error E has
||E||_2=||C||_2. Conditional algebra gives absolute projected bounds
||G-G_r||_2 <= ||C||_2/eta^2 and
||U_PP-U_r,PP||_2 <= min(2, |t| ||C||_2).
Record these bounds and check them numerically without calling a float result
an interval certificate. They can be loose; do not use them to overwrite a
measured held-out error. No sampling-error probability or model-selection
p-value is claimed.

## Verification and deliverables

The independent route starts with the prior full real26 operator and
localization projector, extracts its complex branch, and obtains fresh
orthonormal sector bases. It independently reconstructs the three candidate
families and their training choices, using Schur solves to score response.
Compare basis-independent response errors, selected subset indices, and
ambient complement projectors (phase changes are allowed), with tolerance1e-9.
Full-rank controls must reproduce the full response and propagation to1e-9.
Selected r<5 approximations are judged against the stated percentage budgets,
not a numerical identity tolerance. Independent checks share upstream data.

Separate exact symbolic controls must establish Hermitian reduction,
embedding/bound identities, and a counterexample to selecting importance only
by instantaneous B coupling. Preserve failed budgets and model families.
Save this protocol, inputs and parent provenance, selection freeze, scripts,
pointwise errors, independent checks, proof, scientific figures, and replay.

Pre-execution notation clarification: the phrase "column_stack vertically"
above means VERTICAL block stacking, C=[B Qd; Qr*D Qd], implemented by
`np.vstack`. Its shape is (8+r) by (5-r). No measurement or selection rule
is changed by this clarification.
