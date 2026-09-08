# Hecke-block geometry of the active/fixed split

Frozen 8 September 2026 before reading the new overlap and coupling observables.
This continues MTFT_Active_Module_v0.1.0 on the unchanged supplied release 0.26.2.

## Fixed objects

Use the previous direct C_ij active projector P (complex rank8, real rank16)
and Q=I-P, with the independent localization projector for verification.
Use the exact T2 primary factors x, x+2, g4, h6 to define rational arithmetic
projectors E_ell,E_old,E_q4,E_q6 of real ranks2,4,8,12. Construct these by
CRT polynomials in T2 and independently from the package's exact block bases.
Labels: ell=143a1; old=level11 oldspace; q4=quartic new orbit; q6=sextic new orbit.
The real rank8 q4 block is NOT the complex rank8 local active sector.

Primary arithmetic operators: T2,T3,T5,T7,W11,W13,W143,STAR.
Replication operators: T17,T19. No prime, triangle, or projector optimization.
All comparisons use the same Hodge-orthonormal frame. Exact arithmetic checks
will establish the algebraic block identities and which operators preserve
the blocks. The projector P remains numerical and period-dependent.

## Geometry

For each block b measure tr(P E_b), its complement d_b-tr(P E_b), and the
eigenvalues of U_b^T P U_b for an orthonormal real block basis U_b. These are
squared principal cosines. Report complete spectra, not only sums. Check
0<=values<=1 to numerical tolerance and J-pairing. Eigenvalues within1e-8 of
0 or1 are consistent with the corresponding intersection; values between
1e-8 and1e-6 from an endpoint are ambiguous; values farther away are resolved
mixed directions. These are numerical intersection diagnostics, not exact
dimension certificates. Record overlap sums16 and10.

Also measure the off-diagonal arithmetic-block components E_b P E_c and the
idempotence defect of the diagonal truncation sum_b E_b P E_b. This tests
whether P can be replaced by a sum of arithmetic-block restrictions.

## Coupling and attribution

For A, record the singular spectrum and rank of F^T A U, where U,F are
orthonormal bases of active and fixed spaces. Normalize singular values by
||A||_F: <=1e-9 is small, >=1e-7 is retained, intermediate values make the
rank ambiguous. A numerically zero coupling must not acquire a spurious rank
by normalization to its own roundoff norm.

When A preserves each E_b, let C=[P,A]. Its block components satisfy
E_b C E_c = E_b P E_c A_c - A_b E_b P E_c, A_b=E_b A E_b.
Record the sixteen squared Frobenius norms divided by ||A||_F^2. Their sum
is ||C||_F^2/||A||_F^2; the matrix gives an additive attribution in the
arithmetic coordinates. Report diagonal versus off-diagonal totals. This
does NOT mean A transports vectors between distinct invariant Hecke blocks;
it records the mismatch of A with the local projector P.

As a registered secondary description, decompose QAP=sum_b K_b with
K_b=Q A_b P and record the signed Gram table
G_bc=tr(K_b^T K_c)/||A||_F^2. Its entries sum to ||QAP||_F^2/||A||_F^2.
Individual ||K_b||^2 are not additive contributions: interference terms may
cancel. Report this explicitly if it occurs. No entry is a probability.

## Verification and limits

Check exact E_b identities, CRT/basis agreement, and arithmetic commutation.
Use the independent localization projector and freshly orthonormalized exact
block bases to verify overlaps/spectra. Transport the entire block calculation
through the other two registered coordinate seeds2026,31415 from the previous
study; compare scalar observables with seed143 at50-digit period-input precision.
All matrix arithmetic is double precision, despite the upstream precision.

Full-space block off-diagonal residuals and independent-route agreements are
diagnostic consistency checks at1e-9, not interval bounds. Source remains
unchanged. Do not rerun or repair the known ambiguous Lie closures; no new
claim here depends on them. No theta evaluation, Clifford search, physical
mass/charge assignment, or Standard Model identification is part of this run.
Record exploratory follow-ups separately with their motivation and status.
