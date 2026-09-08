# A rank bound and sign cancellation for W143

8 September 2026. Explanatory follow-up to the registered coupling-rank
measurement. Elementary linear algebra; no novelty or priority claim.

## General involution bound

Let W be an involution on a finite-dimensional vector space over a field of
characteristic different from two. Its sign projectors are
R_+=(I+W)/2 and R_-=(I-W)/2. Let P be any idempotent and Q=I-P.
No orthogonality assumption is needed for

\[
QWP=2QR_+P=-2QR_-P.
\]

Indeed, W=2R_+-I=I-2R_- and QP=0. The rank of a product cannot exceed the
rank of any factor, so

\[
\operatorname{rank}(QWP)\le
\min\{\operatorname{rank}R_+,\operatorname{rank}R_-,
       \operatorname{rank}P,\operatorname{rank}Q\}.
\]

The bound is independent of how P was constructed.

## Application to the packaged W143

Exact rational checks on MTFT 0.26.2 give R_+^2=R_+ and rank_R(R_+)=4.
Its ranks in the four arithmetic blocks are

| Block | ell | old | q4 | q6 |
|---|---:|---:|---:|---:|
| R_+ rank, real | 2 | 2 | 0 | 0 |

Thus the active-to-fixed coupling rank is at most four real dimensions for
any complementary projectors P,Q. The selected local split attains rank
four numerically, with normalized nonzero singular values approximately
0.0677251 (twice) and 0.0370214 (twice). Other singular values are below
7e-16. The repeated pairs reflect the compatible complex structure; the
observed coupling rank is two over the complex numbers.

Although R_+ is supported in ell plus old, multiplication by P and Q can
spread the endpoints of Q R_+ P through other arithmetic coordinates.
The factorization is a rank statement and does not confine every coupling
endpoint to ell plus old.

## Equal signs remove a commutator component

Suppose arithmetic projectors E_b,E_c commute with W, and W acts on them
by scalars s_b,s_c. For any matrix P,

\[
E_b[P,W]E_c=(s_c-s_b)E_bPE_c.
\]

For W143, both q4 and q6 have sign -1. Their two cross-commutator blocks
therefore vanish exactly for every P. In contrast, W11 and W13 have
opposite signs on q4 and q6, so their corresponding squared commutator
blocks are four times the squared projector blocks.

The sign identities and rank-four upper bound are EXACT for the packaged
matrices. The selected projector, its coupling rank saturation, and the
measured magnitudes are DIAGNOSTIC with numerical Hodge data. These signs
are eigenvalues of an involution, not a spacetime signature.
