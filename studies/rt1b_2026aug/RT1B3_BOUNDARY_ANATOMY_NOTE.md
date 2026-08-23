# RT-1B.3 — Boundary Anatomy, and the four artifact fixes

Filed 2026-08-18. Supersedes both scripts staged earlier today.

## Artifact fixes — all four conceded

**1. Obsolete header.** `rt1b_crosslevel.gp` printed "boundary dim is 4
for every squarefree N = p*q" while the note retracted it. Rewritten; the
header now states the true formula, b_N = nu_infinity + nu_2 + nu_3, and
flags the previous claim as superseded.

**2. The hyperelliptic gate is now a gate.** The old script computed
rank mu_2, printed it beside 3g-3, and proceeded. It now halts the level
and prints `INVALID_FOR_CANONICAL_LEAKAGE` before the leakage section.
Verified: N = 35 and 39 are skipped, N = 55 onward proceed.

**3. Provenance, one layer deeper.** The sweep read
`studies/Bmat.txt`, my local scratch file, absent from the v0.17.0 sdist.
It now parses `src/mtft/canonical/_data/X0_143_AL_adapted_basis.txt` —
which *is* in the release — with a `readB` parser in the script, and
prints det B = -1078272 as a self-check on load. Verified running from the
release directory.

**4. Count.** "Three-level, twenty-prime" was wrong; it double-counted the
cross-level reruns of (143, 2) and (143, 3). Correct at that time:
**18 distinct good (N, p) pairs, 13 distinct primes**. With the levels
added below it is now **31 distinct good (N, p) pairs across 7 valid
levels**.

Also caught while rewriting: **N = 39 has nu_3 = 2**, so its boundary is
6, not 4. Another instance of the same wrong assumption.

## The anatomy

B_N splits along the filtration H^0(2K) < S_4^(2) < S_4, where S_4^(2) is
the double-vanishing subspace. Since Atkin-Lehner acts freely transitively
on the four cusps of a squarefree N = pq, membership in S_4^(2) is
testable as a_1(f|W) = 0 over the four AL elements — no cusp-expansion
machinery. Then S_4/S_4^(2) is the cusp channel (dim nu_infinity = 4) and
S_4^(2)/H^0(2K) is the elliptic channel (dim nu_2 + nu_3).

Every operator gets a signature (r_cusp, r_ell). Computed with PARI's own
`mfheckemat` and `mfatkininit` — an independent implementation from the
q-expansion convolutions, and it **reproduces every earlier total exactly**.

| N | b_N | (nu_inf, nu_2, nu_3) | good T_p | U_p | U_q |
|---|---|---|---|---|---|
| 55 = 5·11 | 4 | (4, 0, 0) | (4, 0) | U_5 (1, 0) | U_11 (1, 0) |
| 65 = 5·13 | 8 | (4, **4**, 0) | (4, 4) | U_5 (2, **1**) | U_13 (1, 0) |
| 77 = 7·11 | 4 | (4, 0, 0) | (4, 0) | U_7 (1, 0) | U_11 (1, 0) |
| 85 = 5·17 | 8 | (4, **4**, 0) | (4, 4) | U_5 (2, **1**) | U_17 (1, 0) |
| 91 = 7·13 | 8 | (4, 0, **4**) | (4, 4) | U_7 (2, **1**) | U_13 (1, 0) |
| 95 = 5·19 | 4 | (4, 0, 0) | (4, 0) | U_5 (2, 0) | U_19 (1, 0) |
| 143 = 11·13 | 4 | (4, 0, 0) | (4, 0) | U_11 (2, 0) | U_13 (1, 0) |

### What this establishes

**Good primes saturate both channels simultaneously**, at all seven valid
levels: r_cusp = nu_infinity and r_ell = nu_2 + nu_3, always. Sol's
conjecture that "good primes saturate the boundary" is in fact **two**
saturation statements, and both hold everywhere tested.

**U_7 at N = 91 is (2, 1)** — neither (3,0,0) nor (1,0,2). It splits
across both channel types, taking 2 of 4 cusp channels and 1 of 4
elliptic.

**At every level with elliptic channels, exactly one of the two level
primes reaches them, and it takes exactly one.** The other takes none.

**The two level primes are never symmetric.** One always gives r_cusp = 1;
the other gives 1 or 2. That asymmetry is the finding.

### What this does NOT establish

A rule for r_cusp of the level primes. The values are 1 or 2 and I tried
to explain the split by which prime is larger, by the presence of
oldspace, and by the genus of X_0(q); **each attempt failed on at least
one of the seven levels**. Seven levels and fourteen level-primes are not
enough, and I am not proposing a formula. Registered as the open item.

## Status

| item | status |
|---|---|
| b_N = nu_inf + nu_2 + nu_3 | verified, 7 levels |
| good primes saturate both channels | computation, 31 (N,p) pairs, 7 levels |
| exactly one level prime reaches the elliptic channel, taking 1 | computation, 3 levels |
| r_cusp rule for level primes | **open, no rule claimed** |
| cross-check of totals by an independent operator implementation | passed |
