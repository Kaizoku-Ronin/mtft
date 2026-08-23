# RT-1B round 4 — artifact fixed, sweep widened, and the pattern tested
# across levels (where it corrected me twice)

Filed 2026-08-18. Supersedes the artifact staged earlier the same day.

## 1. The provenance mismatch, fixed

Sol is right and this one is entirely mine. I staged `sweep.gp` — the
**uncorrected** version, which applies the T_p formula at p = 11 and 13 —
while the note reported the corrected U_p numbers computed in a separate
file I never staged. The artifact would have reproduced the discarded
calculation. That is precisely the failure the gate discipline exists to
catch, and I was the one who introduced it two messages after auditing
someone else for the same class of problem.

`rt1b_leakage_sweep.gp` now selects the operator by divisibility in one
place, carries the supersession notice in its header, and reproduces the
reported numbers:

    T_p(f)_n = a_{pn} + p^(k-1) a_{n/p}   p not dividing N
    U_p(f)_n = a_{pn}                     p dividing N

## 2. Good primes, widened to p <= 47

Thirteen good primes, all four AL channels, at q-depth 2700:

| p | (+,+) | (+,-) | (-,+) | (-,-) | total |
|---|---|---|---|---|---|
| 2, 3, 5, 7, 17, 19, 23, 29, 31, 37, 41, 43, 47 | 1 | 1 | 1 | 1 | **4** |

Thirteen primes, no drop. Keeping Sol's phrasing: this is computation,
not a theorem. RT-1B.2 stands open.

## 3. Bad primes, with the partial grading Sol proposed

U_p does not commute with W_p, so per-sector ranks are undefined. But
U_11 commutes with W_13 and U_13 with W_11, so a two-block grading
survives — and unlike the bogus sector sums from the wrong operator
(6 and 4 against ungraded 2 and 1), these are consistent:

| | ungraded | partial grading |
|---|---|---|
| **U_11** | **2** | W_13 = +1 block: **1**; W_13 = -1 block: **1** |
| **U_13** | **1** | W_11 = +1 block: **1**; W_11 = -1 block: **0** |

So im l_11 meets each W_13-half of B in a line, while im l_13 lies
**entirely inside the W_11 = +1 half** — the span of the cusp characters
chi_++ and chi_+-. That is a structural fact about which cusp channels
the level primes can reach, and it is the sharpest form of the "why
4 -> 2 -> 1" question.

## 4. Cross-level test — which corrected me twice

Three data points are not a pattern, so I ran the same computation at
other squarefree N = pq.

**Correction A. "The boundary is 4-dimensional for every squarefree
N = pq" is FALSE**, and I asserted it in the previous note and in a
script header. The codimension of H^0(2K) in S_4 is the number of cusps
**plus contributions from elliptic points**. N = 143 has no elliptic
points (11 = 2 mod 3 kills nu_3, 11 = 3 mod 4 kills nu_2), so there the
boundary is 4 = nu_infinity. **N = 91 has nu_3 = 4**, so its boundary is
**8**, and the good-prime leakage is 8, not 4.

**Correction B.** At N = 35 and 39 the product span came out **5**, not
3g - 3 = 6. The multiplication map is not surjective there, so the space
being tested is not H^0(2K) and the leakage numbers are not comparable.
Both levels are hyperelliptic — and both are on Ogg's list, while 55, 91
and 143 are not. The computation re-derives that classification for free,
and the check should be a gate on any future level added to the sweep.

Valid (non-hyperelliptic) levels:

| N | genus | boundary | good primes | U_p, p \| N |
|---|---|---|---|---|
| 55 = 5·11 | 5 | 4 | **4** | U_5 = 1, U_11 = 1 |
| 91 = 7·13 | 7 | **8** | **8** | U_7 = 3, U_13 = 1 |
| 143 = 11·13 | 13 | 4 | **4** | U_11 = 2, U_13 = 1 |

**The structural statement survives in a better form.** Good primes
saturate the *entire* boundary at every valid level tested, whatever its
dimension — that is now a three-level, twenty-prime statement rather than
a 143-specific observation. Bad primes never saturate it.

The bad-prime ranks are genuinely level-dependent: U_11 gives 1 at N = 55
but 2 at N = 143, so the value depends on the pair, not on the prime.
Three levels is not enough to claim a formula, and I am not proposing one.

## 5. Conceded: the finite-Aut gate as I stated it

Sol is right that the lattice argument does not work. The Petersson Gram
matrix on the q-expansion lattice need not be an integral form, and an
automorphism need not preserve that particular lattice, which is defined
using the cusp at infinity. Withdrawn.

Sol's replacement is the correct target: compute
#Stab_{PGL_13}(I_2) and ask whether it is 4 — the finite analogue of
RT-0, on the object already certified. Current status:

* **>= 4 is certified** — I_2 is W_11- and W_13-stable, verified in
  session 2 and gated in v0.17.0 (`gate_sector_grading` depends on it).
* **Lie algebra = 0 is certified** — RT-0.
* **<= 4 is open.** It is a polynomial system in 169 unknowns, not a
  linear one, and I do not think it falls to the tooling in the container.
  It wants a real algebraic-geometry engine. Flagging honestly rather
  than proposing a method I cannot execute.

## 6. Agreed

RT-1A closed. The "canonical-ring geometry vs Hecke-stable automorphic
geometry" framing is right and better than the ghost-M_3 line it replaced:
H^0(2K) is a distinguished subspace of the Hecke module S_4, of
codimension nu_infinity + (elliptic terms), and the boundary map measures
the mismatch. Pausing RT-2C for one arc on the boundary layer is the right
call. Terminology firewall unchanged.
