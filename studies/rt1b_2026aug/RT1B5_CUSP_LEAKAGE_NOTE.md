# RT-1B.5 — the bad-prime boundary maps, exactly

Filed 2026-08-18. Script: `rt1b5_cusp_leakage_matrices.gp`.

## Sol's identity, verified — and sharpened

For weight 4, a_1(T_l G) = a_l(G), since the second Hecke term cannot
reach n = 1. Verified as `v . T_l == (a_l row)` for l = 2, 3, 5, 7.

**But it also holds for U_11 and U_13.** a_1(U_p F) = a_p(F) just as
directly. So the identity itself is not where good and bad primes differ.
The difference is entirely that **T_l commutes with W_d and U_p does not**:

    c_d(T_l F) = a_l(F | W_d)   for all four d          <- good primes
    c_d(U_p F) != a_p(F | W_d)  in general              <- level primes

That is the precise locus of the phenomenon, and it is a small
strengthening of the derivation rather than a correction to it.

## The images, in cusp coordinates

Cusps indexed by d | 143 via c_d(F) = a_1(F | W_d). Computing the image
of the leakage map on H^0(2K) and converting from the AL character basis:

| operator | rank | image | exact relations on the image |
|---|---|---|---|
| T_2, T_3, T_5, T_7 | 4 | all of C^4 | none |
| **U_11** | 2 | span{ e_1, e_13 } | **c_11 = c_143 = 0** |
| **U_13** | 1 | span{ (1,1,0,0) } | **c_13 = c_143 = 0 and c_1 = c_11** |

So the rank statistics become equations.

**U_p annihilates exactly the cusp channels whose index is divisible
by p.** U_11 kills the cusps indexed 11 and 143; U_13 kills 13 and 143.
That is the local statement one would want: the level prime degenerates
precisely at the cusps where the level structure is p-ramified.

U_13 satisfies one relation beyond that — c_1 = c_11 — which is why its
rank is 1 rather than 2. U_11 satisfies no extra relation. The asymmetry
between the two level primes is now an explicit equation rather than the
gap between two integers.

## Consistency with the earlier partial grading

Two rounds ago the AL partial grading gave im l_11 splitting 1 + 1 across
the two W_13-halves, and im l_13 lying entirely in the W_11 = +1 half.
Both fall out of the equations above:

* span{e_1, e_13} is W_13-stable (W_13 swaps 1 and 13) and splits into
  e_1 + e_13 and e_1 - e_13 — one line in each parity. Matches (1,1).
* (1,1,0,0) is fixed by W_11, which swaps 1 with 11 and 13 with 143.
  Matches "entirely in the W_11 = +1 half".

Two independent computations, one from rank arithmetic on AL-graded
subspaces and one from explicit cusp functionals, agreeing.

## Split of Q1, adopted

Sol's separation into Boundary Saturation A (cusps) and B (elliptic) is
right, and A now has a concrete formulation: for good l, cusp saturation
is exactly the linear independence of the four functionals
F -> a_l(F | W_d) on H^0(2K), equivalently the non-vanishing of all four
AL-character components. B still needs elliptic local functionals, which
I have not built; that is the harder half and the honest blocker.

## What is next

`rt1b4_boundary_census.gp` now has all four gates hardened as errors;
verified running clean. The disjoint-range pre-registration is filed as
`PREREG_census_disjoint_range.md`, including Sol's strengthened
U_q-invariance form of the N = 3q candidate. Nothing in the 160 < N <= 400
range has been computed.
