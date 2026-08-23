# The larger-prime dependency is an equality, not a proportionality

Filed 2026-08-18. Discovery levels only (N <= 160). Script:
`rt1b_dependency_constant.gp`. **The holdout was not touched.**

Sol's second theorem target was

    det [ c_1 . U_q ; c_p . U_q ] = 0

i.e. the two surviving cusp functionals of the larger level prime are
dependent. The constant turns out to be forced:

| N = pq | U_q surviving pair | U_p pair |
|---|---|---|
| 55 = 5·11 | c_5 . U_11 = **1** · c_1 . U_11 | c_11 . U_5 = **1** · c_1 . U_5 |
| 65, 85, 95, 115, 145, 155 (p = 5) | c_5 . U_q = **1** · c_1 . U_q | independent (rank 2) |
| 77 = 7·11 | c_7 . U_11 = **1** · c_1 . U_11 | c_11 . U_7 = **1** · c_1 . U_7 |
| 91, 119, 133 (p = 7) | c_7 . U_q = **1** · c_1 . U_q | independent (rank 2) |
| 143 = 11·13 | c_11 . U_13 = **1** · c_1 . U_13 | independent (rank 2) |

**Twelve for twelve, and the constant is 1 every time.** Where the smaller
prime's pair is also dependent (only 55 and 77), the constant is 1 there
as well.

So the sharper target is not a vanishing determinant but

    c_p . U_q  =  c_1 . U_q     on H^0(2K)

— the two functionals are **equal**, not merely proportional. Combined
with the vanishing support this gives a complete description of the larger
level prime's cusp leakage at every discovery level:

    c_q(U_q F) = c_pq(U_q F) = 0        and        c_1(U_q F) = c_p(U_q F)

which is exactly the N = 143 result (c_13 = c_143 = 0, c_1 = c_11)
holding across the whole family.

**Normalization, per Sol's caveat.** The value 1 is pinned to
`mfatkininit`'s normalization of W_d and must not be quoted as invariant.
What *is* invariant here is that the ratio is **level-independent** — had
it varied with N, that variation would have been arithmetic content. It
does not vary, so the whole effect is a single uniform relation rather
than a family of constants to explain.

The two theorem targets are now cleanly separated:

| target | controls |
|---|---|
| im l^cusp_{U_p} contained in span{e_d : p does not divide d} | which channels can exist |
| c_p . U_q = c_1 . U_q | why the surviving pair collapses to one dimension |
