# I_2(X0(143)) in the Atkin-Lehner adapted basis

**Session 2026-08-16, continuation. Verified in adapted coordinates:
Route C worst relative residual 4.52e-41 at dps 40.**

---

## 1. The adapted basis

Every Galois block of S_2(Gamma_0(143)) is Atkin-Lehner pure (established
last session, structurally forced). So a single basis can carry BOTH
labels. Built as the simultaneous (block, sector) eigenlattice, one
`matkerint` per piece:

| coords | block | sector | dim |
|---|---|---|---|
| y1 | f1 = 143a1 | (+,+) | 1 |
| y2..y7 | f3 orbit | (+,-) | 6 |
| y8 | oldspace, W13 = +1 | (-,+) | 1 |
| y9..y12 | f2 orbit | (-,+) | 4 |
| y13 | oldspace, W13 = -1 | (-,-) | 1 |

W_11 and W_13 are **simultaneously diagonal**:
W11 = diag(+1 x7, -1 x6), W13 = diag(+1, -1 x6, +1 x5, -1). det B =
-1078272; the ideal was recomputed natively from the adapted
q-expansions rather than transformed, so no inverse and no denominators
enter. Recovered dim I_2 = 55 with **max exact residual 0** — the
adapted lattice and the transformed lattice agree, which is a free
consistency check on the change of basis.

## 2. What the grading buys

Because W acts diagonally, each quadric sits in one isotypic class and
its monomial support is confined to that class. Support drops from 91 to:

| class | monomials | dim (I_2)_psi | shape of a quadric in the class |
|---|---|---|---|
| (+,+) | 38 | 26 | a·y1^2 + Sym^2(f3)[21] + Sym^2(f2,old+)[15] + b·y13^2 |
| (+,-) | 11 | 5 | y1·L(f3)[6] + y13·M(f2,old+)[5] |
| (-,+) | 11 | 4 | y1·M(f2,old+)[5] + y13·L(f3)[6] |
| (-,-) | 31 | 20 | c·y1·y13 + f3 (x) f2 bilinear[24] + f3 (x) old+ [6] |

That is the readable payoff. The (+,+) class carries only **within-block**
quadratics. The two small classes, dim 5 and dim 4, are the ones in which
f1 and the (-,-) oldform act as connectors — every quadric there is
literally y1·(linear) + y13·(linear). The (-,-) class is where f2 and f3
meet.

**AG-D5 note.** The (-,-) class being 24-of-31 dominated by f2·f3
monomials is forced by dimension counting (4 x 6 = 24 versus 1 x 4 = 4
for f1·f2). It is arithmetic bookkeeping, not a dynamical statement about
generation mixing, and must not be read as one without a mechanism.

## 3. Sparsity: an honest negative

Beyond the symmetry-forced support, sparsity and coefficient size are in
genuine tension.

| presentation | max support per class | max coefficient |
|---|---|---|
| LLL-reduced (saturated) | 38 / 11 / 11 / 31 | 55 / 4684 / 10008 / 72 |
| HNF, best of 300 monomial orderings | 15 / 8 / 9 / 13 | 5e31 / 2e9 / 4e10 / 8e29 |
| echelon lower bound m - d + 1 | 13 / 7 / 8 / 12 | — |

LLL nearly saturates the class (dense, small entries); HNF nearly attains
the echelon support bound but at 30-digit coefficients. Suspecting the
adapted basis was badly scaled, I rebalanced it by `qflll` on the first 25
q-coefficients of each piece: max|E| fell 1664 -> 1296 and **the tradeoff
did not move** (coefficients 55/4684/10008/72 versus 55/4776/10080/104).
So this is intrinsic, not an artifact of basis choice. Filed as a negative
result: **the Atkin-Lehner grading is the whole of the available sparsity.**

The LLL presentation is the one shipped, since small coefficients are the
usable property.

## 4. The projection table — the actual finding

A global sign is invisible in projective space. So if every coordinate in
a subset S carries the same eigenvalue under some g in <W11, W13>, then
phi_S is constant on g-orbits and **factors through the quotient curve
X/H**, H the subgroup acting by a global scalar. The quadric count is then
pure Riemann-Roch on the quotient, with no reference to X itself.

Predictions were formed this way and then checked:

| S | H | g(X/H) | deg | h0(2H) | Sym^2 | predicted | computed |
|---|---|---|---|---|---|---|---|
| f3 | full group | 1 | 6 | 12 | 21 | **9** | **9** |
| f2 | full group | 1 | 6 | 12 | 10 | **0** | **0** |
| old | <W11> | 7 | 12 | 18 | 3 | **0** | **0** |
| f1+f3 | <W11> | 7 | 12 | 18 | 28 | **10** | **10** |
| f1+f2 | <W13> | 6 | 12 | 19 | 15 | 0 | 1 |
| f2+f3 | <W143> | 2 | 12 | 23 | 55 | **32** | **32** |
| newspace f1+f2+f3 | 1 | 13 | 24 | 36 | 66 | 30 | 33 |
| f2+f3+old | 1 | 13 | 24 | 36 | 78 | 42 | 44 |
| all 13 | 1 | 13 | 24 | 36 | 91 | **55** | **55** |

Six of nine exact. The headline case: **the six f3 coordinates are all
(+,-), so both W11 and W13 act on them by a global scalar; phi_f3
therefore factors through the full Atkin-Lehner quotient X0(143)*, which
has genus 1.** The image is an elliptic normal sextic in P^5 (degree
24/4 = 6), whose ideal is the classical n(n-3)/2 = 9 quadrics. The f3
sector sees the curve only through its genus-1 quotient. Same mechanism
gives f2 (degree-6 genus-1 curve in P^3, too few monomials for any
quadric), f2+f3 through the genus-2 quotient X/W143, and f1+f3 through
the genus-7 quotient X/W11.

**Three honest excesses**, where the restriction Sym^2(S) -> H^0(2K) fails
to have maximal rank:

- **newspace, excess 3.** The quadratic monomials in the 11 newform
  coordinates span only **33 of the 36 dimensions** of H^0(2K). The
  oldforms are not redundant: a 3-dimensional piece of H^0(2K) is
  reachable only with the level-11 ghost. This is the sharpest new fact
  in the table and it is worth its own study.
- **f2+f3+old, excess 2** — rank 34 of 36; dropping f1 alone costs 2.
- **f1+f2, excess 1** — the genus-6 image in P^4 lies on exactly one
  quadric where Riemann-Roch alone predicts none.

## 5. Status and next

Everything above is EXACT (integer/rational arithmetic, residual 0) except
Route C, which is CERTIFIED at 4.52e-41 relative, dps 40. No physics
reading is claimed; the pre-registration for this arc explicitly excludes
it.

Registered open items:
1. **The newspace deficiency of 3.** Which 3 dimensions of H^0(2K) require
   the oldspace, and are they Hecke-stable? Pre-register before computing.
2. The excess-1 quadric on the f1+f2 projection: identify it explicitly.
3. `mtft.canonical` — ship the adapted basis, the adapted q-expansions and
   the four class matrices, with the projection table as gates.

## Files

- `X0_143_AL_adapted_basis.txt` — B, with (block, sector) labels per coordinate.
- `X0_143_AL_adapted_qexpansions.txt` — e_1..e_13 to q^140, integral.
- `X0_143_I2_by_AL_sector.txt` — the ideal, split by class, LLL presentation.
- `X0_143_canonical_ideal_cert.json` — updated with the grading, the
  sparsity tradeoff and the projection table.
