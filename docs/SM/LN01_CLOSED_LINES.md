# LN-01 — the closed lines of the curve's unit tangent bundle: thirty integers, crossings and linking numbers (2026-11-16)

Modules `mtft.lines.darts`, `rademacher`, `closed`, `periods`, `involutions`, `cutting`, `crossings`, `lorenz`, `table`; tests
`tests/test_lines_darts.py`, `tests/test_lines_closed.py`, `tests/test_lines_periods.py`, `tests/test_lines_crossings.py`,
`tests/test_ln01_closed_lines.py`. Ported from the Modular Plane lab's records MP12, MP13 and MP14 (12 November 2026), whose scripts and
number files are kept under `studies/lines_2026nov/`. The frozen record `src/mtft/lines/_data/lines_143.json` (SHA-1 in `PROVENANCE.txt`)
holds the lab's values; the package computes every one of them again with its own code and the gates compare.

What this is: mathematics of M = Γ₀(143)\PSL(2, ℝ), the unit tangent bundle of X₀(143), whose geodesic flow is the frozen frame of the
lab's instrument 12.7. What it is not: a statement about gravity. Whether any of the integers below is a charge of a frozen field is
not known.

## 1. The class of a closed line (LN-1, EXACT)

Γ₀(143) has 168 cosets in SL(2, ℤ), the points of P¹(ℤ/143). With S = [[0, −1], [1, 0]], T = [[1, 1], [0, 1]] and R = S⁻¹T they carry a
graph: a vertex is an orbit of R, an edge an orbit of S, a face an orbit of T. It has 56 vertices, 84 edges (one a loop) and 4 faces of
lengths 1, 11, 13, 143, the widths of the cusps; b₁ = 29 and the genus is 13. The graph is a retract of the punctured curve, so
H₁(Γ₀(143)\ℍ, ℤ) = ℤ²⁹, and each edge outside a spanning tree gives one generator of the group: 29 matrices, the same as the lab's, in
the same order.

An element of the group is a word in S and T; read from left to right it moves a dart, and its 29 signed crossings of the edges outside
the tree are its class, `Graph.cls`. The class is additive and unchanged by conjugation (200 random pairs). M is a circle bundle over the
punctured curve, so H₁(M, ℤ) = ℤ³⁰, and a closed line has the class (Ψ; cls): Rademacher's Ψ along the fibre, as at level 1, and cls.
They are bound by one congruence, Ψ ≡ χ(cls) modulo 6 (χ(T) = 1, χ(S) = 3), and the classes of closed lines fill a lattice of index 6
in ℤ³⁰ (60 random hyperbolic elements: rank 30, Smith invariants with product 6).

## 2. The 28 lines (LN-2, EXACT)

A hyperbolic matrix A of SL(2, ℤ) permutes the 168 darts, and each cycle of length ℓ through the coset of g gives one closed line,
γ = g A^ℓ g⁻¹, of length ℓ · 2 arccosh(|tr A|/2). Three motions: the cat map [[2, 1], [1, 1]] (cycles 7⁴ 35⁴, Ψ = 0), Sphereland 10's
motion [[84, −47], [143, −80]] (1⁴ 5⁴ 6⁴ 30⁴, Ψ = −1) and the Pell unit [[12, 1], [143, 12]] (1, 11, 13, 143, Ψ = −21).

| # | motion | ℓ | dart | Ψ | (φ₁, φ_A, φ_B) | 143a1 | oldspace | W₁₁ W₁₃ W₁₄₃ | mirror |
|---|---|---|---|---|---|---|---|---|---|
| 1 | cat | 7 | (1 : 3) | 0 | (0, 0, 0) | (0, −1) | (−9, 5, 2, −1) | · · · | + |
| 2 | cat | 7 | (1 : 7) | 0 | (0, 0, 0) | (0, −1) | (9, −5, −2, 1) | · · · | + |
| 3 | cat | 7 | (1 : 14) | 0 | (0, 0, 0) | (0, 1) | (−9, 5, 2, −1) | · · · | + |
| 4 | cat | 7 | (1 : 29) | 0 | (0, 0, 0) | (0, 1) | (9, −5, −2, 1) | · · · | + |
| 5 | cat | 35 | (0 : 1) | 0 | (5, −84, 0) | (2, −1) | (27, −13, 20, −5) | · · · | − |
| 6 | cat | 35 | (1 : 0) | 0 | (−5, 84, 0) | (−2, 1) | (−27, 13, −20, 5) | · · · | − |
| 7 | cat | 35 | (1 : 2) | 0 | (−5, −84, 0) | (−2, 1) | (27, −13, 20, −5) | · · · | − |
| 8 | cat | 35 | (1 : 4) | 0 | (5, 84, 0) | (2, −1) | (−27, 13, −20, 5) | · · · | − |
| 9 | s10 | 1 | (0 : 1) | −1 | (0, −1, 0) | (0, 0) | (−1, 0, 0, 0) | · · − | · |
| 10 | s10 | 1 | (1 : 25) | −1 | (0, −1, 0) | (0, 0) | (−3, 1, 0, 0) | · · − | · |
| 11 | s10 | 1 | (11 : 2) | −1 | (0, −1, 0) | (0, 0) | (−1, 0, 0, 0) | · · − | · |
| 12 | s10 | 1 | (13 : 6) | −1 | (0, −1, 0) | (0, 0) | (−3, 1, 0, 0) | · · − | · |
| 13 | s10 | 5 | (1 : 12) | −5 | (1, 7, 0) | (0, 1) | (−7, 3, −4, 1) | · · · | · |
| 14 | s10 | 5 | (1 : 38) | −5 | (−1, −5, 0) | (0, −1) | (11, −4, 4, −1) | · · · | · |
| 15 | s10 | 5 | (13 : 1) | −5 | (1, 7, 0) | (0, −1) | (−7, 3, −4, 1) | · · · | · |
| 16 | s10 | 5 | (13 : 4) | −5 | (−1, −5, 0) | (0, 1) | (11, −4, 4, −1) | · · · | · |
| 17 | s10 | 6 | (1 : 3) | −6 | (0, −6, −2) | (−1, 2) | (12, −5, 11, −3) | · · · | · |
| 18 | s10 | 6 | (1 : 14) | −6 | (0, −6, 2) | (1, −2) | (12, −5, 11, −3) | · · · | · |
| 19 | s10 | 6 | (11 : 1) | −6 | (0, −6, 2) | (1, 1) | (19, −9, 9, −2) | · · · | · |
| 20 | s10 | 6 | (11 : 4) | −6 | (0, −6, −2) | (−1, −1) | (19, −9, 9, −2) | · · · | · |
| 21 | s10 | 30 | (1 : 0) | −30 | (−9, −30, −10) | (2, −1) | (−11, 4, −4, 1) | · · · | − |
| 22 | s10 | 30 | (1 : 1) | −30 | (−11, 42, 10) | (0, 0) | (−20, 10, −16, 4) | · · · | − |
| 23 | s10 | 30 | (1 : 5) | −30 | (11, −30, 10) | (−2, 1) | (−11, 4, −4, 1) | · · · | − |
| 24 | s10 | 30 | (1 : 8) | −30 | (9, 42, −10) | (0, 0) | (−20, 10, −16, 4) | · · · | − |
| 25 | pell | 1 | (0 : 1) | −21 | (0, −15, −1) | (0, 0) | (0, 0, 0, 0) | · · − | − |
| 26 | pell | 11 | (13 : 1) | −231 | (0, 15, −11) | (0, 0) | (0, 0, 0, 0) | · · − | − |
| 27 | pell | 13 | (11 : 1) | −273 | (0, −195, 1) | (0, 0) | (0, 0, 0, 0) | · · − | − |
| 28 | pell | 143 | (1 : 0) | −3003 | (0, 195, 11) | (0, 0) | (0, 0, 0, 0) | · · − | − |

The last two columns give the sign of each involution on the line's 26 periods and of the mirror on its class (·: neither kept nor
negated). The coordinates on f2 (8) and f3 (12) are returned by `mtft.lines.periods.coordinates`. Gated: Ψ of a lift is ℓΨ(A); the
classes of all the lifts of a motion add up to zero, with their windings and their periods block by block; their fibre numbers add up to
168Ψ(A) (0, −168, −3528); each cat lift and its partner through S have opposite classes; no line is null-homologous; the 28 classes have
rank 16 (17 with Ψ; by motion 4, 13, 3). Every integer equals the lab's, and so do the 28 matrices, up to 199 digits long.

Sphereland 10's motion has trace 4, so its four lines of length 1 are systoles of X₀(143) (length 2 arccosh 2 = 2.6339; ledger M05-3).
Matter-05's own example, [[19, −2], [143, −15]], has the class and fibre number of line 12.

## 3. Thirty integers with a meaning (LN-3 EXACT, LN-4 NUMERICAL)

**Three windings.** A unit ∏ η(δτ)^{r_δ} has no zeros or poles on Γ₀(143)\ℍ, and the number of times it winds round 0 along γ is
(1/24) Σ r_δ Φ(γ_δ), γ_δ = [[a, δb], [c/δ, d]], Φ Rademacher's function (Dedekind sums). The three units are
u₁ = η(τ)η(143τ)/(η(11τ)η(13τ)), u_A = (η(τ)/η(11τ))¹² and u_B = (η(τ)/η(13τ))². Their windings are integers and integer linear forms
of cls (300 random elements). Two further routes: on the loop round each cusp they equal the unit's order there by Ligozat's formula
(widths 1, 143, 13, 11: φ₁ = 5, 5, −5, −5; φ_A = −5, 65, −65, 5; φ_B = −1, 11, 1, −11), and on six small elements they equal the
integral of the logarithmic derivative, done numerically (1e-8).

**Twenty-six periods.** H₁ of the closed curve carries the Hecke algebra with four blocks, the kernels of the factors of T₂'s
characteristic polynomial x²(x + 2)⁴g₄²h₆²: 143a1 (2), the level-11 oldspace (4), f2 (8), f3 (12). Two routes:

1. the lab's integer coordinates of the 29 generators on each block (frozen), applied to cls;
2. the modular symbol {∞, γ∞} in Manin symbols by the continued fraction of γ∞, resolved in the package's own model `mtft.hecke`.

The bases differ, so the gate compares what no basis changes: on the 28 lines the same components vanish, the ranks agree (2, 4, 5, 8 by
block, 14 in all) and so do the linear relations, block by block and whole. On the 29 generators the two are related by an invertible
matrix in every block, and the lab's coordinates span all of ℤ^dim (Smith invariants all 1): one lattice. H₁ of the closed curve has index
1327104 = 2¹⁴·3⁴ in the product of the four block lattices.

**The lattice of 143a1** (CERTIFIED, 1e-8). Integrating 2πi f dz for the newform of y² + y = x³ − x² − x − 2 (coefficients from point
counts, 6000 terms) along twelve small elements gives Ω₁ = −0.984962 + 2.015772i and Ω₂ = −1.969923, with residual 6 × 10⁻¹⁵ and the
curve's j, −141.0134. Second route: the 26 periods of that form stored in `mtft.periods` (from PARI/GP) lie in ℤΩ₁ + ℤΩ₂ to 2 × 10⁻¹⁵
and generate it. The Pell unit's period vanishes. Re Ω₁ = Ω₂/2, so the mirror is complex conjugation of the period.

## 4. The involutions, the fibre number and a miss (LN-5, LN-6, EXACT)

W₁₁ = [[66, 5], [143, 11]], W₁₃ = [[78, 7], [143, 13]] and W₁₄₃ = [[0, −1], [143, 0]] normalise the group. On ℤ²⁹ they are commuting
involutions with W₁₁W₁₃ = W₁₄₃; on the blocks they are the signs 143a1 (+, +, +), f2 (−, +, −), f3 (+, −, −), and W₁₁ = −1 on the
oldspace. The mirror τ ↦ −τ̄ reverses Ψ and the three windings and acts on the 143a1 block by (n₁, n₂) ↦ (−n₁, n₁ + n₂).

**The fibre number is not kept** (found during this build). A test asserted Ψ(W_Q γ W_Q⁻¹) = Ψ(γ) and failed on its first run. The
reason: W_Q = g · diag(Q, 1) with g in SL(2, ℤ), so the image is conjugate in SL(2, ℤ) to γ_Q = [[a, Qb], [c/Q, d]], and

    Ψ(W_Q γ W_Q⁻¹) = Ψ(γ) + Φ(γ_Q) − Φ(γ) = Ψ(γ) + the winding of (η(Qτ)/η(τ))²⁴ along γ:

W₁₁ gives Ψ − 2φ_A, W₁₃ gives Ψ − 12φ_B, W₁₄₃ gives Ψ + 24φ₁ − 2φ_A − 12φ_B. Derived after the failure, then checked exactly on 328
elements for each involution; the fibre number moved for 310, 228 and 322 of them. For the Pell unit, which W₁₄₃ carries to its inverse,
−21 + 30 + 12 = +21. So on the thirty integers an involution is (Ψ; cls) ↦ (Ψ + shift(cls); cls · W), with a shift that is odd under it.

**The miss.** Before its run the lab conjectured, with low confidence (MP12, P13), that every line with no 143a1 component is negated in
H₁ of the closed curve by one of the three involutions. It is false. Ten lines have none (9 to 12, 22, 24, 25 to 28). W₁₄₃ negates eight:
the four systoles and the four Pell lines. No involution negates lines 22 and 24. The converse is a theorem and holds. The two zeros have
no mechanism yet; they are flagged under AG-D5 and recorded as negative N23 of the puzzle map.

## 5. Every crossing, exactly (LN-7, EXACT)

Each motion's closed geodesic on the modular surface is cut along the fundamental domain F into arcs of semicircles
a(x² + y²) + bx + c = 0 with integer a, b, c: 4 arcs a period for the cat map (twice through the order-2 point i), 3 for Sphereland 10's
motion (once through ρ), 25 for the Pell unit. At level 143 the arc of frame m of the line through the dart x lies in the sheet x·m, and two
arcs in one sheet cross at a point with rational x and y². A crossing on the boundary of a sheet lies in two sheets, or six at a corner;
one half-open convention counts it once, and the mirrored convention gives the same matrices.

The 28 lines trace 24 geodesics (each cat lift and its partner run one geodesic in opposite directions). Between different geodesics
there are 3214 crossings, and 2918 self-crossings; 103 of the 276 pairs are disjoint and 15 geodesics are simple, the four
systoles among them, which are also pairwise disjoint. The crossings sit over i (84), over ρ (168), on the left edge (2184) and
inside (3696). The Pell line of length 143 crosses itself 2794 times. The matrices equal the lab's entry for entry (MP13), and with the three
motions UV³, U²V² and UV⁵ added, MP14's 38 × 38.

The signed counts (+1 when δ crosses γ from its right to its left) are a form on the 26 periods: every one of the 10 linear relations
among the period vectors of the 24 geodesics annihilates the signed matrix, which has rank 12. Its four Pell rows are zero, although the
Pell line of length 143 crosses every other geodesic. The mechanism is MP14's and is gated: the mirror reverses orientation, so two
geodesics that it reverses have signed count zero (10 of the 24, the Pell lines among them), and the Pell line of length 1 never meets
the other 14. The zero belongs to this set of 24: UV⁵'s four short lines cross the Pell line of length 1 once each, with signs
+1, +1, −1, −1.

## 6. Linking numbers, and a second route to their size (LN-8, STRUCTURAL, caveat)

The classes of all the lifts of a motion add up to zero but their fibre numbers add up to 168Ψ(A), so the union bounds only for the cat
map. What bounds for every line is the line with its reverse, K_γ. The lab derived, before its run (MP13), from the half-turn annulus
over γ in the orientation dx ∧ dy ∧ dθ,

    Lk(K_γ, K_δ) = −I(γ, δ) ≤ 0,

minus the number of crossings on the curve. For the four cat unions:

| | K(1, 4) | K(2, 3) | K(5, 6) | K(7, 8) |
|---|---|---|---|---|
| K(1, 4), length 7 | — | 0 | −6 | −8 |
| K(2, 3), length 7 | 0 | — | −8 | −6 |
| K(5, 6), length 35 | −6 | −8 | — | −24 |
| K(7, 8), length 35 | −8 | −6 | −24 | — |

This block has a route of its own with no cutting (the lab's Route B): the cat geodesics cross only over the order-2 point, once at
each of its 84 lifts, and at the lift {y, yS} the two branches are the cat geodesics through the darts y and yT⁻¹. Counting those 84 pairs
gives the same 4 × 4 block, self-crossings included.

**The second route** shares no step with the cutting. Modular knots are Lorenz knots (Ghys): a hyperbolic matrix is conjugate to a
positive word in x = U and y = V, and two periodic orbits of the template link half as many times as their strands cross where the left
strip passes the right. The words: cat: xy; S10: xxy and xyy; Pell: x²²y and xy²²; UV³: xxxy and xyyy; U²V²: xxyy; UV⁵: xxxxxy and xyyyyy. M covers level 1 with 168 sheets, so the
total of crossings between the geodesics over two motions must be 168 times the sum of the template's linking numbers over the distinct
words of the two (a motion conjugate to its inverse, like the cat map and U²V², has one word). The lab checked three pairs. Here all
fifteen pairs of the six motions are checked, and out of sample at levels 11 and 23, where the factor is 12 and 24:

| pair | template sum | × 168 | cutting, level 143 | level 11: × 12 / cutting | level 23: × 24 / cutting |
|---|---|---|---|---|---|
| cat × S10 | 2 | 336 | 336 | 24 / 24 | 48 / 48 |
| cat × Pell | 2 | 336 | 336 | 24 / 24 | 48 / 48 |
| cat × UV³ | 2 | 336 | 336 | 24 / 24 | 48 / 48 |
| cat × U²V² | 2 | 336 | 336 | 24 / 24 | 48 / 48 |
| cat × UV⁵ | 2 | 336 | 336 | 24 / 24 | 48 / 48 |
| S10 × Pell | 6 | 1008 | 1008 | 72 / 72 | 144 / 144 |
| S10 × UV³ | 6 | 1008 | 1008 | 72 / 72 | 144 / 144 |
| S10 × U²V² | 4 | 672 | 672 | 48 / 48 | 96 / 96 |
| S10 × UV⁵ | 6 | 1008 | 1008 | 72 / 72 | 144 / 144 |
| Pell × UV³ | 8 | 1344 | 1344 | 96 / 96 | 192 / 192 |
| Pell × U²V² | 4 | 672 | 672 | 48 / 48 | 96 / 96 |
| Pell × UV⁵ | 12 | 2016 | 2016 | 144 / 144 | 288 / 288 |
| UV³ × U²V² | 4 | 672 | 672 | 48 / 48 | 96 / 96 |
| UV³ × UV⁵ | 8 | 1344 | 1344 | 96 / 96 | 192 / 192 |
| U²V² × UV⁵ | 4 | 672 | 672 | 48 / 48 | 96 / 96 |

The caveat: the sign of the linking formula assumes that the orientation used here is the one of Ghys's left-handedness. That was not
checked; the size is what the two routes confirm.

## 7. Forecasts and what happened to them

Fifteen forecasts were written with the clock before the code they concern was run; the full list with times is in the build record
(`MTFT_v110_Build_Record.md`, wave 2). The gates forecast for the port held on their first runs: the graph and its generators, the 28
lines, the periods by both routes, the Atkin–Lehner group law, the 24 × 24 and 38 × 38 crossing matrices with signs, the transfer at
levels 11 and 23, and the Lorenz route on all fifteen pairs, the cat × U²V² pair included, which had been given 0.4. The forecast that at
least one gate would fail through a fault of the port or the test (0.75) also held, and every fault was in a test, none in the port.
Five tests failed on their first runs, each through something its author had written without computing it: a loop count at level 11
(0 written, 1 true at every level), the block ranks of the 28 lines, one entry of the crossing matrix, a "new" motion that was the cat
map squared (the package refused it, as it should), and the claim that the involutions keep Ψ. The first three were corrected to the
computed values and say so; the fourth was replaced by a primitive motion; the fifth was a wrong belief, and what replaced it is in
section 4.

## 8. What is not here

- A mechanism for the two zeros of section 4.
- The lab's formula for a union L with zero periods, Lk(K_γ, L) = −½N(γ, L), is stated in `mtft.lines.crossings` and used once (the cat
  lifts against all the lifts of another motion: −168); it has no gate of its own beyond that.
- Levels with elliptic points. The graph and the counting are refused there (13, for instance), so MP14's shadows on X₀(13) are a recorded
  result of the lab's scripts and not reproduced.
- Any physical reading. The helicity of 12.7's static frame has not been computed against these classes.
