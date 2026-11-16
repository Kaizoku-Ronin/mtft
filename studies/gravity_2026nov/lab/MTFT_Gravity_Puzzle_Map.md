# The gravity programme — puzzle map

Written by hand on 10 November 2026, in the layout of `MTFT_Matter_Puzzle_Map.md`. The matter map is generated from a claims ledger; gravity has no ledger yet, so the last column points to records (notes, scripts, instruments, papers) instead of claim numbers, and the first job on this list is to give gravity a ledger of its own. Each piece is one question; it carries its own sub-puzzles. Statuses: stands, partial, open, closed (an honest negative); sub-puzzles are open, answered, failed or moot. A layer tag in brackets says what kind of statement an answer is: [structure] established mathematics or physics, reproduced; [derived] this programme's derivation, unreviewed; [model] an assumption with numbers attached; [hypothesis] not yet tested.

11 pieces, 58 sub-puzzles (27 open), 10 pre-registered tests (1 pending, T-G3; the other nine done), 16 recorded negatives. Last changed 14 November 2026, for T-G7 (MP21), T-G8 (MP22), and T-G1 and T-G2 (MP23): IV.9, III.4 and III.5 answered; two earlier statements about the gravitational hopfion corrected (N-G15); the square shown not to hold for potentials, and a hand value of the helicity ratio missed (N-G16). Before that, 12 November 2026, for T-G10 (MP14): II.10 answered; in piece X, Paper 33 §5.4's zero period explained and the cusps' scattering matrix noted as a lead. Earlier the same day: T-G9 (MP13), which answered II.8 with a correction to MP12 and opened II.10 and N-G14; T-G4 (MP12), which answered II.4 and opened II.8, II.9 and N-G13; and the note on Sphereland's two charts, with leads under I.6 and IV.9, additions to III.6 and IX.2, and N-G12. **Newton's constant is not derived**, and nothing below changes that.

## The pieces in dependency order

| Layer | Piece | Status | Open sub-puzzles | Depends on |
|---|---|---|---|---|
| 0 | I The coupling: Newton's constant | open | 3 of 6 | — |
| 0 | II The arena: where gravity sits on X₀(143) | open | 4 of 10 | — |
| 1 | III Gravity as a square (the coin) | partial | 1 of 6 | — |
| 1 | IV The frozen field | partial | 2 of 9 | III |
| 2 | V Breaking and reconnecting | open | 3 of 4 | IV |
| 2 | VI The field equations | open | 2 of 5 | I, II |
| 3 | VII The vacuum and the expansion | open | 4 of 5 | I, VI |
| 3 | VIII Dark matter | open | 2 of 2 | I, VI |
| 3 | IX Gravitational waves | open | 2 of 3 | III, IV |
| 3 | X Black holes | open | 2 of 3 | II, VI |
| 3 | XI Is gravity quantum? The low-mass frontier | open | 2 of 5 | I, III |

## I — The coupling: Newton's constant (open)

**Question.** Can G, or the Planck mass, be an output of the arithmetic, as α is claimed to be?

**Verdict.** Not yet. The compilation lists it as the "third anchor" problem (Part V, target A) and the Arithmetica Generale says outright that the five-primitive pipeline does not produce G (§13.5 iii). Paper 31's formula G = (ħc/mₑ²)|η(iy_c)|^{2k} is a model: it gives +3.7 % at y_c = M ln 2 and +61 % at the lab's depth 0.18174, with d ln G/dy = 909, so it is a seventh-decimal statement about a depth the programme has not fixed (R&D note §27 iv). What is exact is dimension counting and the 2 + 1 theory, where G is the level of a Chern–Simons theory.

**Records.** MTFT_Complete_Compilation.docx (Part V); Arithmetica Generale §13.5; RD_Notes §27; `gravity_knot_20261021.py`; Sphereland 18.

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| I.1 A derivation of G from the arithmetic | open | A route from the weights or the curve to M_P that fixes its own inputs; the compilation's Fredholm-determinant route (Paper 3) is not quantitative. | Compilation V.A |
| I.2 Is Paper 31's formula a derivation? | answered [model] | No: a model. +3.7 % at y_c = M ln 2 = 0.181256, +61 % at 0.18174; sensitivity 909. | RD_Notes §27 iv |
| I.3 Which depth does the programme mean, to seven digits? | open | An audit of y_c: the two depths in use differ by a factor 1.6 in G. | RD_Notes §27 |
| I.4 The bridge mₑ/M_P = η(iy_c)^k (Paper 22; η is real and positive on the imaginary axis) | open [hypothesis] | A test of the bridge that does not use G itself: it is the whole content of the formula. | RD_Notes §27 iv |
| I.5 What does the dimension of G say? | answered [structure] | [G] = M⁻¹L^{d−1}T⁻² is Kepler's cusp r^{d−1} = cT², whose link is T(2, d − 1): the trefoil in 3 + 1. G rescales the cusp without changing the knot. | RD_Notes §27 i; Sphereland 18 |
| I.6 Where is G literally a knot parameter? | answered [structure] | In 2 + 1 dimensions: gravity with Λ < 0 is Chern–Simons theory at level k = ℓ/4G, and Wilson loops are Jones polynomials at q = e^{2πi/(k+2)}. In 3 + 1 there is no level in general. A lead (12 November 2026, from the sphere charts; the literature cited from memory): the theory depends on the framing of the 3-manifold through e^{2πic/24} with c = 3ℓ/2G (Witten 1989; Brown and Henneaux 1986), and Witten's 2007 conjecture c = 24k would make ℓ = 16Gk, the 24 being that of the third stable stem and of η. A conjecture with known difficulties, not a result. | RD_Notes §27 iii; `sphere_charts_20261112.md` §4 |

## II — The arena: where gravity sits on X₀(143) (open)

**Question.** The curve is used for gauge forces and masses. What object on it, or over it, carries gravity?

**Verdict.** Open (compilation V.F). Three candidates are on record and none is established: the boundary (the early identification of gravity with the Feigenbaum point of the Mandelbrot set, and Paper 33's cusps as horizons); the 2 + 1 reading in which the modular knots are BTZ horizons (§27 iii, exact in 2 + 1); and the 3-fold Γ₀(143)\PSL(2, ℝ), on which a static frame's frozen lines are the closed geodesics (Room 12.7). The third is a model space with exact counts, not a claim about our universe. On it the Rademacher number of a closed line becomes a class in ℤ³⁰ (T-G4, 12 November): thirty integers, sorted by level and by Hecke block, and still no statement about gravity. A line never bounds alone; with its reverse it always does, and two such pairs link by minus the crossings of their shadows (T-G9). The Pell lines' zero signed counts come from the mirror, from shadows that miss each other one level down, and from homology (T-G10).

**Records.** Compilation V.F; Paper 33; RD_Notes §26, §27; `frozen_gravity_20261110.md` (E); The Modular Plane 12.7; `mp8_verify.py` (G1 to G8); MP12, `mp12_tg4_verify.py`, `mp12_tg4_mtft.py`, `test_core_tg4.js` (T-G4); MP13, `mp13_tg9_verify.py`, `test_core_tg9.js`, `mp13_lorenz.py` (T-G9); MP14, `mp14_tg10_verify.py`, `mp14_lattice.py` (T-G10).

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| II.1 "Gravity lives at the boundary" made precise on X₀(143) | open [hypothesis] | A map from the Mandelbrot boundary statement to the cusps or the boundary of the fundamental domain, with something computed on both sides. | Compilation V.F |
| II.2 Is there a natural frozen field over the curve? | answered [structure, model space] | Yes: on ℝ × Γ₀(143)\PSL(2, ℝ) the static frame is ideal; its legs' field lines are the geodesic flow, the fibre circles and the horocycles (hyperbolic, elliptic, parabolic). | 12.7; note E |
| II.3 How many closed lines lie over a given closed geodesic at level 143? | answered [structure] | One for each cycle of the motion on the 168 points of P¹(ℤ/143). Cat map: four of length 7L and four of 35L (predicted by hand, confirmed). Pell unit 12 + √143: 1, 11, 13, 143, the cycle type of T, because 12² ≡ 1 mod 143. Sphereland 10's systole: four each of 1, 5, 6, 30. | 12.7; `mp8_verify.py` G4 to G7 |
| II.4 What replaces the Rademacher linking number at level 143? | answered [structure] (T-G4, 12 November 2026) | At level 1 the 3-fold is the trefoil's complement, its first homology is ℤ, and Ψ(M) is the linking number (Ghys). At level 143 the 3-fold is a product (a surface of genus 13 with four punctures, times a circle) with first homology ℤ³⁰, and what replaces Ψ is the class of the closed line there: the fibre number, which is still Rademacher's, and a vector in ℤ²⁹, bound by one congruence, Ψ ≡ χ modulo 6. In coordinates with a meaning: Ψ; three windings round the cusps, read from the units η(τ)η(143τ)/(η(11τ)η(13τ)), (η(τ)/η(11τ))¹² and (η(τ)/η(13τ))²; and 26 periods, 2 + 4 + 8 + 12 along the Hecke blocks. Computed for 12.7's 28 closed lines: none has class zero, so no two of them have a linking number in Ghys's sense. Sorted by level [derived, post hoc]: Ψ comes from level 1, one winding and the four ghost periods from X₀(11), one winding from X₀(13), and one winding with 22 periods from level 143 alone. Three kinds of zero are symmetries: W₁₄₃ blinds the systoles and the Pell lines to 143a1, the mirror removes the windings of the cat lines of length 7, and W₁₁ their f2 part. This is mathematics of the model space; no statement about gravity follows from it yet. | MP12; `mp12_tg4_verify.py` (25 of 26 pre-registered checks, 6 post hoc), `mp12_tg4_mtft.py` (8), `test_core_tg4.js` (5); the lead was `sphere_charts_20261112.md` §4 (P17) |
| II.5 The BTZ entropy of the Pell knot as a Γ₀(143) statement | open | S = ℓ · 6.3526/4G against VE-01's shortest-geodesic term. | RD_Notes §27 |
| II.6 Is the volume 112π² a new number? | answered | No: 2π × area = 2π · (π/3) · 168. The index again (AG-D5). | 12.7; G8 |
| II.7 The Gödel reading of the fibration | open [hypothesis] | If the modular-time plane is the base of a Gödel fibration, the clock's return is a closed timelike curve: to be tested against causality, not adopted. | RD_Notes §26 |
| II.8 The linking numbers of unions that bound | answered [derived; the counts EXACT] (T-G9, 12 November 2026) | The classes cancel for all the lifts of one motion together, and for a cat lift with its partner through S (MP12, P5 and P6). Wanted: a linking number for such unions on the 3-fold, and what it depends on; on a product with second Betti number 29 the surface a union bounds is not unique. Correction (12 November, MP13): the lifts of a motion together bound only for the cat map; for the others they are homologous to 28Ψ(A) fibres. Answer: a closed line with its reverse always bounds (the cat's partner through S is its reverse), and Lk(K_γ, K_δ) = −I(γ, δ), minus the number of crossings of the two geodesics on the curve; Lk(K, L) needs only the periods of L to vanish. The four cat pairs: 0, −6, −8, −8, −6, −24, all crossings over the 84 lifts of i, found from exact arcs and again from the darts alone. Summed over the 168 sheets the crossings are 168 times the linking numbers of Ghys's Lorenz template (336, 336, 1008 for cat × Sphereland 10, cat × Pell, Sphereland 10 × Pell). Of the 276 pairs among the 24 geodesics, 103 are unlinked, and homology accounts for little of the rest. | MP13; `mp13_tg9_verify.py` (8 checks), `test_core_tg9.js` (4), `mp13_lorenz.py` (5) |
| II.9 Why are two of Sphereland 10's lines of length 30 blind to 143a1? | open | Lines 22 and 24 of MP12's table have no 143a1 component and no symmetry found forces it: the Atkin–Lehner involutions carry them to lifts of other motions, and the finer transfer and the mirror leave one integer free. In the labels ε₁₁, ε₁₃ of MP12 the four components are −ε₁₃(1 − ε₁₁)/2 times (2, −1). Lead [hypothesis]: a formula of Waldspurger's kind for the characters of the ring class group of conductor 143 of ℤ[√3]. | MP12 (the miss; H1 to H3) |
| II.10 Why do the Pell lines meet everything with algebraic intersection zero? | answered [derived; the counts EXACT] (T-G10, 12 November 2026) | In MP13 the signed crossings of each Pell line with all 24 geodesics add up to zero, though the line of length 143 crosses every other one 6 to 506 times; the signed matrix has rank 12 on a span of rank 14. Answer, in three steps, each checked on 38 geodesics. (1) The mirror τ ↦ −τ̄ reverses the surface's orientation, so two geodesics it reverses have signed count zero: the class of a reversed geodesic lies in the mirror's −1 eigenspace, which is isotropic. Line 25 and 9 more of the 24 are reversed. (2) The other 14 have shadows of length 1 on X₀(11) or X₀(13), and the Pell shadow misses them, so line 25 never meets the 14 (the plan's conjecture, made with low confidence, held). (3) Lines 26 to 28 are ±X in homology (MP12, H1), so their rows are ± line 25's. The zero belongs to that set of 24: UV⁵'s four lines of length 7 cross line 25 once each, signed ±1. MP13's lead named a reflection; the one that does the work is the mirror itself, and τ ↦ 1/(143τ̄) was not needed. | MP14; `mp14_tg10_verify.py` (A0 to A4) |

## III — Gravity as a square (the coin) (partial)

**Question.** Roger's coin: are electromagnetism and gravity rewrites of each other?

**Verdict.** In established physics, yes, as a square and not as a flip (8D Box 6.2 for Kerr–Schild; The Modular Plane 12.3 for radiation): a scalar, light and linearised gravity sit on one null vector as S = 1/d, F = n/d³, W = n ⊗ n/d⁵ with S·W = F ⊗ F, and n·σ = 8ψψ̃ᵀ is one spinor squared. Kepler's orbit is the square of an oscillator in the same sense (12.2). The square forgets a sign in both. What is missing is any modular origin for the square: nothing on X₀(143) is yet known to do the squaring. 14 November 2026: the squared hopfion has a frame frozen into its null flow (III.4), whose helicity is two thirds of the light's (III.5); the square holds for the fields and not for these potentials (N-G16).

**Records.** RD_Notes §35; 8D Box Room 6 and its scorecard; 12.2, 12.3; `mp8_verify.py` (B, C); `frozen_gravity_20261110.py` (B). MP23, `mp23_tg1.py` (T-G1, T-G2).

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| III.1 Black holes as squares | answered [structure] | Kerr in Kerr–Schild form is the double copy of Coulomb's field moved to z = ia, with g = 2 (Monteiro, O'Connell and White 2014). | 8D Box 6.2 |
| III.2 Radiation as a square | answered [structure] | The hopfion triple on one spinor; linearised vacuum equations only. Known physics (Swearngin et al. 2013; Thompson, Swearngin and Bouwmeester 2014). | 12.3; C1 to C7 |
| III.3 Kepler as a square | answered [structure] | z = w²; G = MMᵀ has trace 2a, det b², discriminant (2c)²; e = tanh d; one revolution is half a turn of the oscillator. | 12.2; B1 to B6 |
| III.4 A full metric for the squared hopfion | answered [derived; tested in space; unreviewed] (T-G1, 14 November 2026) | The tetrad e^μ = dx^μ + B^μ, B_μ = Re[(i a² (∂_μb) dγ − b² (∂_μa) dβ)/d²] with γ = 1 − α, is frozen into the null flow: k^μ𝓕^a_μν = 0 to 3·10⁻¹⁶ of |𝓕|, and the curvature of h = B + Bᵀ is the squared hopfion to 2·10⁻¹⁶ of |W|, at 60 points. Its field strength is 𝓕_μ = Re[G_μ dα ∧ dβ], G_μ = (i/d)(a ∂_μb − b ∂_μa): the light's own 2-form, with weights carried by the rays. It is the only solution among holomorphic weights up to degree 4. The square of the potentials, d·σ ⊗ σ, is not a potential for it (N-G16). | MP23; `mp23_tg1.py` |
| III.5 Its gravitational helicity against the light's | answered [derived; tested in space; unreviewed] (T-G2, 14 November 2026) | In the frame of III.4, for W = n ⊗ n/d⁵ and F = n/d³: H_g = ∫ B_a ∧ dB^a = (2/3)π² and H_em = π², both constant in time, equal on the label sphere and on two slices in space to nine figures. The ratio is 2/3. The number belongs to that frame (IV.2), and it has no mechanism beyond one integral (AG-D5). A hand value of 3/4 was registered and missed (N-G16). | MP23 |
| III.6 Does the square have a modular origin? | open [hypothesis] | A statement on the curve that sends a weight-k object to its square with gravitational meaning. None is known; 6.1's flip τ → −1/τ is modular, the square is not (the scorecard's "mixed"). 12 November 2026: the hermitian square ψψ† is the Hopf map at each position of Room 8's dial, the plane is its first rung (X = ggᵀ), and the ray that carries the knotted field is ψψ† = (ρ/2)(1 + v·σ) for the same spinor whose other square is the field. That places the square on the dial, not on the curve: still open. | 8D Box scorecard; `sphere_charts_20261112.md` §2 |

## IV — The frozen field (partial)

**Question.** Asenjo, Winkler and Comisso (2026) show that the tetrad's field can be frozen into a flow, with a conserved helicity of twist, writhe and linking. Is that the survival into 3 + 1 of §27's knot reading?

**Verdict.** A candidate, with two weaknesses on record. Their helicity is the translational Chern–Simons form e_a ∧ de^a, so the 2 + 1 Lagrangian reappears in 3 + 1 as a charge that counts linking while the ideal condition holds. But the charge belongs to a choice of frame (flat space in a twisted frame carries it), and in an observer's own freely falling frame the condition holds only where there is no tide. The identification Roger proposed on 10 November, that the oscillator's hidden circle and the frozen lines are the same objects, does not hold as stated: they are the two sides of one matrix multiplication, mirror families linked −1 and +1 (12.1). T-G6 (11 November) adds a caution of another kind: the tidal lines of force of the squared field coincide with the light field's carried lines only at t = 0, so the null flow that carries the light's lines does not carry gravity's. Those are not the papers' lines, which need a tetrad (T-G1). 14 November 2026: T-G1 supplies such a tetrad for the squared hopfion, frozen into the null flow (III.4), and T-G8 shows that the hopfion's later lines of force are integrable, on tori between two closed lines (IV.9).

**Records.** `frozen_gravity_20261110.md` and `.py` (39 checks); 12.1, 12.4, 12.5, 12.8; `mp8_verify.py` (A, D, E); MP08; `mp81_verify.py` and MP09 (T-G5); `mp10_tg6_verify.py` and MP10 (T-G6); 12.9, `mp82_verify.py` and MP11 (T-G6 on the page). MP21 and `mp21_tg7.py` (T-G7); MP22 and `mp22_tg8.py` (T-G8); MP23 (T-G1, T-G2).

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| IV.1 What does the ideal condition say geometrically? | answered [derived] | The field is Lie-dragged by the flow: a normal form e^a = dχ^a(s, q) + a^a_i(q)dq^i. | note A1 |
| IV.2 Is the helicity a property of the spacetime? | answered [derived] — a negative | No: of the frame. Minkowski space with a triad turning along z is ideal and has helicity density −2θ′(z). | note A2 |
| IV.3 What breaks the freeze for an observer? | answered [derived] | The rate of strain R of the observers; Ṙ + R² = K with K the tidal matrix of 1.5, so rigidity cannot last in a tide. Rain into a black hole: R = (⅓, −⅔, −⅔)/T, Kasner's exponents at every radius. | 12.5; E1 to E7 |
| IV.4 Are the hidden circle and the field lines the same thing? | answered [structure] — a negative | No: left and right flows of SU(2), linked +1 and −1. They share the group, not the circles. | 12.1; A5, A6 |
| IV.5 Which frames and flows freeze the field of a real spacetime? | open | A worked example beyond static frames: none exists in the papers or here. 14 November 2026: one worked example now exists, in linearised theory and for a null flow: the squared hopfion (III.4, MP23). A full spacetime and a timelike flow are still wanted. | — |
| IV.6 Waves: frozen into what? | open (answered in part by T-G1) | A plane-fronted wave's frame is transported only by a null vector, which the papers' timelike assumption leaves out; extend their theorem or find the obstruction. 14 November 2026: T-G1 gives a radiating, twisting example frozen into a null flow (III.4). Whether it is frozen into any timelike flow is not tested, and the theorem is not extended. | note A4 |
| IV.7 What is the 2 × 2 matrix of the field lines round a core line, and the tube's physical twist? | answered [derived; tested on the field in space; unreviewed] | Elliptic for every pair (p, q): in the torus frame a rotation by 2πν a circuit, ν = √(2pq(p + q)), the identity exactly when 2pq(p + q) is a perfect square (the hopfion; (1,8); for p = 1 Pell's equation). The lines keep Re(α^p β^q) and lie on tori nested round the core line. Against the Seifert framing they wind ν − pq times a circuit at every time (1, 1.7460, 1.7460, 1.8322 for (1,1), (2,3), (3,2), (2,5)); the physical twist is that minus the writhe. This is the light field's tube. 14 November 2026 (T-G7, MP21): the shear at the core line is dι/du = ν[3pq(p + q) − 2(p² + pq + q²)]/(12pq(p + q)), u = 1 − h/h_max, zero only for the hopfion; confirmed on three routes. | 12.8; MP09; `mp81_verify.py` (P1 to P7) |
| IV.8 Does the gravitational wave have such tubes? | answered [derived; tested in space; unreviewed; for the square W = F ⊗ F/S only] | At t = 0, yes: the tendex lines are the light field's E and B lines and the vortex lines are those lines turned by 45° of duality, four families of tubes with the same matrix, ν = √(2pq(p + q)) (a vortex line of (2,3): det 1.00000000, trace −0.0506786). After t = 0 the lines of force are not carried by the light rays: they turn in the screen by arg(d)/2 (8° to 31° along the (2,3) core line at t = 0.8) and are the lines of the closed 2-form λ_t(α)^{1/2} d(α^p) ∧ d(β^q), λ_t = 1/(1 − it(1 − α)). The tubes survive round displaced closed lines, flux-preserving and elliptic, with ν_t = ν + O(t²): 7.745944 at t = 0.1 and 7.751166 at t = 0.8 against 7.745967. Another choice of gravitational (p, q) field would have other lines. | 12.9; MP10; `mp10_tg6_verify.py` (P1 to P5); MP11: the light's carried core line circles the displaced axis on the matrix's ellipse, as J predicts to 1.4 % for t ≤ 0.2 |
| IV.9 What do the gravitational hopfion's lines of force do after t = 0? | answered [derived; tested in space; unreviewed] (T-G8, 14 November 2026) | They are the orbits of a particle in a plane two-well quartic potential: H = |P|² + |s² − σ²|²/t² on H = 1, with s = (1 − it(1 − α))^{1/2}, σ² = 1 − it and P = β̄ (magnetic type) or iβ̄ (electric type). That quartic is the integrable one of ratio 1 : 2 : 1, so there is an exact second invariant I₂, constant along lines integrated in space to 3·10⁻¹⁰. The lines lie on tori nested between two closed lines, an axial great circle of the label sphere and a cross line, linked once; every torus with a rational winding is a torus of closed lines. The winding is W = T_ν/T_μ, two periods of a curve of genus 2, and W − 1 = 5t²/(48(1 + t²)) to leading order on every torus alike; the whole picture depends on t through ε = t/(2√(1 + t²)). The lines are never chaotic. **Two corrections** (N-G15): MP10's "one turn in about 29 circuits" at t = 0.4 was a miscount, the figure is one turn in 63.9 circuits (the sixty returns were 131 circuits); and the line "lost" at t = 0.8 (MP11) keeps its invariant and returns, the search plane was what failed. Open leads: the integer invariants of the closed Legendrian lines (Thurston–Bennequin, rotation), untested; the other pairs (p, q). On the page since version 8.3: instrument 12.10, with 12.9 corrected (MP24). | MP22; `mp22_tg8.py`; earlier 12.9, MP10, MP11; `sphere_charts_20261112.md` §2.2 (P21) |

## V — Breaking and reconnecting (open)

**Question.** What would it take for "gravity can break and reconnect" to be more than the form of the equations?

**Verdict.** A computed example. The magnetic case is on the instruments in its simplest exact form (12.6): two islands, the rate E_X at the X-point, the merger as a pitchfork with the private flux vanishing as (t* − t)², the exponent predicted before the run. The gravitational papers define the corresponding rate and give no spacetime and no number.

**Records.** 12.6; `mp8_verify.py` (F1 to F7); Comisso and Asenjo (2026).

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| V.1 Reconnection where it can be computed exactly | answered [model: magnetic] | A = Σ± (σ₀²/s)exp(−((x ∓ x₀)² + y²)/2s); E_X = 2ησ₀²e^{−x₀²/2s}(x₀² − 2s)/s³; the X-point's generator goes hyperbolic, parabolic, elliptic. | 12.6 |
| V.2 A three-dimensional example with helicity | open (pre-registered T-G3) | Two linked rings with resistivity: dH/dt = −2∫E·B to tolerance, and the linking lost per unit of the papers' rate. | T-G3 |
| V.3 A gravitational example | open | A spacetime, a frame and a flow in which the non-ideal term acts and the connectivity of the tetrad's field changes. | — |
| V.4 Anything observable | open | The papers' energy release is a possibility they raise, not a result; no signature is derived. | — |

## VI — The field equations (open)

**Question.** Does the programme have its own field equations, and do they have solutions?

**Verdict.** Proposed, unsolved. The compilation's modular Einstein equation with discrete labels n has no Schwarzschild or Friedmann solution on record (V.D). What the instruments carry is ordinary general relativity read on the strip: geodesic deviation as Room 1 with κ = −K (1.5) and the strain rate as the point that matrix moves (12.5).

**Records.** Compilation V.D; 1.5; 12.5; RD_Notes (the elasticity entry of 22 October, part iv); `poisson_mass_20261022.py`.

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| VI.1 Solve the modular field equations for a point mass | open | A solution, its Newtonian limit, and what the label n does to it. | Compilation V.D |
| VI.2 Solve them for a homogeneous universe | open | The Friedmann case, against VII. | Compilation V.D |
| VI.3 The tide on the strip | answered [structure] | ξ̈ = Kξ is Room 1's equation with κ = −K; the catalogue runs from the Earth's surface to Λ. | 1.5 |
| VI.4 The strain rate on the strip | answered [structure] | Ṙ = κ − R²: a Möbius action of exp(τ[[0, κ], [1, 0]]), with de Sitter on the attracting fixed point. | 12.5 |
| VI.5 Gravity's kinetic term as an elastic medium | answered [structure] | The ADM kinetic term is the σ₂ part of an isotropic quadratic energy (DeWitt's supermetric, the labile aether); F_P = c⁴/G is its stiffness scale. | RD_Notes, 22 October |

## VII — The vacuum and the expansion (open)

**Question.** Can the cosmological constant and the expansion history come out of the programme?

**Verdict.** DIAGNOSTIC at best. The dark-energy formula ρ_Λ/M_P⁴ = δ⁻⁶e^{−2/α} waits for a route from an action (CC-34), and its exponent is itself undecided (compilation V.E). The matter programme's own terms give a classical runaway and a one-loop energy with no minimum (VE-01; matter map I.6). What 12.5 shows about the expansion is standard cosmology with Planck's numbers, not an MTFT result.

**Records.** Compilation V.E; CC34_DARK_ENERGY_FORMULA_AUDIT.md; VE01_VACUUM_ENERGY.md; matter map I.6, I.7; 12.5; 3.2.

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| VII.1 The dark-energy formula from an action | open | A derivation, or its retirement to PHENO. | CC-34 |
| VII.2 The exponent: 4, 6 or the fold formula | open | One of them from a calculation rather than a fit. | Compilation V.E |
| VII.3 A vacuum energy from the parent action | open | A stationary point in the volume; today there is none. | VE-01; matter map I.6 |
| VII.4 Where the universe sits on the line of strain rates | answered [structure] | de Sitter is the fixed point R/√κ = 1; today 1.3769 with Ω_Λ = 0.685; matter era √2 on the squeezing branch. Inputs are Planck 2018's. | 12.5; E5 |
| VII.5 A background H(t) for the big-bang-to-matter run | open | The gate of PG-01: a declared parent, a stabilised volume, then standard physics. | matter map I.7 |

## VIII — Dark matter (open)

**Question.** Is the programme's information-geometric dark matter, G_eff(r)/G_N = 1 + η|R_F(r)|, a prediction?

**Verdict.** Not yet: η is free and R_F has not been computed for a real mass distribution (compilation V.B). By the standard of the programme's rule of 17 October on new particles (say why it has not been seen), a proposal with no particle should say what it predicts where particle dark matter and modified gravity differ.

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| VIII.1 Fix η from the arithmetic and compute R_F for a galaxy | open | A rotation curve with no fitted parameter. | Compilation V.B |
| VIII.2 A sealed test | open | A criterion registered before the comparison: rotation curves, lensing of a merging cluster, or the CMB's third peak. | — |

## IX — Gravitational waves (open)

**Question.** Does the programme say anything about gravitational waves that general relativity does not?

**Verdict.** Not so far. The compilation's corrections to the primordial spectrum are uncomputed (V.C). The knotted waves on the instruments are solutions of linearised general relativity, known since 2013.

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| IX.1 The primordial spectrum's modular corrections | open | A number. | Compilation V.C |
| IX.2 Knotted gravitational radiation | answered [structure] | The (p, q) fields and their squares; core lines linked +pq along B, carried along light rays; at t = 0 the ribbon's twist is −√(pq) and the writhe √(pq) − pq [derived]. The squares' own lines of force coincide with the light's only at t = 0 (IV.8). The ribbon's −pq is the Hopf invariant of f = α^p/β̄^q, whose fibres the core lines are: the class of f in π₃(S²) = ℤ, checked three ways (12 November 2026); in four space dimensions the corresponding group is ℤ₂ and only the parity would survive. | 12.3, 12.4; D1 to D8; MP10; `sphere_charts_20261112.md` §2.3 (P7 to P10) |
| IX.3 Anything specific to X₀(143) | open | A reason for particular (p, q), or a selection rule; (2, 3) on the instruments is our choice, and every torus knot is available. | — |

## X — Black holes (open)

**Question.** Paper 33 identifies black holes with the cusps of X₀(143). Does that identification compute anything a horizon is known to do?

**Verdict.** Not audited under the present protocol. Paper 33's computation is that the arithmetic entropy S(y) = Σ n w_n e^{−2πny} grows as C/y² toward a cusp; reading that as the area law needs the map from y to r, which the paper itself lists as postulated (its open problem i), and the constant C = 1/4G rests on Paper 31's formula, which is a model (I.2). On the instruments there are two exact, standard facts about real black holes: the tide (1.5) and the rain frame's Kasner strain (12.5). One overlap is worth a test: the four cusp widths 143, 13, 11, 1 that Paper 33 calls four horizon classes are the cycle type of T on P¹(ℤ/143), and of the Pell unit (II.3). MP14 (T-G10, 12 November) adds two facts [structure]. Paper 33 §5.4's zero period of 143a1 from ∞ to the cusp 0 is (i/2π)L(143a1, 1) up to sign, and it vanishes because W₁₄₃ reverses that path and 143a1 has root number −1. And the four cusps scatter into one another with a matrix that is unitary on the line Re s = 1/2 and equals −I at the bottom of the continuum: whatever goes into a cusp comes back out of the cusps. A horizon absorbs; whether that difference counts against the identification is a lead [hypothesis], not a result.

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| X.1 The cusp–Schwarzschild map y ↔ r | open | A derivation; today it is postulated through the η bridge. | Paper 33 §9.3 i |
| X.2 Rain into a hole | answered [structure] | R = (½, −1, −1)√(2GM/r³) = (⅓, −⅔, −⅔)/T; from a solar-mass horizon T = 4GM/3c³ = 6.57 µs. | 12.5; E1 to E6 |
| X.3 Do the cusps and the 3-fold's parabolic lines say the same thing? | open [hypothesis] | Paper 33's horizon classes against 12.7's cycles of T. A note for the audit: the lab and Sphereland use the standard labels (the cusp ∞ of Γ₀(143) has width 1, the cusp 0 width 143); Paper 33 §2.3 has those two labels the other way round. The set of widths is the same. | Paper 33 §2.3; 12.7 |

## XI — Is gravity quantum? The low-mass frontier (open)

**Question.** Tabletop experiments are pushing gravity towards masses small enough to be put in superposition. Does anything in the programme meet them?

**Verdict.** Not yet. An et al. (2026) measure Newton's law between two 7.6 mg wafers, classically and to 5 %; the quantum experiment they point to needs masses 8 × 10⁸ times lighter. What the instruments already hold is the mathematics of that experiment: its witness is 10.3's determinant, its entangling phase is the tide of 1.5 and 12.5, and its comparison with Coulomb is the square of piece III. All three are established physics. The programme itself predicts nothing here: no short-range force law and no statement on whether gravity entangles.

**Records.** `gravity_milligram_20261111.md` and `.py` (5 checks); An et al., Nat. Commun. (2026), doi:10.1038/s41467-026-78187-4.

| Sub-puzzle | Status | Answer, or what would settle it | Records |
|---|---|---|---|
| XI.1 Does Newton's law hold between milligram masses? | answered [structure: experiment] | Yes, from 0.654 to 1.554 mm: 471 ± 23 aN at the closest, G = 6.67(35) × 10⁻¹¹, no dependence on the modulation frequency. A null test. | An et al. 2026 |
| XI.2 What would show that gravity is quantum at these energies? | answered [structure: a proposal, not done] | Entanglement of two masses by gravity alone (Bose et al.; Marletto and Vedral, 2017). On the lab's matrix: ψ = ½[[e^{iφ_LL}, e^{iφ_LR}], [e^{iφ_RL}, e^{iφ_RR}]], and the concurrence is twice the modulus of det ψ, which equals the modulus of sin(Δφ/2). Whether entanglement proves quantisation is itself debated. | note; 10.3 |
| XI.3 What sets the entangling phase? | answered [derived here; standard physics] | The tide: Δφ = −(m Δx_A Δx_B t/ħ) e_Aᵀ K e_B with K = (Gm/d³)(3nnᵀ − I), the tidal matrix of 1.5. For the proposal's values Δφ = −0.314 rad, C = 0.156. By the coin, Coulomb does the same with e²/4πε₀ for Gm²; equal at the Stoney mass, 1.86 µg. | note P1 to P4 |
| XI.4 Does the programme predict a deviation (α, λ) these balances could see? | open [hypothesis] | A force law at 10 to 100 µm with a mechanism, sealed before the data. The dark-energy length (ħc/ρ_Λ)^{1/4} ≈ 88 µm is a coincidence of scales until then (AG-D5; VII.1). | — |
| XI.5 Does the programme say whether gravity entangles? | open [hypothesis] | A quantum version of the square with a modular origin (III.6). Nothing today. | — |

## Pre-registered tests

| Test | Registered | Criterion | Outcome |
|---|---|---|---|
| T-G1 A metric perturbation for the double-copy hopfion, frozen into its null flow | 2026-11-10 (note) | h_μν with k^μ𝓕^a_μν = 0 to machine precision, or a proof none exists | **done 2026-11-14** (MP23, predictions before the run): exhibited. The tetrad of III.4 is frozen into k to 3·10⁻¹⁶, its curvature is the squared hopfion to 2·10⁻¹⁶ and its Ricci tensor 4·10⁻¹⁶, at 60 points. Found by a fit in the class 𝓕 = Re[G dα ∧ dβ] with holomorphic G, where it is unique up to degree 4. One forecast missed (the square of the potentials), one held in part |
| T-G2 Its gravitational helicity against the electromagnetic one | 2026-11-10 (note) | Both numbers in the frame of T-G1 | **done 2026-11-14** (MP23): H_g = (2/3)π² and H_em = π², constant in time. **The registered hand value 3/4 missed** on an algebra slip, which is recorded |
| T-G3 A three-dimensional resistive example | 2026-11-10 (note) | dH/dt = −2∫E·B to tolerance on two linked rings; linking lost per unit rate | pending. A lead by hand, not checked (MP23): with a uniform resistivity the first half is an identity between two double line integrals; the second half needs a localised resistivity and traced field lines |
| T-G4 The linking invariant at level 143 | 2026-11-10 (note); predictions in MP12, written before the run | What replaces Ghys's Rademacher function on Γ₀(143)\PSL(2, ℝ) | **done 2026-11-12**: the class of the closed line in ℤ³⁰ (Ψ, three windings, 26 periods as 2 + 4 + 8 + 12), computed for 12.7's 28 lines; none bounds. 25 of 26 pre-registered checks on Route A, 8 on mtft's own model, 5 on the page's core; **one pre-registered conjecture missed and recorded** (Atkin–Lehner signs account for eight of the ten lines with no 143a1 component, not for the other two); one fault of the plan's sample recorded (P9) |
| T-G5 The 2 × 2 matrix of field lines round a core line | 2026-11-10 (MP08, added after the run; predictions in MP09) | Elliptic, parabolic or hyperbolic for (p, q) ≠ (1, 1); it gives the tube's physical twist | **done 2026-11-10**: elliptic, a rotation by 2π√(2pq(p + q)); the identity for the hopfion; winding √(2pq(p + q)) − pq against Seifert. 10 checks on the field in space, no prediction missed; one printed decimal of the plan was wrong and is recorded (MP09) |
| T-G6 The gravitational tube | 2026-11-10 (MP09, opened by T-G5's result) | The return matrix of the lines of force of W round the core lines, against ν = √(2pq(p + q)) | **done 2026-11-11** (predictions in MP10, written before the run): at t = 0 the same matrix, for tendex and vortex lines alike; afterwards the lines of force are not carried by the light rays, and the tubes persist, displaced, with det 1 and ν_t = ν + O(t²). 7 checks in space passed; **two pre-registered windows missed and recorded** (the size of the term after t², for the shift of ν and for the hopfion's defect) |
| T-G7 The shear at the core line | 2026-11-10 (MP09) | dι/dh at the core line in closed form, checked against the profile ι(h) | **done 2026-11-14** (MP21, predictions before the run): dι/du = ν[3pq(p + q) − 2(p² + pq + q²)]/(12pq(p + q)) at the core line, u = 1 − h/h_max. Confirmed by the profile integral (5·10⁻¹⁰ of ν), by lines in space at two times (3·10⁻⁸) and by the page's own function (3·10⁻⁸). No prediction missed; three faults of mine recorded (two hand decimals of the plan, a root finder, a wrong first repair) |
| T-G8 The gravitational hopfion after t = 0 | 2026-11-11 (MP10, opened by T-G6's result) | An exact invariant of its lines of force at t > 0, or evidence that there is none; the closed lines; the slow winding in closed form against the measured one turn in about 29 circuits at t = 0.4 | **done 2026-11-14** (MP22, predictions before the run): an exact invariant exists (IV.9). Thirteen checks held; one forecast missed (how far the line at t = 0.8 travels: 4.7, not over 10); one clause of an added test printed MISS on a fault of the test (its counter counted the starting point). Two earlier statements corrected (N-G15) |
| T-G9 The linking numbers of the lines that bound | 2026-11-12 (MP13, predictions before the run; opened by T-G4) | Lk of the unions K_γ = γ̂ ∪ γ̂ʳ from the crossings of the geodesics on the curve, against a formula on the darts and against Ghys's Lorenz template through the 168-sheeted cover | **done 2026-11-12**: Lk(K_γ, K_δ) = −I(γ, δ) ≤ 0; the cat block 0, −6, −8, −8, −6, −24; totals 336, 336, 1008 = 168 × the template's 2, 2, 6. 8 checks on Route A (exact rational arithmetic), 4 on the page's core, 5 on the template; nothing missed |
| T-G10 The Pell zeros, and the gap at infinite horodepth | 2026-11-12 (MP14, predictions before the run; opened by II.10 and by Roger's idea of 3 November) | A mechanism for the Pell lines' zero signed counts, tested on new closed geodesics (the mirror rule, the Pell classes crossing by crossing, a non-zero pairing, the shadows one level down); the scattering matrix of Γ₀(143), the threshold s = 1/2 and the Eisenstein series on the Pell line, with a second route that sums over the group | **done 2026-11-12**: the mirror, shadows that miss one level down, and homology account for the zeros (II.10). Φ(1/2) = −I, so every Eisenstein series vanishes at the bottom of the continuum; on the Pell line E_∞ = E₀ and E_{1/11} = E_{1/13}; the gap is exactly 1/4, set by the cusps, with nothing below it (a theorem at this level, from the literature). 9 checks on Route A, 7 on Route B; nothing missed, the low-confidence conjecture included; two faults of the apparatus recorded and repaired (a residue check's precision; rounding near two cusps) |

Sealed and run in MP08 (10 November): linking +1 for left flows and −1 for right flows, pass; the private flux exponent 2 at the merger, pass (2.022); the cycle types at 143 for T, S, the order-3 motion, the cat map and the Pell unit, pass; the sign of the ribbon's linking in 12.4, **miss** (predicted +, found −).

Written and run in MP10 (11 November): the field equations and eigen-directions of the square, the t = 0 tubes, the turn by arg(d)/2, the displaced tubes, the hopfion's circles breaking, all pass; the window for the shift of ν between t = 0.1 and 0.2, **miss** (ratios 1.60, 1.09, 2.58 against 0.75 to 1.25); the window for the hopfion's defect, **miss** (2.80 against 3.2 to 4.8).

Written and run in MP12 (12 November): the graph and its homology, the three units and their windings, Ψ and its congruence modulo 6, the transfer, reciprocity for the cat map, the Fricke involution on the Pell unit, the Manin model and its four blocks, the periods of 143a1, the rank 30 and the index 6, the agreement with mtft's model, the ghost of the Pell lines (P12) and the Atkin–Lehner signs on the blocks (P13), all pass; the conjecture that every line with no 143a1 component is negated by an Atkin–Lehner involution, **miss** (eight of ten). The plan's six small elements could not carry P9 and twelve were used: a fault of the sample, recorded. Six post hoc checks (H1 to H5) are counted apart.

Written and run in MP13 (12 November): the exact cutting sequences of the three orbits, the 24 geodesics behind the 28 lines, the convention-free crossing counts, the cat block from the darts, the totals from the Lorenz template, the homological relations among the signed counts, and the invariance under the mirror and the Atkin–Lehner involutions, all pass. Nothing missed; a sign error in the frame at the lifts of ρ was caught on reading, before the first run.

Written and run in MP14 (12 November): the mirror rule on 14 new geodesics, the Pell classes crossing by crossing, a non-zero pairing of the Pell class, the shadows on X₀(11) and X₀(13), the scattering matrix (residues, functional equation, unitarity, −I at s = 1/2, the Dirichlet series), the vanishing at the threshold, the Pell line, and Route B's sums over the group, all pass. Nothing missed. B1's first run failed on the precision of its residue check (40 digits where 80 were needed), and Route B's label check, added with the route, on rounding near the cusps 6/11 and 6/13; both were repaired and recorded. The plan's amendment carried a time not read from the clock (12:10; the file times say 12:08:49), corrected by appending.

## Honest negatives

| | Piece | What was closed | Record |
|---|---|---|---|
| N-G1 | I | G is not derived; Paper 31's formula is a model, +3.7 % or +61 % by the depth | RD_Notes §27 iv |
| N-G2 | I | In 3 + 1 dimensions gravity has local degrees of freedom and no Chern–Simons level: the 2 + 1 chain does not carry over by itself | RD_Notes §27 |
| N-G3 | III | The coin holds for gravity as a square, not as a flip; "gravity on the strong force's face" is mixed | 8D Box scorecard |
| N-G4 | IV | The oscillator's hidden circle and the frozen field's lines are not the same objects: mirror families | 12.1; note |
| N-G5 | IV | The frozen-in helicity is a property of the frame, not of the spacetime | note A2 |
| N-G6 | IV | In the observers' own freely falling frame the freeze cannot hold where there is a tide | 12.5 |
| N-G7 | IV | The pre-registered sign of the ribbon's linking was wrong: −pq, as the 22 October record for the trefoil's torus framing already had it (−6.00000) | MP08; RD_Notes, 22 October |
| N-G8 | V | No gravitational reconnection has been computed, in the papers or here | 12.6 |
| N-G9 | IV | The lines of force of the squared field W = F ⊗ F/S are not frozen into the flow of light rays: they agree with the light field's carried lines only at t = 0 | MP10 (P3) |
| N-G10 | IV | Two pre-registered windows in T-G6 were wrong: the shift of ν and the hopfion's defect do start at order t², but the next term is many times larger than assumed | MP10 (P4w, P5w) |
| N-G11 | XI | The milligram experiment tests nothing the programme predicts: no force law at short range, no statement on entanglement; its 5 % value of G adds nothing to CODATA for Paper 31's formula | `gravity_milligram_20261111.md` |
| N-G12 | I, III, IX | Sphereland's two charts, joined to the matrix algebra, organise and predict nothing for the programme: no pair (p, q), no dimension, no G. Three space dimensions are singled out twice (an integer knot charge for lines; the oscillator's square is Newton's potential) and derived by neither, and the dial still offers four space-times | `sphere_charts_20261112.md` |
| N-G13 | II | The pre-registered conjecture of T-G4 was wrong: Atkin–Lehner signs do not account for every closed line with no 143a1 component. They account for eight of the ten. The involutions carry the other two, lines of length 30, to lifts of other motions, and those two zeros have no mechanism (II.9) | MP12 (P13; H2) |
| N-G14 | II | MP12's Result, and II.8 as first written, said that all the lifts of one motion together bound. Only the cat map's do: the fibre numbers of the lifts add up to 168Ψ(A), so for Sphereland 10's motion and the Pell unit the union is homologous to −28 and −588 fibres. Corrected by appending, in MP12 and II.8 | MP12 (correction); MP13 |
| N-G15 | IV | Two statements about the gravitational hopfion were wrong. "One turn in about 29 circuits" at t = 0.4 (MP10) was a miscount: the return function took the crossing nearest the start within an arc of 60, ten times the eighth one, so sixty returns were 131 circuits; the figure is one turn in 63.9. And the line "lost" at t = 0.8 (MP11) was not lost: it keeps its invariant and returns, and it was the search plane that its loops left for an arc of 69.8. The "half a turn round the mean" printed at t = 0.1 and 0.2 was the angle of a short arc seen from its own middle. | MP22 |
| N-G16 | III | The coin does not work at the level of potentials: the square d·σ ⊗ σ of the light's potential over the scalar is not a potential for the squared hopfion; its own Ricci tensor is 5.9 times |W| (forecast at 0.75, missed). And the hand value H_g/H_em = 3/4, registered before the integral was computed, was wrong: it is 2/3. | MP23 |

## Chronology

| Date | Record | Piece | What it added |
|---|---|---|---|
| 2026-04 | Paper 33, Black Holes as Modular Cusps | X | The cusp identification and S(y) ~ C/y² |
| 2026-10-21 | RD_Notes §27; Sphereland 18 | I, II | Kepler's cusp and the trefoil; the oscillator; G as a level in 2 + 1; the audit of Paper 31 |
| 2026-10-22 | RD_Notes (elasticity); Sphereland 20 | VI | DeWitt's supermetric as the labile aether; the trefoil's ledger on a torus |
| 2026-10-27 | RD_Notes §35; 8D Box Room 6 | III | The coin: the flip, the double copy, one face heavier |
| 2026-11-08 | The Modular Plane 1.5 | VI | The tidal reading |
| 2026-11-10 | `frozen_gravity_20261110.md` | IV | The papers read; 39 checks; T-G1 to T-G4 |
| 2026-11-10 | The Modular Plane Room 12; MP08 | II to V, IX, X | Seven instruments; 51 checks and one recorded miss; T-G5 |
| 2026-11-10 | The Modular Plane 12.8 (version 8.1); MP09 | IV | T-G5 answered: the flux tube's matrix is elliptic, ν = √(2pq(p + q)); T-G6 and T-G7 opened |
| 2026-11-11 | MP10; `mp10_tg6_verify.py` | IV, IX | T-G6 answered: the gravitational square has the same tubes at t = 0, and afterwards its lines of force are not carried by the light rays; two windows missed and recorded; T-G8 opened |
| 2026-11-11 | `gravity_milligram_20261111.md` | XI, I, III | An et al. 2026 read: Newton holds between 7.6 mg masses; the entangling phase of the proposed quantum test is the tide of 1.5; nothing for the programme to claim |
| 2026-11-11 | The Modular Plane 12.9 (version 8.2); MP11 | IV | T-G6 on the page: the turn, the displaced closed line, its return matrix, and the returns of the light's core line on the matrix's ellipse (6 checks in space, 13 by the page's core; nothing missed) |
| 2026-11-12 | `sphere_charts_20261112.md`; the glossary's section on Sphereland's two charts | I, II, III, IV, IX | The holes chart and the sphere table joined to the matrix algebra: the charge rule, the Hopf ladder on the 𝕂 dial, the ray as the spinor's hermitian square, the ribbon's integer as a Hopf invariant, the strange numbers as Bernoulli numbers. Leads for I.6, T-G4 and T-G8; N-G12. 49 checks in space, 8 on the page's core, 9 on mtft; one pre-registered clause missed and recorded (the sign convention of the test's linking integral) |
| 2026-11-12 | MP12; `mp12_tg4_verify.py`, `mp12_tg4_mtft.py`, `test_core_tg4.js` | II | T-G4 answered: thirty integers for each closed line of 12.7 in place of the Rademacher number, sorted by level and by Hecke block; none of the 28 lines bounds. One conjecture missed and recorded (N-G13); II.8 and II.9 opened |
| 2026-11-12 | MP13; `mp13_tg9_verify.py`, `test_core_tg9.js`, `mp13_lorenz.py` | II | T-G9 answered: a line with its reverse bounds, and two such pairs link by minus the crossings of their geodesics on the curve. The cat block from exact arcs and from the darts; the totals 168 times Ghys's Lorenz template. 17 checks, nothing missed; a correction to MP12 (N-G14); II.10 opened |
| 2026-11-12 | MP14; `mp14_tg10_verify.py`, `mp14_lattice.py` | II, X | T-G10 answered: the Pell zeros come from the mirror, from shadows that miss one level down, and from homology (II.10). Alongside, the gap of X₀(143) at infinite horodepth: the scattering matrix of the four cusps, −I at the bottom of the continuum, where every Eisenstein series vanishes, and the Pell line as the nodal line of the σ-odd half. 16 checks on two routes, nothing missed; two faults of the apparatus recorded. Paper 33 §5.4's zero period explained (X) |
| 2026-11-14 | MP21; `mp21_tg7.py`, `mp21_core.js` | IV | T-G7 answered: the shear of a flux tube at its core line in closed form, zero only for the hopfion; three routes; no prediction missed |
| 2026-11-14 | MP22; `mp22_tg8.py`, `mp22_after.py`, `mp22_p7_check.py` | IV, IX | T-G8 answered: the gravitational hopfion's lines of force are the orbits of an integrable quartic; the invariant, the two closed lines, the winding as a ratio of periods; MP10's 29 circuits and MP11's lost line corrected (N-G15); one forecast and one clause of an added test missed. Put on the page the same day as instrument 12.10, version 8.3 (MP24: 8 checks on the page's core, nothing missed) |
| 2026-11-14 | MP23; `mp23_tg1.py`, `mp23_kernel.py` | III, IV | T-G1 answered: a tetrad for the squared hopfion frozen into its null flow, in closed form, to machine precision. T-G2 answered: H_g = (2/3)π² against H_em = π². One forecast missed (the square of the potentials) and one registered hand value missed (3/4 for 2/3) (N-G16). T-G3 left pending |
