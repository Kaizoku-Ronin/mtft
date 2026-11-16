# MP22 · T-G8: the gravitational hopfion's lines of force after t = 0

14 November 2026. Written at 08:47:55 New York time by the clock, read by the command that wrote this file, before any run. The puzzle map's criterion for T-G8: "An exact invariant of its lines of force at t > 0, or evidence that there is none; the closed lines; the slow winding in closed form against the measured one turn in about 29 circuits at t = 0.4."

Status of what follows: derived here by hand, unreviewed. Nothing below has been computed by code. Before writing this I read MP10's and MP11's results and scripts, and I used three recorded numbers of MP10 (defect/t² = 0.196, 0.191, 0.179 at t = 0.0125, 0.025, 0.05); where that bears on a prediction it is said.

Records: `mp22_tg8.py`, `mp22_numbers.json`.

## What MP10 established, and where this starts

At time t the lines of force of the gravitational hopfion are, on the label sphere |α|² + |β|² = 1, the lines of the closed 2-form λ_t(α)^{1/2} dα ∧ dβ, with λ_t(α) = 1/(1 − it(1 − α)): its imaginary part for the magnetic type (the lines of Im G, G = √d·F₁₁), its real part for the electric type.

## The derivation

1. **A change of label turns the 2-form into dA ∧ dβ.** Put s = (1 − it(1 − α))^{1/2}, the branch with s = 1 at α = 1. The branch point is α = 1 + i/t, outside the sphere for every t, so s is single-valued on it. With A = 1 + 2(s − 1)/(it), dA = λ_t^{1/2} dα, and the old label is a quadratic polynomial of the new one: **α = A + (it/4)(A − 1)²**.
2. **The lines of force are the orbits of a particle in a plane.** The lines of Im(dA ∧ dβ) on the hypersurface |α(A)|² + |β|² = 1 are its characteristic lines, which are the orbits of the Hamiltonian |β|² + |α(A)|² on its level 1, with Im(dA ∧ dβ) as the symplectic form. In the variable s:

   **H = |P|² + |s² − σ²|²/t²,  σ² = 1 − it,  on H = 1,**

   with the position s in the plane and the momentum P = β̄ for the magnetic type, P = iβ̄ for the electric type. The potential is |α|². It is a quartic with two wells, at s = ±σ; the lines live in the well at s = σ, and the level 1 lies below the saddle, whose height is 1 + 1/t².
3. **That quartic is an integrable one.** Turn the plane so the wells lie on the real axis: ζ = x + iy = s·σ̂*, p = P·σ̂*, σ̂ = σ/|σ|, and m = |σ|² = √(1 + t²). Then V = [r⁴ − 2m(x² − y²) + m²]/t², r² = x² + y². This is the quartic with the ratio 1 : 2 : 1, which separates in elliptic coordinates; here the foci are at ζ = ±2√m, twice as far out as the wells. (That this potential is integrable is classical, Garnier's system; quoted from memory, not checked. The reduction of the lines of force to it is what is new here.) The second invariant, with L = xp_y − yp_x:

   **I₂ = L² + 4m p_x² + (4m/t²) x² (r² − 2m).**

   In the labels, for the magnetic type: I₂ = (Im sβ)² + 4m (Re σ̂β)² + (4m/t²)(Re sσ̂*)²(|s|² − 2m). For the electric type the first two terms are (Re sβ)² and 4m (Im σ̂β)².
4. **The tori and the two closed lines.** With elliptic coordinates ζ = 2√m cosh(μ + iν) and K = −I₂/(4m), the motion separates: p_μ² = 4m·Π(cosh²μ) and p_ν² = −4m·Π(cos²ν), with the cubic Π(X) = K + X − (4/ε²) X (X − ¾)² and

   **ε = t / (2√(1 + t²)).**

   After a scaling the whole problem depends on t only through ε: H = |π|² + |u|²·|1 + εu/2|². Π has three roots 0 < X₁ < X₂ < 1 < X₃; cos²ν swings between X₁ and X₂ and cosh²μ between 1 and X₃. Every orbit is of this one kind, so the label sphere is filled by tori nested between exactly two closed lines:
   - **the axial line**, K = 1/t², y = 0: on the label sphere it is the great circle (α, β) = (−iσ̂² cos φ, σ̂* sin φ);
   - **the cross line**, where X₁ = X₂ = c₀ = ½ − ¼√(1 + 4ε²/3): it crosses the axis at x₀ = 2√(m c₀) at right angles.
5. **The winding.** On the torus K the number of μ-swings per ν-swing is

   **W(K, ε) = T_ν/T_μ,  T_ν = ∫ from X₁ to X₂ of dX/√R,  T_μ = 2 ∫ from 1 to X₃ of dX/√R,  R(X) = X(X − 1)Π(X).**

   These are periods of a curve of genus 2. The slow winding round a closed line is W − 1 turns per circuit.
   - On the axial line: T_μ = (4πε/3)(1 − 4ε²/9)^{−1/2} and T_ν = (ε/2) ∫ dφ / [(¾ − y)(¼ + y)^{1/2}], y = (ε/2) sin φ, φ from −π/2 to π/2.
   - To leading order, by the second-order normal form of the 1 : 1 resonance: the averaged Hamiltonian is (ε²/32)[−12 n_a² − (32/3) n_a n_b + (4/3) n_b²] in the two actions; the terms that would couple the phases cancel. So the two frequencies differ by (5/6)ε² on every torus alike:

   **W − 1 = (5/12) ε² + O(ε⁴) = 5t² / (48(1 + t²)) + …, the same on every torus.**

## What this says about the recorded numbers

- **"One turn in about 29 circuits" should be wrong as a count of circuits.** At t = 0.4, ε² = 0.0345 and the leading law gives W − 1 = 0.0144; on the axial line my hand value of the exact ratio is 1.0156. So one turn takes about 60 to 70 circuits. MP10's sixty "returns" were each the crossing that came closest to the starting point within an arc of 60, which is up to eight circuits; so sixty returns were more than sixty circuits.
- **The defect.** After one circuit a line misses its start, across the line, by 2π(W − 1)·|a||b| on the label sphere, where |a|² = (Im α)² + (Re β)² and |b|² = (Re α)² + (Im β)² are the shares of the two swings (magnetic type, t → 0). In space that is multiplied by (1 + |X|²)/2. For MP10's point (0.6, −0.3, 0.5): |a||b| = 0.36258, so **defect/t² → (5π/24)·0.85·0.36258 = 0.2017 as t → 0**. MP10's three recorded values, 0.196, 0.191, 0.179, extrapolate to about 0.201; I saw them before writing this, so this clause is a check of the derivation and not a blind prediction.
- **The lost line at t = 0.8.** If the invariant holds, that line lies on a torus like every other. MP11's search gave up after an arc of 60 in space; a torus that passes near the label of infinity (α = 1) carries its lines far out.

## Predictions and what counts as a miss

- **E1 (EXACT)** α = A + (it/4)(A − 1)² and (dα/dA)² = 1 − it(1 − α).
- **E2 (EXACT)** {H, I₂} = 0 for symbolic t, with the V and I₂ of step 3.
- **E3 (EXACT)** The normal form of step 5: the coefficients −12, −32/3, 4/3, the cancellation of the coupling terms, and the frequency difference (5/6)ε².
- **P0 (the reduction, in space)** At 60 random points of space and t = 0.1, 0.4, 0.8, 1.5, a step along the field itself (Im G, then Re G) moves (s, P) parallel to (2P, −4(s² − σ²)s*/t²), to 10⁻⁶ in angle. Miss: any point outside that.
- **P1 (the invariant)** Along lines of force integrated in space on the field itself, I₂ changes by less than 10⁻⁸ of its range over 20 circuits: five lines of each type at each of t = 0.1, 0.2, 0.4, 0.8, 1.5, one of them the magnetic line through (0.6, −0.3, 0.5). **Miss: a drift above that. This is the test of "an exact invariant".**
- **P2 (the closed lines)** At t = 0.1, 0.4, 0.8, 1.5 the axial great circle and the cross line, carried to space, are lines of force (tangent to the field to 10⁻⁷) and close (10⁻⁸). They are linked once. On 20,000 random points of the label sphere I₂ lies between its values on the two lines and reaches neither elsewhere. Miss: any of these.
- **P3 (the periods)** Along a line integrated in space, the separated time (dτ′ = 2dτ/|ζ² − 4m|, τ the Hamiltonian time) between successive upward crossings of y = 0 is constant and equals T_μ/(2√m); between successive upward passages of ν through the middle of its range it is T_ν/(2√m); both to 10⁻⁷, for three lines at each of t = 0.1, 0.4, 0.8. So W = T_ν/T_μ is the measured winding. Miss: a period that varies along the line or differs from the integral.
- **P4 (the leading law)** (W − 1)/ε² from the integrals lies within 2 % of 5/12 at t = 0.05 on five tori across the range of K, and W on the axial line at t = 0.4 is 1.0156 ± 0.0002. Miss: either.
- **P5 (MP10's numbers re-read)**
  - (a) defect/t² at t = 0.0125, 0.00625, 0.003125 for MP10's line, by MP10's own return function, extrapolates to 0.2017 ± 0.0010.
  - (b) At t = 0.4 the true first returns of that line turn (W − 1) per circuit, W from its own K: over 200 circuits the turns counted about the centre of the curve of returns equal 200(W − 1) to ±0.05. I expect W − 1 between 0.0125 and 0.0175.
  - (c) MP10's sixty closest-returns took between 115 and 165 circuits, not sixty.
  - Miss: any of (a), (b), (c). If (b) gives 0.034 per circuit the winding formula is wrong.
- **P6 (the lost line)** At t = 0.8 the line through (0.6, −0.3, 0.5), integrated at a speed that stays bounded on the label sphere, keeps I₂ to 10⁻⁸ over 100 circuits and comes back within 0.3 of its start. Forecast, not a test (0.6): it reaches more than 10 from the origin.

## Amendment 1, before the run of P7

Written at 08:57:22 New York time by the clock. Sections E, P0 to P6 have been run; their result is written below, after this amendment and after P7 is run. This adds one test, on the part of the criterion that asks for "the closed lines". Nothing of it has been computed.

- **P7 (a torus of closed lines).** Besides the two closed lines there should be a whole torus of closed lines wherever W is a ratio of whole numbers. At t = 0.4, W runs from 1.01442 to 1.01577 over the tori (P4's table, already computed), so W = 65/64 is taken on one torus. Prediction: a magnetic line started on that torus, integrated in space, closes after 64 circuits to 10⁻⁷, having crossed y = 0 upward 65 times; a line on the torus halfway between that one and the torus W = 66/65 comes no nearer to its start than 10⁻³ in its first 70 circuits. Miss: no closure, closure at another count, or the neighbour closing too.

## Result

Written at 08:59 New York time on 14 November, by the clock. Numbers are in `mp22_numbers.json`.

**T-G8 is answered: there is an exact invariant.** At every time the lines of force of the gravitational hopfion are the orbits of a particle in a plane quartic potential that is integrable. They lie on tori nested between two closed lines, and the winding is a ratio of two periods. One forecast missed and one clause of the added test printed MISS; both are below. Status: derived here, tested in space on the field itself, unreviewed; prior art not checked.

### The checks

- **E1 to E4** exact: the change of label; {H, I₂} = 0 for symbolic t; the normal form and its −12, −32/3, 4/3; W = 1 + (5/12)ε² on the axial line.
- **P0** held. A step along the field moves (s, P) along the Hamiltonian flow to 2.1·10⁻⁸ in angle (60 points, four times, both types); H = 1 to 2·10⁻¹⁴.
- **P1** held. **I₂ is constant along the lines to 2.7·10⁻¹⁰ of its range**: fifty lines integrated in space, both types, t = 0.1 to 1.5, twenty circuits each (tolerance 10⁻⁸).
- **P2** held. The axial great circle and the cross line are lines of force (tangent to 3·10⁻⁸) and close (3·10⁻¹¹); the cross line stays on its hyperbola to 7·10⁻¹³. Their linking number is −1 for the magnetic type and +1 for the electric type, at all four times. On 20,000 random points I₂ lies strictly between its values on the two lines; the nearest came within 2.6·10⁻⁶ of the range of an end.
- **P3** held. The two separated periods are constant along a line to 5·10⁻¹¹ and equal the integrals to 4.7·10⁻¹¹ (nine lines; tolerance 10⁻⁷). So the measured winding is T_ν/T_μ to ten figures.
- **P4** held. At t = 0.05, (W − 1)/ε² is 0.41730, 0.41717, 0.41701, 0.41685, 0.41672 across the tori, against 5/12 = 0.41667. On the axial line at t = 0.4, W = 1.015770; my hand value was 1.0156 ± 0.0002, so it held by 0.00003.
- **P5** held, all three clauses.
  - (a) defect/t² = 0.196321, 0.199046, 0.200386 at t = 0.0125, 0.00625, 0.003125, which extrapolate to **0.201712**. The derivation's limit, from the labels of the point alone, is **0.201712**.
  - (b) At t = 0.4 the line through (0.6, −0.3, 0.5) has W − 1 = 0.015660: **one turn in 63.9 circuits**. Counted over 200 returns: 3.1319 turns against 3.1320.
  - (c) MP10's sixty returns took **131 circuits**: 49 were the first crossing, one the second, ten the eighth. 131 × 0.015660 = 2.05 turns, which is MP10's 2.04.
- **P6** held. At t = 0.8 the line keeps I₂ to 5.7·10⁻¹¹ over 100 circuits and comes back within 0.013 of its start. Its torus is at 0.43 of the range, with W − 1 = 0.0478, one turn in 21 circuits.
  - **The forecast missed** (0.6): I expected it to reach more than 10 from the origin. It reaches 4.7.
- **P7** (Amendment 1): **printed MISS, on a fault of the test.** The line on the torus W = 65/64 closes after 64 circuits to 1.6·10⁻¹¹, and the line halfway to W = 66/65 comes no nearer than 0.0197; those two clauses held. The count of crossings of y = 0 came out 66, not 65. The cause: the starting point lies on y = 0 (y = −1.4·10⁻¹⁶ there) and the integrator counted it as a crossing, at label length 2·10⁻¹⁵. With the start left out the count is 65, the last one at the closure to 5·10⁻¹¹. The MISS stays in the record of the run.

### What was wrong in the earlier records

- **"One turn in about 29 circuits" (MP10; the puzzle map's IV.9 and T-G8) is wrong.** It is one turn in 63.9 circuits. MP10's return function took, from each point, the crossing nearest the starting point within an arc of 60, and ten times that was the eighth crossing. The 2.04 turns were real; the sixty circuits were 131.
- **The line was not lost at t = 0.8 (MP11).** Followed for an arc of 1,500 in space, it crosses MP11's plane 106 times after the start, with gaps of 9.0, 9.1, 9.4, 9.7, 10.2, 11.0, 12.0, 13.4, 15.4, 17.9 and then 69.8. MP11's search looked no further than 60, so it stopped after ten returns. For part of each slow turn the line's loops do not cross that plane. (Found after the run, not predicted; `mp22_after.py`.)
- The "half a turn round the mean" that MP11 printed at t = 0.1 and 0.2 was the angle of a short arc seen from its own middle. The true windings there are one turn in 964 and in 244 circuits (axial line).

### The winding, in numbers

W − 1, turns per circuit, from the axial line to the torus at 0.999 of the way to the cross line:

| t | ε | (5/12)ε² | axial line | middle torus | beside the cross line |
|---|---|---|---|---|---|
| 0.1 | 0.0498 | 0.00103 | 0.00104 | 0.00103 | 0.00103 |
| 0.2 | 0.0981 | 0.00401 | 0.00411 | 0.00406 | 0.00401 |
| 0.4 | 0.1857 | 0.01437 | 0.01577 | 0.01506 | 0.01442 |
| 0.8 | 0.3123 | 0.04065 | 0.05483 | 0.04686 | 0.04109 |
| 1.5 | 0.4160 | 0.07212 | 0.13994 | 0.09521 | 0.07368 |
| 3.0 | 0.4743 | 0.09375 | 0.29630 | 0.13883 | 0.09659 |

- At t = 0.4, 1/(W − 1) runs from 63.4 to 69.4. So the tori with W = 65/64, 66/65, … 70/69 are tori of closed lines, closing after 64 to 69 circuits; P7 followed one.
- As t grows ε tends to ½, where the level 1 reaches the saddle between the two wells. The saddle is the label 1 + i/t, which tends to the label of infinity. On the axial line T_ν then grows without bound, like a logarithm (derived; the table shows the growth, the limit is not tested).

### What T-G8 now says

1. **The invariant.** I₂ = L² + 4m p_x² + (4m/t²)x²(r² − 2m), in the variables of step 3; in the labels it is the expression given there. It is quadratic in β.
2. **The closed lines.** Two for each type at each time: the axial great circle (−iσ̂² cos φ, σ̂* sin φ) and the cross line; linked once. Besides them, every torus on which W is a ratio of whole numbers.
3. **The winding.** W = T_ν/T_μ, two periods of the genus-2 curve Y² = X(X − 1)Π(X); to leading order 5t²/(48(1 + t²)) on every torus.
4. **One parameter.** The picture at time t depends on t only through ε = t/(2√(1 + t²)).
5. The lines of force of the gravitational hopfion are never chaotic.

### Not done

- The page's instrument 12.9 still says that the closed line at the middle of the torus has not been found, and that the line was lost at t = 0.8. I have not changed the page; both sentences should be corrected in its next version.
- The other gravitational (p, q) fields. The same change of label works for q = 1, but the potential is then not a polynomial and I have no invariant for it.
- Prior art: the integrable quartic is classical; whether the reduction is known for the gravitational hopfion I have not checked.
- The lead in IV.9 about the integer invariants of closed Legendrian curves is untouched.
- Housekeeping: importing the script for the two after-run checks wrote an empty entry into the numbers file, which I removed.
