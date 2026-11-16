# MP21 · T-G7: the shear at the core line, in closed form

14 November 2026. Written at 08:27:06 New York time by the clock, read by the command that wrote this file, before any run. T-G7 was opened in MP09: "the shear dι/dh at the core line in closed form, and whether (1,8)'s first-order closure has any consequence." The puzzle map's criterion is "dι/dh at the core line in closed form, checked against the profile ι(h)".

Status of what follows: derived here by hand, unreviewed, prior art not checked. Nothing has been computed by code.

Records: `mp21_tg7.py`, `mp21_core.js`, `mp21_numbers.json`.

## Conventions (MP09's)

Labels α = cos η e^{iξ₁}, β = sin η e^{iξ₂}, Θ = pξ₁ + qξ₂. A magnetic line of F_pq keeps H = cos^p η sin^q η cos Θ. The core lines are cos²η = p/(p + q), Θ = 0 or π, where H = ±h_max, h_max = (p/(p + q))^{p/2} (q/(p + q))^{q/2}. On the torus H = h a line winds ι(h) times round the core line per circuit, in the torus frame; ι → ν = √(2pq(p + q)) at the core line and 2pq at the edge. Write u = 1 − h/h_max, so u = 0 is the core line.

## The derivation

1. **The cross-section is a Hamiltonian system.** With x = ln cos²η and ξ₁ as the time, MP09's line field gives dx/dξ₁ = 2 tan Θ and dΘ/dξ₁ = p − q cos²η/sin²η. These are Hamilton's equations for the pair (Θ, x) and the function 𝓗 = 2 ln H = λ(x) + 2 ln cos Θ, λ(x) = px + q ln(1 − eˣ).
2. **So the period is the derivative of an area.** One turn round the core line takes ξ₁ the time P = dA/dε, where A(ε) is the area of the cross-section inside the level 𝓗 = 𝓗_max − ε, and ε = 2 ln(h_max/h) = 2u + O(u²). One circuit is ξ₁: 0 → 2πq, so ι = 2πq/P.
3. **The area to second order.** Near the core line ε = a y² + b y³ + c y⁴ + Θ² + Θ⁴/6 + …, with y = x − x₀ and a, b, c from the second, third and fourth derivatives of λ at x₀ = ln(p/(p + q)):
   a = p(p + q)/(2q), b = p(p + q)(2p + q)/(6q²), c = p(p + q)((p + q)² + 4p(p + q) + p²)/(24q³).
   For ε = a y² + b y³ + c y⁴ + Θ² + βΘ⁴ the area is A = (πε/√a)·[1 + ε·(15b²/(32a³) − 3c/(8a²) − 3β/8)] + O(ε³).
4. **The result.**

   **dι/du at the core line = ν · [3pq(p + q) − 2(p² + pq + q²)] / (12pq(p + q)),  so  dι/dh = −(that)/h_max.**

   It is zero for the hopfion, where ι = 2 on every torus, and positive for every other pair: for p = 1 the bracket is (q + 2)(q − 1).

| (p, q) | ν | dι/du at the core | dι/dh at the core |
|---|---|---|---|
| (1,1) | 2 | 0 | 0 |
| (1,2) | 3.464102 | 0.192450 | −0.500000 |
| (2,3), (3,2) | 7.745967 | 1.118862 | −6.018516 |
| (2,5) | 11.832160 | 1.859339 | −15.2404 |
| (1,8) | 12 | 0.972222 = 35/36 | −4.4244 |

The last column is the hand value; the first two columns are what is tested.

## What it says about (1,8)

For (1,8) the matrix of one circuit is the identity, because ν = 12 is a whole number. The shear is 35/36, not zero. So a line at the small level u turns 2π·(35/36)·u round the core line per circuit, to first order in u: the neighbours close to first order in their distance and miss at third order. Nothing else follows inside this field, because H is exact and every line stays on its torus.

Recorded as a hypothesis, not tested here (**H-G7**): the core line of (1,8) is degenerate, since its return map is the identity, so the usual argument that a closed line survives a small change of the field does not apply to it. With a shear that is not zero I expect the closed lines near it to move by the cube root of a perturbation's strength, against the first power for the pairs with ν not a whole number.

## Predictions and what counts as a miss

- **E1 (EXACT)** The equations of step 1 follow from MP09's line field, and are Hamilton's equations for 𝓗 = 2 ln H, for symbolic p and q.
- **E2 (EXACT)** The coefficients a, b, c of step 3 and the bracket of step 4, for symbolic p and q.
- **P1 (the profile integral)** (ι(u) − ν)/u from MP09's integral at 30 digits, at u = 10⁻³, 2·10⁻³ and 4·10⁻³, extrapolated to u = 0, equals the closed form to 10⁻⁸ of ν for (1,2), (2,3), (3,2), (2,5), (1,8), (3,5) and (1,49); for the hopfion ι − 2 is zero to 10⁻²⁰ at the three levels. Miss: any pair outside that.
- **P2 (lines in space)** The same slope from magnetic lines integrated in space on the field itself, at t = 0 and t = 0.8, at u = 0.004, 0.008 and 0.016, extrapolated: within 10⁻⁵ of ν for (2,3), (3,2), (2,5) and (1,8); for the hopfion |ι − 2| < 10⁻⁸. Miss: any outside that, or a difference between the two times above the tolerance.
- **P3 (the page)** The page's own profile function, taken from version 8.1's core, gives the same slope to 10⁻⁶ of ν for the five pairs of 12.8. Miss: any outside that.
- **F1 (forecast, not a test of the formula)** The next term is of ordinary size: |ι(u) − ν − S·u| < 3·S·u² at u = 0.016 for the four knotted pairs of P2, S being the closed form (0.8).

## Result

Written at 08:34 New York time on 14 November, by the clock. **The closed form holds on every route. No prediction was missed.** Three faults of mine are recorded below; none touches a prediction.

- **E1, E2** exact, for symbolic p and q: the cross-section is a Hamiltonian system for 2 ln H; the coefficients a, b, c; the level curve to second order and its area; the bracket of the closed form.
- **P1** held. From the profile integral at 30 digits the slope agrees with the closed form to 4.9·10⁻¹⁰ of ν at worst over the seven pairs (tolerance 10⁻⁸). The hopfion's ι − 2 is 5·10⁻²⁵. Taken from levels a hundred times closer to the core line (not a registered test), the slope agrees to 2.7·10⁻¹⁵ for (2,3) and 1.9·10⁻¹⁵ for (1,8).
- **P2** held. From magnetic lines integrated in space on the field itself, the slope agrees to 2.6·10⁻⁸ of ν at worst (tolerance 10⁻⁵), and the two times t = 0 and t = 0.8 agree to 1.0·10⁻⁹. The hopfion's |ι − 2| is 7·10⁻¹²; H drifts by 3·10⁻¹² along the lines.
- **P3** held. The page's own profile function, taken as it stands from version 8.1's core, gives the slope to 2.7·10⁻⁸ of ν at worst over the five pairs of 12.8 (tolerance 10⁻⁶).
- **F1** held (forecast 0.8). The next term is 0.615, 0.615, 0.625 and 0.567 times S·u² for (2,3), (3,2), (2,5) and (1,8). No closed form is claimed for it.

| (p, q) | ν | dι/du at the core | dι/dh at the core |
|---|---|---|---|
| (1,1) | 2 | 0 | 0 |
| (1,2) | 3.464102 | 0.192450 | −0.500000 |
| (2,3), (3,2) | 7.745967 | 1.118862 | −6.018519 |
| (2,5) | 11.832160 | 1.859339 | −15.092000 |
| (1,8) | 12 | 35/36 = 0.972222 | −4.671936 |

### Faults, recorded

1. **Hand slips in the plan's last column.** The plan printed −15.2404 for (2,5) and −4.4244 for (1,8); the formula gives −15.0920 and −4.6719. My arithmetic for h_max was wrong. For (2,3) the sixth decimal was wrong (−6.018516 for −6.018519). The plan said that column was a hand value and not what is tested; the printed decimals stay in the text above.
2. **The root finder.** The first run of P1 stopped on a type error, before any number was produced. The cause: for (1,8) at u = 0.004 my first guesses sent both roots of the level to the same one. The roots are now bracketed, one each side of the core line.
3. **A wrong first repair.** Before finding that, I guarded the two ends of the integrand (the last sliver by its tangent, as in MP09), taking rounding there to be the cause. It was not; the second run stopped the same way. The guard stays in and is harmless.

### What T-G7 now says

- The winding of a flux tube rises away from its core line at the rate dι/du = ν·[3pq(p + q) − 2(p² + pq + q²)]/(12pq(p + q)), u = 1 − h/h_max. Status: derived here, tested by three routes, unreviewed; prior art not checked.
- It is zero only for the hopfion. So every knotted tube has shear at its core, which is the condition under which most of its tori survive a small change of the field (the twist condition of the standard theorem; not tested here).
- For (1,8) the neighbours of the core line close to first order and miss at third: a line at level u turns 2π·(35/36)·u per circuit. H-G7, the cube-root response of that core line to a perturbation, stays a hypothesis.
- Not done: MP11 hoped the shear would turn the reported errors of its P3 (the matrix's miss at t = 0.4 and 0.8) into predictions. That needs the level u of the carried core line inside gravity's displaced tube, which I have not computed.
