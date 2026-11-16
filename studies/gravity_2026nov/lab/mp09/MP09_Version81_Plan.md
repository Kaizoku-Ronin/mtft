# MP09 · The Modular Plane 8.1 · T-G5: the 2 × 2 matrix of field lines round a core line

Pre-registration written 10 November 2026, 23:35 New York time, before any run. T-G5 was stated in MP08 after its run: the linearised field-line flow in a frame carried round a core line of the (p, q) field, is it elliptic, parabolic or hyperbolic? For the hopfion it should be the identity. For (p, q) ≠ (1, 1) it would give the physical twist of the flux tube, which the torus ribbon of 12.4 does not.

## Conventions (those of MP08)
a = x − iy, b = t − i − z, d = r² − (t − i)², α = 1 − 2ib/d, β = 2a/d, F_pq = pq α^{p−1} β^{q−1} n/d³ = E + iB. Hopf coordinates on the labels: α = cos η e^{iξ₁}, β = sin η e^{iξ₂}, and Θ = pξ₁ + qξ₂. The core lines of 12.4 are cos²η = p/(p + q) with Θ = 0 and Θ = π. "One circuit" is φ: 0 → 2π along (e^{iqφ}α, e^{−ipφ}β).

## The derivation (by hand, to be tested; status: derived here, unreviewed, prior art not checked)
1. The flux 2-form of B on a slice of constant time is the imaginary part of d(α^p) ∧ d(β^q) pulled back by the labels. A field line is the kernel of that 2-form, which needs no metric, and the 2-form is the same on the label sphere at every time. So the lines are carried with the labels (12.4's ride) and everything below is independent of t.
2. On the label sphere the kernel is, up to the positive factor pq cos^{p−1}η sin^{q−1}η,
   η′ = −sin Θ sin η cos η,  ξ₁′ = cos Θ sin²η,  ξ₂′ = −cos Θ cos²η,  hence  Θ′ = cos Θ (p sin²η − q cos²η).
3. H = cos^p η sin^q η cos Θ = Re(α^p β^q) is constant along every magnetic line (and Im(α^p β^q) along every electric line). The level sets of H are tori nested round the two core lines, where H is extreme: ±h_max, h_max = (p/(p + q))^{p/2} (q/(p + q))^{q/2}. I believe this invariant is already in Kedia et al. (2013); not checked.
4. Linearised at the core line Θ = 0, in the frame (δη, δΘ), which is the torus frame (δη = 0 is the neighbour on the same torus, the ribbon of 12.4) and closes after one circuit:
   **L = [[0, −s], [2(p + q)s, 0]], s = √(pq)/(p + q)**: trace 0, det L = ω² = 2pq/(p + q) > 0, an elliptic generator.
   One circuit takes the parameter time T = 2π(p + q), so the matrix of one circuit is
   **M = exp(TL) = [[cos 2πν, −sin 2πν/√(2(p + q))], [√(2(p + q)) sin 2πν, cos 2πν]], ν = √(2pq(p + q))**.
5. Against the Seifert framing the torus frame itself turns −pq times (12.4's ribbon), so neighbouring lines wind round the core line **ν − pq = √(2pq(p + q)) − pq** times per circuit; and against the rotation-minimising frame, by Călugăreanu, **Tw = ν − pq − Wr(t)**. The sign + is fixed by the hopfion, whose magnetic lines link +1 (12.1).
6. Off the core line, on the torus H = h, the winding per circuit in the torus frame is
   **ι(h) = 2πq / ∫ du/√(e^{λ(u) − 2 ln h} − 1)**, u = ln cos²η between the two roots of λ(u) = 2 ln h, λ = ln(x^p (1 − x)^q).
   Limits by hand: ι → ν at the core line; ι → 2pq at the edge H = 0; for the hopfion ι = 2 on every torus (the integral is elementary).

## Predictions and what counts as a miss
- **P1 (EXACT + NUMERICAL)** H is constant along a magnetic line integrated in space, at t = 0 and at t = 0.8, to 1e-8. Miss: drift above that.
- **P2 (kind)** The return map round the core line, from the field in space (variational equation along the closed line), has determinant 1 and trace 2 cos(2π√(2pq(p + q))): 2 for (1,1) (the identity), −0.0507 for (2,3) and (3,2), 0.9868 for (2,5); the same at t = 0 and t = 0.8, to 1e-6. So: elliptic for the three knotted pairs, the identity for the hopfion. Miss: any trace outside the tolerance, or a hyperbolic one.
- **P3 (rotation number in the torus frame)** ν = +2, +7.745967, +7.745967, +11.832160, right-handed in the lab's chart, counted continuously in the closed torus frame. Miss: another integer part, or the sign −.
- **P4 (the physical numbers)** Winding against the Seifert framing ν − pq = 1, 1.745967, 1.745967, 1.832160, at both times (the torus frame's own linking is the integer −pq of 12.4, recounted). Physical twist Tw = ν − pq − Wr: at t = 0 (Wr = √(pq) − pq) 1, 5.296477, 5.296477, 8.669882; at t = 0.8 (MP08's writhes) 1.176192, 5.751232, 5.703596, 9.289322.
- **P5 (profile)** ι(h) from lines integrated in space at h/h_max = 0.9, 0.5, 0.1 agrees with the formula of step 6 to 1e-6; ι is 2 at all three levels for the hopfion; ι(10⁻⁶ h_max) is within 2 % of 2pq. Miss: any of these.
- **P6 (when the neighbours close)** M = I exactly when 2pq(p + q) is a perfect square. For p = 1 that is Pell's equation (2q + 1)² − 2m² = 1: q = 1, 8, 49, 288. Among coprime p ≤ q ≤ 60 I expect only (1,1), (1,8), (1,49) and (32,49). Test: brute force, and the trace for (1,8) from the field in space is 2 to 1e-6 with ν = 12 (winding 4 against Seifert). Miss: another pair in range, or (1,8) not the identity.
- **P7 (a closed line off the axis)** For (2,3) the torus with ι = 8 carries closed magnetic lines: integrated in space from that torus a line returns to its start after one circuit (to 1e-6) and links the core line +2 (integer, by crossings, both along B). Miss: no closure, or a linking number other than +2.
- Stated, not tested: the torus-frame matrix is first order. M = I for (1,8) means the neighbours close to first order only, because ι(h) is not constant there (it runs from 12 to 16).

Amendment before the run (10 November, about 23:40): where M is the identity the angle is recovered from the trace as a square root, so for (1,1) and (1,8) P3 and P4 are tested to 1e-4, not 1e-6.

## Instrument 12.8, "The flux tube's matrix" (version 8.1)
Chips (1,1), (2,3), (3,2), (2,5), (1,8); sliders: the torus (h/h_max), t, view. Figures: the core line in space with one magnetic line on the chosen torus; the cross-section (cos²η across, Θ up) with the level curves of H; the profile ι(h) with ν and 2pq marked. The strip reads M. Routes: A `mp81_verify.py` (space side, numpy/scipy/sympy, reuses the helpers of `mp8_verify_part1.py`), B `build/test_core81.js` on the page's core against A's numbers, and the browser look test.

## Result (10 November 2026, 23:50)
No prediction was missed. One printed decimal of the plan was wrong and is recorded below; the tolerance amendment above was made before the run.

- **Route A** (`mp81_verify.py`, in space, on the field itself): 10 checks passed, 0 failed (E1 to E3 exact, P1 to P7).
- **P1** drift of H along an integrated line 3.6e-10 of h_max (four pairs, t = 0 and 0.8).
- **P2** traces of the return map from the variational equation: 2.0000000, −0.0506786, −0.0506786, 0.9871990 at both times, against the formula's 2, −0.0506786, −0.0506786, 0.9871991; worst |det − 1| 2.5e-8, worst closure of the core line 3.1e-10. Elliptic for (2,3), (3,2), (2,5); the identity for the hopfion. **Slip, recorded**: the plan printed 0.9868 for (2,5). The formula it stands for, 2 cos(2π√140), is 0.9871991; my hand arithmetic was wrong in the fourth decimal. The prediction as a formula holds; the printed decimal stays in the text above.
- **P3** ν = 2.000000 (t = 0) and 1.999996 (t = 0.8, an identity case, tolerance 1e-4), 7.745967, 7.745967, 11.832160, counted continuously in the closed torus frame (means 2.0000, 7.7460, 11.8322): the sign is +.
- **P4** the torus frame links −1, −6, −6, −10 (−8 for (1,8)); winding against Seifert 1.000000, 1.745967, 1.745967, 1.832160 at both times; twist at t = 0: 1.000000, 5.296477, 5.296477, 8.669882; at t = 0.8: 1.176189, 5.751232, 5.703596, 9.289322 (the hopfion's is 3e-6 from the pre-registered 1.176192, inside the identity-case tolerance).
- **P5** lines integrated in space against the integral: worst 2.0e-10; 2 on every torus of the hopfion; at 1e-6 h_max: 2.0, 11.9994, 11.9994, 19.9798 (2pq = 2, 12, 12, 20); ι falls from the edge to the axis without a turn on a grid of 40 levels.
- **P6** perfect squares for coprime p ≤ q ≤ 60: (1,1), (1,8), (1,49), (32,49), as predicted; Pell q = 1, 8, 49, 288; for (1,8) the return map from the field has trace 2.00000000 and |M − I| = 7.1e-8, with ν = 12 and winding 4 against Seifert.
- **P7** the torus ι = 8 of the (2,3) field is at h/h_max = 0.800844; a line integrated in space from it closes to 1.4e-10 after one circuit and links the core line +2 (16 crossings; Gauss integral 2.00000000).
- **Repairs during the run that touch no prediction**: the profile integral divided by zero at a level within 1e-9 of h_max; its last sliver at each end is now taken by its tangent, and the axis check uses 1 − 1e-6.
- **Route B** `test_core81.js`, 17 of 17 on the page's extracted core; `test_core8.js` 44/44, `test_core71.js` 15/15, `test_core7.js` 24/24, `test_core6.js` 30/30. Found by Route B and fixed in the page: the drawn line ended with a linear interpolation (drift in H 5e-5 at the last point); it now ends with a true step (6e-9). One test definition corrected: turns are counted at upward crossings.
- **Look test** `look81.py` (desk 1280 × 900, phone 390 × 844, dark): 56 instruments, no console errors, no overflow, no clipped or overlapping labels in Room 12 for the five chips at three levels; slider handlers 13 to 26 ms.
- **Page** version 8.1, 842,872 bytes, sha256 b753d8af989f84a8561934f965a5d99b178243463711372676002178b1353299, artifact version 1794372425-1008.

## What T-G5 now says
The matrix is elliptic for every coprime (p, q): a rotation by 2πν, ν = √(2pq(p + q)), in the torus frame. It is the identity exactly when 2pq(p + q) is a perfect square. The winding against the Seifert framing, ν − pq, does not change in time; the physical twist is ν − pq − Wr(t). Off the core line the winding rises to 2pq at the edge, and on the tori where it is a whole number m the lines close, linked m − pq times with the core line. Status: derived here, tested in space by two routes, unreviewed; prior art for ν not checked (the invariant H is very likely in Kedia et al. 2013).

## Open after this
- **T-G6** the gravitational tube: W of 12.3 is built from the same labels. Do its lines of force (the eigen-directions of its electric and magnetic parts) form tubes round the same core lines, and with what matrix?
- **T-G7** second order: the shear dι/dh at the core line in closed form, and whether (1,8)'s first-order closure has any consequence.
- A literature check for ν = √(2pq(p + q)).
