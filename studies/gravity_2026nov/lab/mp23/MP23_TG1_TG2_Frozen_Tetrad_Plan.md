# MP23 · T-G1 and T-G2: a tetrad for the squared hopfion, frozen into its null flow; the two helicities

14 November 2026. Written at 09:11:21 New York time by the clock, read by the command that wrote this file, before any run. The criteria, registered on 10 November in `frozen_gravity_20261110.md`:

- **T-G1.** "A metric perturbation h_μν for the double-copy hopfion with k^μ𝓕^a_μν = 0 for the null flow k: exhibit it to machine precision, or prove none exists."
- **T-G2.** "Its gravitational helicity in that frame against the hopfion's electromagnetic helicity."

Status of what follows: derived here by hand, unreviewed, prior art not checked. Nothing below has been computed by code.

Records: `mp23_tg1.py`, `mp23_numbers.json`.

## Setting (the note's, and MP08's conventions)

The tetrad is a potential, 𝓕^a = de^a. To first order e^a = dx^a + B^a, where B_aν = ½h_aν + ω_aν: the symmetric part is the metric perturbation and the antisymmetric part a local Lorentz rotation. The squared hopfion is W = n ⊗ n/d⁵ = E + iB (12.3), a linearised vacuum field of type N whose repeated null direction is the light field's flow k = (1, v). Labels α = 1 − 2ib/d, β = 2a/d, with a = x − iy, b = t − i − z, d = r² − (t − i)².

## The derivation

1. **A change of coordinates does not touch 𝓕.** Under x → x + ξ, B^a changes by dξ^a. So the frozen condition is a condition on the frame ω alone, and the freedom is six functions against twelve equations: a solution, if there is one, is special.
2. **What a frozen tetrad must look like.** i_k𝓕^a = 0 with d𝓕^a = 0 says each 𝓕^a is carried by the rays: it is a closed 2-form on the label sphere, pulled back (the note's A1). A class of these is plain: **𝓕^a = Re[G^a(α, β) dα ∧ dβ] for four holomorphic functions G^a**, with potential B^a = Re[P^a dα + Q^a dβ], ∂Q^a/∂α − ∂P^a/∂β = G^a. Every such tetrad is frozen into k, because k·∂α = k·∂β = 0. The metric is h_μν = B_μν + B_νμ.
3. **Why a solution should exist for type N.** Take a null frame (k, m) whose Cartesian components are constant along each ray, V^{ab} = k^[a m^b], and ϑ = m·dk. The spin connection ω^{ab} = 2 Re(Φ V^{ab} ϑ) gives a frozen tetrad for any function Φ, and its curvature is 2 Re[V^{ab} (dΦ + 2iΦς) ∧ ϑ], with dm = iς m − ϑ v. That has the type-N shape when Φ is constant along the rays and holomorphic in the sense of the label sphere. So frozen tetrads go with type-N fields along a shear-free congruence; whether the hopfion's is among them is the test.
4. **A candidate from the coin.** The light field has the potential A = σ/2 with σ = α dβ − β dα, and the scalar is S = 1/d. The square "A ⊗ A/S" is the complex tensor H_μν = d·σ_μσ_ν. By hand, d·σ_μ = N_μ − 2βx_μ exactly, with N = 2(0, α, −iα, −iβ) and x_μ = (−t, x, y, z). The first piece, N ⊗ σ, is of the frozen class with G = 3N.

## Predictions, forecasts and what counts as a miss

- **Q0 (checks of the setting, EXACT or to 10⁻³⁰).** k·∂α = k·∂β = 0; ∂α and ∂β are null and orthogonal; d·σ_μ = N_μ − 2βx_μ; W is trace-free and W·v = 0.
- **Q1 (forecast, 0.75).** The curvature of Re(c·d·σ ⊗ σ) is the squared hopfion for one constant c: Ricci-flat, and E + iB proportional to n ⊗ n/d⁵ (or to its conjugate, by the sign convention of B) with the same constant at every point. If this misses, the square does not hold at the level of these potentials.
- **Q2 (the search; forecast 0.65 that it succeeds at degree ≤ 2).** In the holomorphic class of step 2 with G^a a polynomial of degree ≤ 2 in (α, β), 24 complex unknowns, a real least-squares fit of the curvature to the squared hopfion at 40 random points of space-time. Success: a residual below 10⁻¹⁰ of |W|. If it fails at degree 2 the search goes to degree 4 (60 unknowns). A failed search is not a proof that none exists; T-G1 would then stay open with the class recorded as empty.
- **Q3 (the exhibit, if Q2 succeeds).** For the tetrad found, computed again from B itself by a separate routine at 60 fresh points and t from −1 to 2: |i_k𝓕^a| < 10⁻¹² of |𝓕|; the curvature of h equals the squared hopfion to 10⁻¹⁰ of |W|; its Ricci tensor vanishes to 10⁻¹⁰. **That is T-G1's criterion.**
- **Forecast on the form (0.4).** G is of degree 1 and contains N.
- **Q4 (T-G2, if T-G1 is exhibited).** The helicity 3-form of the frame is K = e_a ∧ de^a. Its first-order part is exact; the second-order part is B_a ∧ dB^a, contracted with the Minkowski metric. For a frozen tetrad of the holomorphic class it is a 3-form on the label sphere, so the helicity H_g = ∫ B_a ∧ dB^a over a slice does not depend on the time. Predictions: H_g computed on the slices t = 0 and t = 0.8 in space agrees with the integral over the label sphere to 10⁻⁶; the light's helicity for F = Re(dα ∧ dβ) is **H_em = π²** (by hand: A ∧ dA = ½ of the sphere's volume form, and the volume is 2π²), checked the same way.
- **Q5 (the caution of IV.2, tested).** The helicity belongs to the frame. If the fit of Q2 has a kernel (frozen frames of this class with no curvature), H_g is reported as a function on that kernel and not as one number.

## Amendment 1, after Q0 to Q2 and before Q3 and Q4

Written at 09:15:55 New York time by the clock. Q0, Q1 and Q2 have been run once: Q1's forecast missed, Q2's search succeeded at degree 2 with a single solution (rank 48 of 48, residual 3.4·10⁻¹⁵ of |W|). The full result is written below after Q3 and Q4. This amendment records what I derive from the fitted coefficients by hand, before the exhibit and the helicities are computed.

- **The fitted G in closed form.** The 24 coefficients are, to the printed six decimals, G_μ = −½(1 − α)∂_μa + (i/2)β ∂_μb, with the two constant null covectors ∂a = (0, 1, −i, 0) and ∂b = (1, 0, 0, −1). With γ = 1 − α = 2ib/d that is G_μ = (i/d)(a ∂_μb − b ∂_μa), and a potential for it is

  **B_μ = ¼ Re[γ² (∂_μa) dβ + iβ² (∂_μb) dγ] = Re[(i a² (∂_μb) dγ − b² (∂_μa) dβ)/d²],  e^μ = dx^μ + B^μ.**

  Q3 is run on this closed form.
- **The gravitational helicity, by hand.** B_μ ∧ dB^μ = ⅛ Re[|γ|²γ dβ ∧ dγ̄ ∧ dβ̄], because ∂a·∂a = ∂b·∂b = ∂a·∂b = ∂b·∂ā = 0 and ∂a·∂ā = 2. On the sphere dβ ∧ dᾱ ∧ dβ̄ = 2ᾱ·vol, so the 3-form is −¼ Re[(1 − α)²(ᾱ − |α|²)]·vol, and only the term −3|α|² survives the integral: **H_g = (3/4)π²**, with the orientation in which **H_em = π²**. Prediction for Q4: the ratio H_g/H_em is 3/4, to 10⁻⁶, on the label sphere and on the slices t = 0 and t = 0.8.
- **The normalisation.** Re(dα ∧ dβ) has E + iB = i·n/d³ if n/d³ = ∇α × ∇β, so its helicity is that of the light field F = n/d³ (a duality turn does not change it). Checked in Q4.

## Result

Written at 09:19 New York time on 14 November, by the clock. Numbers are in `mp23_numbers.json`.

**T-G1 is answered: the tetrad exists, and it is exhibited to machine precision. T-G2 is answered: the two helicities are (2/3)π² and π².** One forecast missed (Q1), one prediction missed on a hand slip of mine (the ratio 3/4 in Amendment 1; it is 2/3), and one forecast held only in part. Status: derived here, tested in space, unreviewed; prior art not checked.

### The tetrad

**e^μ = dx^μ + B^μ,  B_μ = Re[(i a² (∂_μb) dγ − b² (∂_μa) dβ)/d²],**

with a = x − iy, b = t − i − z, d = r² − (t − i)², β = 2a/d, γ = 1 − α = 2ib/d, and the two constant null covectors ∂a = (0, 1, −i, 0), ∂b = (1, 0, 0, −1). Its field strength is

**𝓕_μ = dB_μ = Re[G_μ dα ∧ dβ],  G_μ = −½(1 − α)∂_μa + (i/2)β ∂_μb = (i/d)(a ∂_μb − b ∂_μa).**

So each of the four field strengths is the light field's own 2-form, dα ∧ dβ, times a weight that is constant along every ray and linear in the labels. The metric perturbation is h_μν = B_μν + B_νμ.

### The checks

- **Q0** held: the labels are carried by k, their gradients are null and orthogonal, d·σ_μ = N_μ − 2βx_μ, and W is trace-free and across v, all to 10⁻⁴⁰. The curvature routine reproduces a plane wave to 10⁻¹⁵ and gives 1.2·10⁻¹⁰ for a pure gauge term (nested differences).
- **Q1: the forecast missed** (0.75). The square of the potentials, d·σ ⊗ σ, is not a potential for the squared hopfion. Its own Ricci tensor is 5.9 times |W|, and the best constant leaves 0.41 of |W|. The coin does not work at the level of these potentials.
- **Q2: the search succeeded, as forecast** (0.65), at degree 2. The fit leaves 3.4·10⁻¹⁵ of |W| and has full rank, 48 of 48: the solution is the only one of its class. The 24 coefficients are the G above; those of degree 2 are zero.
  - The sign of the magnetic part is the convention B_ij = ½ε_jkl R_0ikl = −Im W. With the other sign there is no solution (0.49 of |W|).
  - **Forecast on the form, held in part** (0.4): G is of degree 1, as forecast. It is not built on N: it is built on the two constant null covectors. I count the forecast as missed.
- **Q3 held: T-G1's criterion.** At 60 fresh points with t from −1 to 2, from the closed form by a separate routine: |i_k𝓕| is 3.4·10⁻¹⁶ of |𝓕|; the curvature of h equals the squared hopfion to 2.0·10⁻¹⁶ (E) and 1.3·10⁻¹⁶ (B) of |W|; its Ricci tensor is 4.2·10⁻¹⁶; R·k is 5.9·10⁻¹⁶ of |R|. The fitted potential and the closed form have the same 𝓕 to 7.6·10⁻¹⁵.
- **Q4: MISS on the ratio; the rest held.**
  - H_em = π² to 1.6·10⁻⁹, on the label sphere and on the slices t = 0 and t = 0.8 in space. ∇α × ∇β = n/d³ to 4·10⁻¹⁶, so this is the helicity of the light field F = n/d³.
  - H_g = 6.579736, the same on the sphere and on both slices to nine figures. **That is (2/3)π², not the (3/4)π² of Amendment 1.**
  - **The slip.** In Amendment 1, |1 − α|²(1 − α)ᾱ is (1 − α)²(ᾱ − ᾱ²); I wrote (1 − α)²(ᾱ − |α|²). The terms that survive the integral are −2|α|² − |α|⁴, with means ½ and ⅓ over the sphere, which gives −(8/3)π² and H_g = (2/3)π². The method was right and the algebra was wrong; the prediction as registered missed.
- **Q5, the caution.** The fit has no kernel at degree 2. Run after Q3 and Q4 and not predicted (`mp23_kernel.py`): at degree ≤ 4 it still has full rank, 120 of 120, and returns the same G with every higher coefficient below 10⁻¹⁵. So within the holomorphic class the frozen frame is unique as far as degree 4, and H_g is one number there.

### What it says

- **III.4, T-G1.** The squared hopfion has a frame that is frozen into its null flow. The frame's field strength is carried by the rays exactly as the light's is.
- **III.5, T-G2.** In that frame H_g/H_em = 2/3, and both are constant in time. For W = n ⊗ n/d⁵ and F = n/d³: H_g = (2/3)π² = 6.5797, H_em = π² = 9.8696.
- **The caution of IV.2 stands.** The helicity belongs to the frame. Frozen frames outside the holomorphic class are not excluded, and another frame would give another number.
- **IV.6, in part.** The papers' theorem assumes a timelike flow. Here is a worked radiating example with a null one. Whether the squared hopfion is frozen into any timelike flow is not tested.
- 2/3 comes out of one integral over the label sphere. It has no mechanism beyond that, and I attach no meaning to it (AG-D5).

### Not done

- A proof that the frame is unique among all frozen frames.
- The other pairs (p, q). By step 2 the same class is there; no fit was run.
- **T-G3 (the three-dimensional resistive example) is not attempted here; it stays pending.** A lead, by hand and not checked: with a uniform resistivity and no flow, two rings with Gaussian cores stay two Gaussian rings, the helicity is 2Φ₁Φ₂ times a Gauss linking integral smoothed over the summed variance S = s₁ + s₂ + 4ηt, and dH/dt = −2∫E·B becomes an identity between two double line integrals. That would do the first half. The second half, the linking lost per unit of the papers' rate, needs the largest voltage along a field line through a bounded resistive region, which means a localised resistivity and traced field lines: a simulation, not a closed form.
- Housekeeping: importing the script for the kernel run wrote an empty entry into the numbers file, which that run removed.
