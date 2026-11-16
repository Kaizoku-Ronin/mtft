# `mtft.gravity` — knotted weak fields: the tube, the squared hopfion's lines of force, a frozen tetrad (2026-11-16)

Modules `mtft.gravity.fields`, `knotted`, `hopfion`, `tetrad`; tests `tests/test_gravity_knotted.py`, `test_gravity_hopfion.py`,
`test_gravity_tetrad.py` (23 tests, 12 s). Ported from the Modular Plane lab's records MP09 and MP21 to MP23 (10 and 14 November 2026),
which are kept under `studies/gravity_2026nov/lab/`. Status of the whole: **model**. Every field is a weak wave on flat space (Maxwell's
equations; the linearised vacuum Einstein equations). Newton's constant is not derived anywhere in MTFT (the gravity map's N-G1), and
nothing here bears on it. No link to X₀(143) is claimed.

## 1. Whose work is whose

The release plan asked for a search of the literature before this module shipped, so that nothing is called new that is not. It was made
on 16 November 2026: nine web searches; four papers read through an automated summary of their text and a fifth by its abstract. That
is a search, not a line-by-line reading. "Not found" means not found in these sources.

| In the lab's records | In the literature | Standing |
|---|---|---|
| The (p, q) fields, their invariant Re(α^p β^q), tori nested round two core lines, h_max | Kedia, Bialynicki-Birula, Peralta-Salas, Irvine, Phys. Rev. Lett. 111, 150404 (2013), arXiv:1302.0342: the same, with helicity over energy 1/(p + q); the lines on a torus "either periodic or quasi-periodic" | known |
| The rotation number ν = √(2pq(p + q)) at a core line, the profile ι(h), the shear (MP09, MP21) | no rotation number in Kedia et al., nor in Thompson, Swearngin, Bouwmeester, J. Phys. A 47, 355205 (2014), arXiv:1402.3806 | not found; a few lines from the invariant |
| Gravity's lines of force as the light's lines turned by a local duality angle, vortex lines 45° from tendex lines (MP10) | Swearngin, Thompson, Wickes, Dalhuisen, Bouwmeester, arXiv:1302.1431, equations 28 to 30; Thompson, Wickes, Swearngin, Bouwmeester, J. Phys. A 48, 205202 (2015), arXiv:1411.2073 | known since 2013 |
| The later lines of force as an integrable system: H, I₂, the winding, two closed lines (MP22) | those papers show the later lines in figures, with no invariant, Hamiltonian or winding number; the potential reached is Garnier's system, which separates in elliptic coordinates (arXiv:1305.5686, arXiv:1810.12028, arXiv:math-ph/0212017) | the reduction: not found; the system: classical |
| A tetrad frozen into the null flow (MP23) | Asenjo, Winkler, Comisso, "Frozen-In Gravitational Fields", arXiv:2609.25240, for a timelike flow; Kerr–Schild potentials for the type-N hopfion in Sabharwal, Dalhuisen, JHEP 07 (2019) 004, arXiv:1904.06030 | not found; its relation to the Kerr–Schild potential is not checked |
| Helicities (2/3)π² and π² in that frame (MP23) | e_a ∧ T^a is the Nieh–Yan 3-form (arXiv:hep-th/9708138). Smołka, Jezierski, Class. Quantum Grav. 35, 245010 (2018), arXiv:1802.01467: a topological charge for linearised gravity, on hopfions. **Huang, Long, Zhao, Zhou, arXiv:2608.30154, §6.4: a Nieh–Yan-type helicity of the type-N hopfion in transverse-traceless gauge, (π²κ²/4)[1 + Re(1 − it)⁻⁶], which changes in time** | a helicity of the gravitational hopfion is prior art, in another gauge and with another answer; the frozen frame's value: not found; the normalisations are not compared |

## 2. The tube round a core line (`knotted`)

On the label sphere, α = cos η e^(iξ₁), β = sin η e^(iξ₂), Θ = pξ₁ + qξ₂. A magnetic line keeps H = cos^p η sin^q η cos Θ; the core lines are
cos²η = p/(p + q), Θ = 0 or π. In the torus frame the lines next to a core line turn by 2πν per circuit, **ν = √(2pq(p + q))**; on the torus
H = h the winding ι(h) rises from ν at the core to 2pq at the edge; and with u = 1 − h/h_max

**dι/du at the core line = ν · [3pq(p + q) − 2(p² + pq + q²)] / (12pq(p + q))**,

zero for the hopfion only. Three routes to it, the third the package's own:

| (p, q) | ν | closed form | from the profile integral (30 digits) | (ι − ν)/u at u = 0.004, by integrating the section | dι/dh |
|---|---|---|---|---|---|
| (1, 1) | 2.000000 | 0.000000 | 0.000000 | 0.0000 | 0 |
| (1, 2) | 3.464102 | 0.192450 | 0.192450 | 0.1929 | −0.500000 |
| (2, 3) | 7.745967 | 1.118862 | 1.118862 | 1.1216 | −6.018519 |
| (2, 5) | 11.832160 | 1.859339 | 1.859339 | 1.8639 | −15.092000 |
| (1, 8) | 12.000000 | 0.972222 | 0.972222 | 0.9744 | −4.671936 |

- `exact_checks()`: for symbolic p and q, the cross-section is a Hamiltonian system for 2 ln H; its linear part at the core line is
  [[0, −s], [2(p + q)s, 0]], s = √(pq)/(p + q); one circuit is a rotation by 2πν; the area of a level to second order gives the bracket.
- ν is a whole number, and the circuit matrix the identity, exactly when 2pq(p + q) is a square. For coprime p ≤ q ≤ 60 that is
  [(1, 1), (1, 8), (1, 49), (32, 49)] (MP09); for p = 1 it is Pell's equation, q = [1, 8, 49, 288]. Without "coprime", which the release plan's gate left
  out, there are 22 pairs, among them (p, p) for every square p, and (18, 32).
- For (1, 8) the matrix is the identity and the shear is 35/36: the neighbours close to first order and miss at third. H-G7, that such
  a core line answers a perturbation by its cube root, is a hypothesis and is **not tested**.

Helicity and energy, which are not in the lab's records. By hand, ∫A ∧ dA over the label sphere is 2π² pq·p! q!/(p + q)! and the energy
is (p + q) times it; the ratio is Kedia et al.'s. By quadrature:

| (p, q) | helicity, closed form | on the sphere | energy, closed form | in space at t = 0 | ratio |
|---|---|---|---|---|---|
| (1, 1) | 1.000000 π² | 1.000000 π² | 2.000000 π² | 2.000000 π² | 1/2 |
| (1, 2) | 1.333333 π² | 1.333333 π² | 4.000000 π² | 4.000000 π² | 1/3 |
| (2, 3) | 1.200000 π² | 1.200000 π² | 6.000000 π² | 6.000000 π² | 1/5 |
| (2, 5) | 0.952381 π² | 0.952381 π² | 6.666667 π² | 6.666667 π² | 1/7 |
| (1, 8) | 1.777778 π² | 1.777778 π² | 16.000000 π² | 16.000000 π² | 1/9 |

## 3. The squared hopfion's lines of force after t = 0 (`hopfion`)

With s = (1 − it(1 − α))^(1/2), σ² = 1 − it and P = β̄ (magnetic type) or iβ̄ (electric), the lines of force at the time t are the orbits of

**H = |P|² + |s² − σ²|²/t² on H = 1**, with the second invariant **I₂ = L² + 4m p_x² + (4m/t²) x² (r² − 2m)**

in the turned variables ζ = sσ̂*, p = Pσ̂*, m = √(1 + t²). So the lines lie on tori nested between two closed lines, the axial great circle
(−iσ̂² cos φ, σ̂* sin φ) and the cross line, and they are never chaotic. On the torus K = −I₂/(4m) the winding is W = T_ν/T_μ, two periods of
Y² = X(X − 1)Π(X); the whole picture depends on t only through ε = t/(2√(1 + t²)); to leading order W − 1 = (5/12)ε².

| t | ε | (5/12)ε² | W − 1, axial line | middle torus | beside the cross line |
|---|---|---|---|---|---|
| 0.1 | 0.0498 | 0.00103 | 0.00104 | 0.00103 | 0.00103 |
| 0.2 | 0.0981 | 0.00401 | 0.00411 | 0.00406 | 0.00401 |
| 0.4 | 0.1857 | 0.01437 | 0.01577 | 0.01506 | 0.01442 |
| 0.8 | 0.3123 | 0.04065 | 0.05483 | 0.04686 | 0.04109 |
| 1.5 | 0.4160 | 0.07212 | 0.13994 | 0.09521 | 0.07368 |
| 3.0 | 0.4743 | 0.09375 | 0.29630 | 0.13883 | 0.09659 |

What is checked, and how:

- `exact_checks()`: the change of label; {H, I₂} = 0 for symbolic t; the normal form at the 1 : 1 resonance; W = 1 + (5/12)ε² on the axial line.
- The reduction, in space, with exact gradients: a step along the field moves (s, P) along the Hamiltonian flow to 10⁻¹⁴ in angle (the lab,
  by differences: 2·10⁻⁸).
- **I₂ along lines of force integrated in space** on the field itself: the line through (0.6, −0.3, 0.5) keeps it to 3e-14 of its range
  over 20 circuits at t = 0.4 and to 1e-13 over 100 circuits at t = 0.8 (the lab, with scipy: 2.7·10⁻¹⁰). The package integrates with its
  own routine, an extrapolation of the midpoint rule, because scipy is an optional extra.
- The winding by a second route: the two periods measured along an orbit of the plane system equal the two integrals to 10⁻¹⁵. For that
  line at t = 0.4, W − 1 = 0.015660: **one turn in 63.9 circuits**.
- The two closed lines: tangent to the field, closed to 10⁻¹³, linked −1 (magnetic type) and +1 (electric) by Gauss's integral.

## 4. A tetrad frozen into the null flow, and its helicity (`tetrad`)

With a = x − iy, b = t − i − z, d = r² − (t − i)², β = 2a/d, γ = 1 − α, ∂a = (0, 1, −i, 0), ∂b = (1, 0, 0, −1):

**B_μ = Re[(i a² (∂_μ b) dγ − b² (∂_μ a) dβ)/d²],  𝓕_μ = dB_μ = Re[G_μ dα ∧ dβ],  G_μ = (i/d)(a ∂_μ b − b ∂_μ a).**

The tetrad e = dx + B is frozen into the light's null flow k (i_k 𝓕 = 0); the curvature of h = B + Bᵀ is the squared hopfion n ⊗ n/d⁵; its
Ricci tensor is zero. `exhibit(x)` checks this at an event with exact derivatives:

| arithmetic | i_k 𝓕 | E − Re W | B + Im W | Ricci |
|---|---|---|---|---|
| double precision, 60 events | 1.5e-15 | 1.1e-15 | 1.3e-15 | 2.0e-15 |
| 40 digits, 8 events | below 2e-40 | | | |
| Gaussian rationals, 10 rational events | exactly 0 | exactly 0 | exactly 0 | exactly 0 |

`certificate()` turns the last row into a statement about every event. Each of the 114 differences is a polynomial of total degree at
most 22 over a denominator with no real zero. By the Schwartz–Zippel lemma, exact zeros at ten events whose coordinates are drawn
from 4001 values each leave a chance below 3·10⁻²¹ that an identity which is not true has passed. That is a randomized certificate
(CERT), not a symbolic proof.

The helicity 3-form of a frame is e_a ∧ de^a; its part of second order, B_a ∧ dB^a, lives on the label sphere for this frame, so it does not
change in time:

| | H_g | H_em | ratio |
|---|---|---|---|
| on the label sphere | 6.5797362674 | 9.8696044011 | 0.6666666667 |
| in space, t = 0 | 6.5797362674 | 9.8696044173 | 0.6666666656 |
| in space, t = 0.8 | 6.5797362674 | 9.8696044173 | 0.6666666656 |
| closed form | (2/3)π² = 6.5797362674 | π² = 9.8696044011 | 2/3 |

- **The number belongs to the frame.** Huang et al. find a helicity of the same kind in transverse-traceless gauge that changes in time.
  The lab's fit found this frame to be the only one of its holomorphic class up to degree 4; frames outside that class are not excluded.
- 2/3 comes out of one integral over the label sphere. It has no mechanism, and no meaning is attached to it (AG-D5).

## 5. Misses and corrections, kept

From the lab's records:

- MP09's plan printed the trace 0.9868 for (2, 5), by hand; the formula gives 0.9871991, as its own result says.
- MP21's plan printed three hand values of dι/dh with slips (−15.2404 for −15.0920, −4.4244 for −4.6719, −6.018516 for −6.018519).
- MP22: the forecast that the line through (0.6, −0.3, 0.5) at t = 0.8 runs more than 10 from the origin missed; it reaches 4.7. One test
  counted the starting point as a crossing. The older records' "one turn in about 29 circuits" was wrong: 63.9.
- MP23: the forecast that the square of the light's potential is a potential for the squared hopfion missed. A prediction of 3/4 for the
  ratio of the helicities missed through a slip of algebra; it is 2/3.

Found in this build:

- MP09 says the line field on the label sphere equals its η′, ξ₁′, ξ₂′ "up to a positive factor". With B = Im F_pq the factor is negative
  everywhere. No result depends on the sense.
- MP23's number file holds the two helicities on the sphere with a minus sign, and its test took absolute values. The file's order of
  the sphere's coordinates, (η, ξ₁, ξ₂), is the orientation opposite to the one space gives the sphere. Here the sphere is oriented as
  space orients it, and sphere and slices agree in sign.

## 6. Open

- **T-G3**, a three-dimensional example with resistivity: not attempted.
- **H-G7**, the cube-root response of the core line of (1, 8): a hypothesis, not tested.
- The frozen tetrad against the Kerr–Schild potential of the same field, and its helicity against Huang et al.'s in one normalisation.
- The other gravitational (p, q) fields: no invariant is known for their later lines of force.
- A proof that the frame is unique among all frozen frames.
