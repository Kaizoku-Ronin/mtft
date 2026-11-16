# The gravity ledger

Generated from `mtft.gravity.ledger` (`to_markdown()`); do not edit by hand. Entries are in the matter ledger's format. The pieces, the
pre-registered tests and the negatives are the lab's gravity map's own words (last changed 14 November 2026; SHA-1 `8730ed451d2fc3f6c85a5c281fa95b2c53ccab11`).
Nothing here is a claim that gravity is derived: N-G1 is the first negative.

12 entries (8 active, 1 with a caveat, 3 retracted and kept); 10 pre-registered tests (9 done, 1 pending); 16 negatives; 6 further recorded misses.

## Entries

| id | wave | class | status | statement | module | test | note |
|---|---|---|---|---|---|---|---|
| GR-1 | GR-01 | EXACT | active | Next to a core line of the (p, q) field the lines turn by 2 pi nu per circuit in the torus frame, nu = sqrt(2pq(p + q)): the cross-section is a Hamiltonian system for 2 ln H and its linear part at the core line has determinant 2pq/(p + q). The circuit matrix is the identity exactly when 2pq(p + q) is a square: among coprime p <= q <= 60, (1, 1), (1, 8), (1, 49), (32, 49). | `mtft.gravity.knotted` | `test_gravity_knotted.py` | the fields, their invariant Re(alpha^p beta^q) and its tori are Kedia, Bialynicki-Birula, Peralta-Salas and Irvine's (2013); nu was not found in the sources read on 2026-11-16, which is not the same as new; derived in MP09, unreviewed |
| GR-2 | GR-01 | EXACT | active | The shear at the core line, with u = 1 - h/h_max: d iota/du = nu [3pq(p + q) - 2(p^2 + pq + q^2)]/(12pq(p + q)), zero for the hopfion only; 35/36 for (1, 8). | `mtft.gravity.knotted` | `test_gravity_knotted.py` | three routes: the closed form (exact for symbolic p, q), the profile integral at 30 digits (4.9e-10 of nu), and the section integrated as a differential equation; derived in MP21, unreviewed |
| GR-3 | GR-01 | NUMERICAL | active | On the torus H = h the winding iota(h) rises from nu at the core line to 2pq at the edge. | `mtft.gravity.knotted` | `test_gravity_knotted.py` | MP09's integral, valid on every torus |
| GR-4 | GR-01 | EXACT | active | The helicity of the (p, q) field on the label sphere is 2 pi^2 pq p! q!/(p + q)! and its energy is (p + q) times that. | `mtft.gravity.knotted` | `test_gravity_knotted.py` | by hand and by quadrature to 1e-14; the ratio 1/(p + q) is Kedia et al.'s published value, the one outside anchor of the module's normalisation; not in the lab's records |
| GR-5 | GR-02 | EXACT | active | At a time t > 0 the lines of force of the squared hopfion are the orbits of H = \|P\|^2 + \|s^2 - sigma^2\|^2/t^2 on H = 1, with the second invariant I_2 = L^2 + 4m p_x^2 + (4m/t^2) x^2 (r^2 - 2m): {H, I_2} = 0 for symbolic t. They are never chaotic. | `mtft.gravity.hopfion` | `test_gravity_hopfion.py` | that the lines of force are the light's lines turned by arg(d)/2 is Swearngin et al.'s (2013); the potential reached is Garnier's integrable quartic; the reduction itself was not found in the sources read; derived in MP22, unreviewed |
| GR-6 | GR-02 | NUMERICAL | active | In space, on the field itself: the reduction holds pointwise to 1e-14 in angle; I_2 is kept to 1e-13 of its range along lines of force over 100 circuits; the axial great circle and the cross line are closed lines of force, linked -1 (magnetic type) and +1 (electric). | `mtft.gravity.hopfion` | `test_gravity_hopfion.py` | integrated by the package's own routine; the lab's runs with scipy kept I_2 to 2.7e-10 |
| GR-7 | GR-02 | NUMERICAL | active | The winding on a torus is W = T_nu/T_mu, two periods of Y^2 = X(X - 1)Pi(X) depending on t only through eps = t/(2 sqrt(1 + t^2)); to leading order W - 1 = (5/12) eps^2 on every torus (exact on the axial line to that order). For the line through (0.6, -0.3, 0.5) at t = 0.4, W - 1 = 0.015660. | `mtft.gravity.hopfion` | `test_gravity_hopfion.py` | MP22's table to 1e-16 of the lab's file; measured along an orbit to 1e-15 of the integrals |
| GR-8 | GR-03 | NUMERICAL | active | The tetrad e = dx + B, B_mu = Re[(i a^2 (d_mu b) d gamma - b^2 (d_mu a) d beta)/d^2], is frozen into the null flow of the squared hopfion; the curvature of h = B + B^T is n (x) n/d^5 and its Ricci tensor is zero. | `mtft.gravity.tetrad` | `test_gravity_tetrad.py` | exactly zero differences over the Gaussian rationals at the events tried; a randomized certificate (Schwartz-Zippel, degree at most 22) bounds a false pass by 3e-21; not a symbolic proof; the degree count is unreviewed. Frozen-in gravitational fields for a timelike flow are Asenjo, Winkler and Comisso's (2026); this null-flow frame was not found in the sources read; derived in MP23, unreviewed |
| GR-9 | GR-03 | NUMERICAL | caveat | In that frame the helicity is H_g = (2/3) pi^2, beside H_em = pi^2 for the light field, on the label sphere and on slices of space, constant in time. | `mtft.gravity.tetrad` | `test_gravity_tetrad.py` | the number belongs to the frame (N-G5): in transverse-traceless gauge Huang, Long, Zhao and Zhou (arXiv:2608.30154) find a helicity of the same kind that changes in time, and the two normalisations have not been compared; the frame is unique only within the lab's holomorphic class up to degree 4; 2/3 has no mechanism (AG-D5) |
| GR-10 | GR-02 | NUMERICAL | retracted | At t = 0.4 the line through (0.6, -0.3, 0.5) makes one turn in about 29 circuits (MP10). | `mtft.gravity.hopfion` | `test_gravity_hopfion.py` | a miscount of returns; one turn in 63.9 circuits (MP22; N-G15) |
| GR-11 | GR-02 | ESTIMATE | retracted | At t = 0.8 the line of force through (0.6, -0.3, 0.5) is lost (MP11's reading); it reaches more than 10 from the origin (MP22's forecast, given 0.6). | `mtft.gravity.hopfion` | `test_gravity_hopfion.py` | it was not lost: the search plane missed it for an arc of 69.8; it keeps I_2 on a torus at 0.43 of the range, returns, and reaches 4.7 (MP22; N-G15) |
| GR-12 | GR-03 | ESTIMATE | retracted | H_g/H_em = 3/4 (a hand value, registered in MP23 before the integral was computed). | `mtft.gravity.tetrad` | `test_gravity_tetrad.py` | a slip of algebra; the ratio is 2/3 (N-G16) |

## Waves

- **GR-01** (2026-11-16): the tube round a core line (lab records MP09, MP21). Report: `docs/gravity/KNOTTED_FIELDS.md`.
- **GR-02** (2026-11-16): the squared hopfion's lines of force after t = 0 (lab record MP22). Report: `docs/gravity/KNOTTED_FIELDS.md`.
- **GR-03** (2026-11-16): a tetrad frozen into the null flow, and its helicity (lab record MP23). Report: `docs/gravity/KNOTTED_FIELDS.md`.

## The map's pieces

- **I** The coupling: Newton's constant: open
- **II** The arena: where gravity sits on X₀(143): open
- **III** Gravity as a square (the coin): partial
- **IV** The frozen field: partial
- **V** Breaking and reconnecting: open
- **VI** The field equations: open
- **VII** The vacuum and the expansion: open
- **VIII** Dark matter: open
- **IX** Gravitational waves: open
- **X** Black holes: open
- **XI** Is gravity quantum? The low-mass frontier: open

## Pre-registered tests (the map's words)

| test | registered | criterion | outcome |
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

## Honest negatives (the map's words)

| id | pieces | what was closed | record |
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
| N-G16 | III | The coin does not work at the level of potentials: the square d·σ ⊗ σ of the light's potential over the scalar is not a potential for the squared hopfion; its own Ricci tensor is 5.9 times \|W\| (forecast at 0.75, missed). And the hand value H_g/H_em = 3/4, registered before the integral was computed, was wrong: it is 2/3. | MP23 |

## Further recorded misses

Among the negatives, N-G7, N-G10, N-G13, N-G14, N-G15, N-G16 are misses of a forecast or corrections of an earlier statement. Beside them:

| id | record | kind | what |
|---|---|---|---|
| M-G1 | MP09 | hand value | the plan printed the trace 0.9868 for (2, 5); the formula gives 0.9871991, as the record's own result has it |
| M-G2 | MP21 | hand value | three hand values of d iota/dh in the plan were wrong (-15.2404 for -15.0920, -4.4244 for -4.6719, -6.018516 for -6.018519); the closed form was right |
| M-G3 | MP22 | fault of a test | a count of crossings included the starting point, so a check printed MISS on a line that closes; corrected in the record |
| M-G4 | MP09 | found in the 1.1.0 build | the line field on the label sphere is said to hold 'up to a positive factor'; with B = Im F_pq the factor is negative; no result depends on the sense |
| M-G5 | MP23 | found in the 1.1.0 build | the number file holds the two helicities on the label sphere with a minus sign and the test took absolute values: the order (eta, xi_1, xi_2) is the orientation opposite to space's; the package uses (xi_1, eta, xi_2) |
| M-G6 | the release plan of 1.1.0 | found in the 1.1.0 build | the plan's gate for the whole-number pairs left out 'coprime'; without it there are 22 pairs for p <= q <= 60, not four |

## Open

- **T-G3**: a three-dimensional example with resistivity: not attempted
- **H-G7**: the core line of (1, 8) answers a perturbation by its cube root: a hypothesis of MP21, not tested
- **the Kerr-Schild comparison**: whether the frozen tetrad's metric is the Kerr-Schild potential of the same field in another gauge: not checked
- **the helicity comparison**: the frozen frame's (2/3) pi^2 against Huang et al.'s transverse-traceless value in one normalisation: not done
- **the other gravitational (p, q) fields**: no invariant is known for their later lines of force
- **the frame's uniqueness**: a proof that the frozen frame is unique among all frozen frames: not attempted
