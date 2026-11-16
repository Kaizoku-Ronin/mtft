# GB-01 — H24 against the lattice: is the glueball spectrum the Maass spectrum of X₀(143)? (2026-11-16)

Module `mtft.maass.glueballs`; tests `tests/test_maass_glueballs.py`, `tests/test_gb01_glueball_comparison.py`. Ported from the Modular
Plane lab's records MP15, MP18, MP19 and MP20 (12 to 14 November 2026), which are kept under `studies/maass_2026nov/lab/`. The spectrum it
uses is MA-01's frozen table (`mtft.maass.spectrum`). Class: PHENO. Everything here is a comparison under a stated dictionary.

**The answer.** H24 as written fails, in both of its readings. No claim about Yang–Mills theory follows from the curve. Two forecasts of
one fitted description are still open; they are sealed and have not been looked up.

## 1. The hypothesis, the lattice and the dictionary

H24 (R&D notes, 16 October 2026): if MTFT's configuration space is X₀(143) with the hyperbolic Laplacian as Hamiltonian, then the gap
is arithmetic, and the glueball mass ratios would be ratios of the Maass spectrum of Γ₀(143).

The lattice: Athenodorou and Teper, JHEP 11 (2020) 172 (arXiv:2007.06422), SU(3), continuum limit, in units of √σ. Six values are used:
0⁺⁺ 3.405(21); 0⁺⁺* 5.855(41); 0⁻⁺ 5.276(45); 0⁻⁺* 7.29(13); A₁⁺⁻ 9.32(28); A₁⁻⁻ 9.78(46). The four with J = 0 identified are from the paper's Table 17; A₁⁺⁻ and A₁⁻⁻ are the lightest states of their channels, from its
Table 14, with spins not identified. **They are the lab's transcription of 12 November.** They were not re-read against the paper when
this module was written, because the same tables would show the excited C = −1 states that section 5's sealed forecasts are about.

The dictionary, fixed before the spectrum was computed (MP15), with one later variant (MP18):

| | |
|---|---|
| D1 | energy: λ = ¼ + R²; the unit is free, so only ratios are tested |
| D2 | spin: weight 0 is J = 0 only; the 2⁺⁺, the second-lightest glueball, has no partner |
| D3 | parity: the mirror τ ↦ −τ̄ |
| D4 | charge conjugation: the Fricke sign ε₁₁ε₁₃ |
| D4′ | charge conjugation (MP18): C = −1 for the η-twisted forms |
| D5 | states: reading a, all of L² with the 0⁺⁺ at the edge ¼; reading b, the cusp forms only |

## 2. Reading b: all five predictions miss (GB-1)

With the cusp forms only and C the Fricke sign, the channels are u.e.mm.1 (0⁺⁺), u.e.pp.1 (0⁺⁺*), u.o.mm.1 (0⁻⁺)
and u.o.pp.1 (0⁻⁺*). `reading_b()` recomputes the five registered predictions from the frozen table:

| | Prediction | Found |
|---|---|---|
| N0 | the lowest cusp form is even with Fricke sign + | it is odd |
| N1 | λ(0⁻⁺)/λ(0⁺⁺) in [1.472, 1.627] | 0.870 |
| N2 | λ(0⁺⁺*)/λ(0⁺⁺) in [1.634, 1.806] | 2.236 |
| N3 | λ(0⁻⁺*)/λ(0⁺⁺) in [2.020, 2.262] | 1.247 |
| N4 | no form with Fricke sign − below 2.60 λ(0⁺⁺) | one at 2.152 |

Reading b is closed for this dictionary.

## 3. Reading a and the η-twist (GB-2, GB-3)

Reading a puts the 0⁺⁺ at the edge of the continuum, so every energy is λ/¼ = 1 + 4R². With D4′ for the C = −1 states
(`comparison()`, in the form of `mtft.falsify`: PASS within 2σ, TENSION to 3σ, FAIL beyond):

| Lattice state | Curve state | Lattice ratio | Curve ratio | σ | Needed shift of λ | |
|---|---|---|---|---|---|---|
| 0⁻⁺ | odd (−,−), R = 0.37385 | 1.549(16) | 1.559 | +0.6 | −0.0024(41) | PASS |
| 0⁺⁺* | even (−,−), R = 0.44498 | 1.720(16) | 1.792 | +4.5 | −0.0181(40) | FAIL |
| 0⁻⁺* | odd (+,+), R = 0.55555 | 2.141(40) | 2.235 | +2.3 | −0.0234(101) | TENSION |
| A₁⁺⁻ | twisted even (+,+), R = 0.48009 | 2.737(84) | 1.922 | −9.7 | +0.2038(210) | FAIL |
| A₁⁻⁻ | twisted odd (−,−), R = 0.46274 | 2.872(136) | 1.857 | −7.5 | +0.2539(341) | FAIL |

- **Reading a fails on structure, whatever the numbers.** Its 0⁺⁺ is the edge of a continuum: no isolated state, and an effective mass
  that falls like Δ + 3/(2t) with no plateau, where the paper reports plateaux from t = a or 2a. That was derived before the run (S1) and
  was tested against the paper's description, not its correlators. With C a symmetry of the curve, a C = −1 threshold state also sits at
  the gap, since all four Eisenstein series are even.
- **One registered prediction passed.** S2: the lowest odd form with Fricke sign + in R ∈ [0.344, 0.396]. It is at 0.37385:
  1.559 against the lattice's 1.549(16), 0.6σ. The forecast had been that it would miss (0.9). The lab's control found no
  such form below R = 0.8 at levels 133 and 145, so the match is specific to 143; it is recorded as a coincidence without a mechanism
  (AG-D5), attached to a reading that fails on structure. The second ratio, 2.235 against 2.141(40) (2.3σ), was found after looking.
- **The 0⁺⁺\*** is 4.5σ off, and the trace formula excludes an even form in R ∈ [0.410, 0.438], where the lattice would need one.
- **The η-twist is a registered miss.** D4′ clears the edge of C = −1 states, as designed, and then puts four twisted forms below the
  registered bound R = 0.6325. The lightest C = −1 states would sit at 1.857 and 1.922 times the edge, 7.5σ and 9.7σ below the lattice.
  (MP18 first printed 9.6σ, from a slip in one lattice error, 0.085 for 0.084; its own list of faults corrects it.)

## 4. Could it be a hidden self-interaction? (GB-4)

A term −Δu + g|u|²u shifts each λ, to first order, by κ·I with I = area·∫u⁴/(∫u²)² the participation ratio. The seven ratios are
1.8905, 2.3288, 2.2821, 2.1920, 1.9035, 1.8884, 2.4722; the package's own quadrature over the Ford domain reproduces the lab's to 10⁻⁸ (a test, 12 s).

| Version | Result |
|---|---|
| one coupling for the three C = +1 states (`self_interaction_fit`) | κ₊ = −0.0059, χ² = 10.2 on 2, p = 0.006: fails, as forecast (0.75). With no interaction χ² = 24.4 on 3 |
| one coupling for all five (`one_coupling_fit`) | χ² = 170 on 4 with independent errors, the lab's figure; 178 on 4 with the shared error of the 0⁺⁺. The needed shifts are negative for C = +1 and positive for C = −1 |
| a coupling for C = −1 fixed by A₁⁺⁻ (`minus_sector_fit`) | the lowest odd twisted state lands at 2.564 against 2.872(136): −1.99σ with the coupling's error carried, −2.26σ against the lattice's alone. Passed by the letter of a plan that did not say which; as a fit χ² = 3.97 on 1, p = 0.046. No support |

The detail about the 170: the lab computed it after its run and kept no code. Rebuilding it showed that it used independent errors,
while the three-state fit used the covariance. Both versions are in the function; the verdict is the same.

## 5. A constant twist energy, and two sealed forecasts (GB-5)

One number on every η-twisted form, μ = 0.2176(179) in λ, fits the two C = −1 ground states: χ² = 1.57 on 1, p = 0.21
(`twist_energy`). It is one fitted number for two measurements with errors of 3 % and 5 %, so it could not have failed, and at
0.87 of the edge it is another operator, not a correction to this one. It is recorded as a description, not as evidence.

It does forecast. With μ fixed, the four lowest twisted states sit at 2.73, 2.79, 2.87, 3.24, each ± 0.07, in units
of the 0⁺⁺. The first two are the fitted ground states. The other two are open:

| | Forecast | Against |
|---|---|---|
| second odd C = −1 state | 2.8658 ± 0.0715 (9.76 in √σ) | the first excited state of the lattice's A₁⁻⁻ channel |
| second even C = −1 state | 3.2399 ± 0.0715 (11.03 in √σ) | the first excited state of its A₁⁺⁻ channel |

The lab registered the two numbers at 08:21 New York time on 14 November 2026, before any lookup, and gave no criterion. The criterion
was fixed on 16 November at 06:54:50, still before any lookup, and the text, the numbers and the rule's code were hash-sealed as GB-S1
(`SEALS`; sha256 `e75106afbfc3e3bd92c29ebcfeea81db48a0c67ab715c50fdc29dad1735500e7`): σ = |forecast − lattice ratio|/√(0.0715² + e²); PASS to 2σ,
TENSION to 3σ, FAIL beyond; a channel with no excited state listed has no verdict; a fail in either channel closes the description.
`evaluate_open(excited_odd, excited_even)` runs the sealed rule and refuses if the seal no longer verifies.

**The lookup has not been made.** It is experiment E5 of the release plan. The six lattice values of section 1 should be checked against
the paper at the same time.

**Added later on 16 November 2026.** The lookup was made that day: `GB02_E5_SEALED_FORECASTS.md`. The odd forecast passed at 1.4σ with little
power; the even one has no lattice value; the six values are right; and the source assigns the A₁⁺⁻ ground state to spin 4, which this report
did not know. Everything above stands as it was written before the lookup.

## 6. Every registered forecast and its outcome

`register()` returns this table; the outcomes that can be recomputed from the frozen data are checked against it in the tests.

| id | registered (New York) | whose | statement | forecast | outcome | result |
|---|---|---|---|---|---|---|
| MP15-N0 | 2026-11-12 13:31 | H24, reading b | the lowest cusp form of Gamma_0(143) is even with Fricke sign +: the 0++ is the mass gap |  | missed | the lowest cusp form is odd |
| MP15-N1 | 2026-11-12 13:31 | H24, reading b | lambda(0-+)/lambda(0++) in [1.472, 1.627] |  | missed | 0.870 |
| MP15-N2 | 2026-11-12 13:31 | H24, reading b | lambda(0++*)/lambda(0++) in [1.634, 1.806] |  | missed | 2.236 |
| MP15-N3 | 2026-11-12 13:31 | H24, reading b | lambda(0-+*)/lambda(0++) in [2.020, 2.262] |  | missed | 1.247 |
| MP15-N4 | 2026-11-12 13:31 | H24, reading b | no cusp form with Fricke sign - below 2.60 lambda(0++) |  | missed | one at 2.152 |
| MP15-S1 | 2026-11-12 13:31 | derived | under reading a the 0++ is the edge of a continuum: no isolated state, no plateau |  | stands | reading a fails on structure; tested against the paper's description of its plateaux, not its correlators |
| MP15-S2 | 2026-11-12 13:31 | H24, reading a | the lowest odd cusp form with Fricke sign + lies in R in [0.344, 0.396] | 0.1 | passed | R = 0.37385: 1.559 against the lattice's 1.549(16), 0.6 sigma; the forecast that it would miss (0.9) was wrong |
| MP15-F1 | 2026-11-12 13:31 | forecast | in every one of the eight sectors the lowest cusp form has R < 1.5 | 0.7 | held | 0.37 to 1.41 |
| MP15-F2 | 2026-11-12 13:31 | forecast | reading b misses at least one of N0 to N4 | 0.9 | held | all five |
| MP15-F3 | 2026-11-12 13:31 | forecast | N4 fails | 0.85 | held |  |
| MP15-F4 | 2026-11-12 13:31 | forecast | S2 fails: the window is narrow | 0.9 | missed | S2 passed |
| MP15-F5 | 2026-11-12 13:31 | forecast | N0 fails: the lowest cusp form overall is odd | 0.55 | held |  |
| MP15-control | 2026-11-12 13:58 | control | at levels 133 and 145 no odd cusp form with Fricke sign + lies below R = 0.8 |  | held | the S2 match is specific to 143: a coincidence without a mechanism; a recorded result of the lab's scripts, not recomputed here |
| MP18-P1 | 2026-11-13 08:52 | D4', the test | no eta-twisted form has R < 0.6325 | 0.1 | missed | four do; the lightest C = -1 states at 1.857 and 1.922 times the edge, 7.5 and 9.7 sigma below the lattice; the forecast of a miss (0.9) held |
| MP18-P2 | 2026-11-13 08:52 | two routes | the twisted trace formula closes on the first full list | 0.7 | missed, then held | three forms were missing from the first list of 39; with all 42 it closes |
| MP18-F1 | 2026-11-13 08:52 | forecast | the lowest twisted form has R < 0.5 | 0.75 | held | 0.4627 |
| MP19-P4 | 2026-11-13 13:39 | forecast | no even form with Fricke sign + lies in R in [0.410, 0.438], where the lattice's 1.720(16) would need one | 0.97 | held | reading a's 0++* stays at 1.792, 4.5 sigma from the lattice |
| MP20-V13 | 2026-11-14 07:59 | check | the quadrature gives the area to 1e-4 and orthogonality to 1e-5 |  | missed, then held | the first rule gave orthogonality 1.4e-4; on panels between the domain's breakpoints 2e-10 |
| MP20-SIU | 2026-11-14 07:59 | decided by sign | one self-coupling for every state |  | excluded | the needed shifts are negative for C = +1 and positive for C = -1; chi2 = 170 on 4 |
| MP20-P1 | 2026-11-14 07:59 | test | one coupling for the three C = +1 states passes at p > 0.05 | 0.25 | failed | chi2 = 10.2 on 2, p = 0.006; the forecast of a fail (0.75) held |
| MP20-P2 | 2026-11-14 07:59 | test | a coupling for C = -1 fixed by A1+- puts the lowest odd twisted state within 2 sigma of A1-- | 0.65 | passed by the letter | -1.99 sigma with the coupling's error, -2.26 with the lattice's alone; no support |
| MP20-F1 | 2026-11-14 07:59 | forecast | all seven participation ratios lie between 1.5 and 4.5 | 0.8 | held | 1.89 to 2.47 |
| GB-S1 | 2026-11-14 08:21 | the twist energy | a second odd C = -1 state at 2.87(7) and a second even one at 3.24(7), in units of the 0++ |  | open | sealed on 2026-11-16 before any lookup |

The "forecast" column is the probability given beforehand that the statement would hold, where one was given.

## 7. A new dictionary

A dictionary not tried here is tested forecast first. `register_forecast(name, statement, predictions, registered)` seals the predicted
ratios with the time they were made; `compare_forecast(record)` takes only such a record, unchanged, and answers in the form of
`mtft.falsify`. Nothing in the package compares an unsealed number with the lattice table.

## 8. What stands

- S1: in reading a the 0⁺⁺ is the edge of a continuum. No twist and no interaction of the kinds tried changes that.
- The coincidence at 1.559 is a coincidence (AG-D5).
- The puzzle map records the whole as negative N24.
