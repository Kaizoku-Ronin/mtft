# mtft v0.9.1 — The m₀-EM Edge Predictor & the Restored Independence Legs

Patch release folding in the R1–R3 findings from Claude Fable's v0.9.0
review — each verified to the digit by the independent auditor before
acceptance (Addendum X; two of the three root causes are the auditor's
own, in the lift from the delivered suite to the module).

## Fixed — `edge_mass`: the m₀-anchored EM predictor (R1)

- v0.9.0's `predicted_asymptote` quietly carried a +f(M)/2 term — in
  neither the note nor Addendum U.4; the auditor's own unannounced
  "improvement" in the lift. Anchored at the real M it is one-sided:
  mean bias f(M)/2 over the fractional part, worst 15% at ε = 0.3.
  Removed. The 0.971-vs-0.9923 spread of the V.3 round was this term,
  not "boundary sensitivity" — V.3's mechanism attribution was wrong
  and is corrected in Addendum X.
- `predicted_asymptote` is now the bare U.4 law itself (mean-unbiased,
  ±f(M)/2 oscillation — the law `predicted_em` converges to).
- New primary predictor **`predicted_em`**: Euler–Maclaurin anchored
  at m₀, the actual first included level — deterministic, no
  fractional-part scatter. mass_plus_tail/pred_em = 1 within
  5×10⁻⁹ at ε = 0.3 (where the old predictor was off 15%), 2×10⁻¹²
  at ε = 0.2, 10⁻¹⁵ at ε = 0.1 — sub-ppm, as advertised in R1.
- New detail fields: `m0`, `theta_frac` (fractional part m₀−M),
  `tail_beyond_nmax` (EM estimate), `mass_plus_tail`,
  `ratio_em`, `ratio_em_corrected`. The nmax = 2×10⁶ truncation
  (7.5×10⁻⁵ of the mass at ε = 0.07) is now reported and corrected
  to sub-ppm — the instrument is ready for the spectral study's
  smaller ε.

## Fixed — the dead-gate class (R2)

v0.9.0's G4 compared `flow_phase` against a re-typed copy of its own
formula; `test_flow_phase_closed_form` did the same;
`test_edge_mass_convention_is_per_level`'s "brute force" was the
function's internals transcribed; and the `kms_check` docstring's
"both paths evaluated independently" covered two float orderings of
one termwise expression. The delivered suite's independent legs did
not survive the lift. Nothing computed was wrong — the labels claimed
an independence the checks no longer had.

- **Restored, stronger**: `flow_phase_matrix(p, t)` computes
  α_t(μ_p) = U μ_p U† by literal 400×400 matrix conjugation with
  U = diag(e^{itK}), K = −log ρ̂ built from the Gibbs weights — no
  spectral-formula input (even the delivered G4's matrix path read
  the closed spectrum; this one derives ΔE from the weights). G4 now
  compares the two genuine paths: worst 1.1×10⁻¹⁴.
- **`kms_check` carries four legs**: matrix (7.9×10⁻¹⁷), spectral
  with the rn-vs-rpn Gibbs form (2.8×10⁻¹⁷), restricted-spectral,
  and the cross-check — the cross evaluated on the matrix leg's own
  index set and normalization (cross-leg comparisons must be
  apples-to-apples; two normalization slips were caught and fixed
  in the fold-in).
- The two transcription tests are replaced: flow vs the matrix path,
  and the per-level edge mass against the exact Hurwitz-zeta tail
  Σ_{n≥a} (log n) n^{−3} = −ζ′(3, a).

## Fixed — `correlator` (R3)

- Now accepts complex t with Im t ≥ 0 — the KMS point t+i included
  (|e^{izΔE}| = e^{−Im z·ΔE} ≤ 1: the tail bound is unchanged).
  Raises ValueError on Im t < 0 (phases would amplify the tail).
- The test evaluates the KMS identity through the API instead of
  feeding t+i and discarding the result.

## Tests

- Suite 440 → 442 (two new edge-predictor tests; the three R2/R3
  replacements keep the count honest — no dead assertions).
