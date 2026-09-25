# mtft v0.9.0 — The Marked Primon Gas & the A.7 Envelope Patch

One release carrying both the planned v0.8.1 envelope patch and the new
marked-gas module (semver: a new module is a minor bump), per the author's
call. Audit trail: Addenda T–V; Addendum W records this release.

## Added — `mtft.marked_gas` (Tier 5d)

The marked primon gas of the July 2026 note (v0.1.1), shipped as a
certified module. The Gibbs state
ρ̂_n = (log n) n^{-(β+1)} / (-ζ′(β+1)) on ℓ²(ℕ≥2), H = log Q, spectrum
E_n = (β+1) log n - log log n, prime-shift isometries μ_p|n⟩ = |pn⟩.

- **`zD_certified_interval(beta, N)`** — the two-sided hyperbola-method
  interval for the Dirichlet-ensemble closed form Z_D(β) = -ζ(β)ζ′(β+1).
  Width law (log N)²/(2N²): 9.6e-11 at N = 1e6, 1.4e-12 at N = 1e7.
  Kahan accumulation plus an explicit float64 slack budget
  (4 · 2.3e-16 · log₂N · (|S|+2)); the reorder identity gate (G3b) guards
  the H_M off-by-one bug class the auditor caught in review.
- **`kms_check(beta, p, t)`** — the modular flow α_t = Ad e^{itK}
  satisfies KMS at **t+i** termwise: ρ̂_n e^{-ΔE_n} = ρ̂_{pn}. Ships the
  wrong-sign control (t-i fails by 4.574, as it must) and the BC twist
  (1 + log p/log n)^{-it} → 1 in the UV.
- **`cold_gas_report(N, lo)`** — Ψ(q) = Π (1-qⁿ)^{-(log n)/n²}. Karamata
  closed form A_n ~ B n^α / Γ(α+1) with α = -ζ′(2) and
  B = e^{-ζ″(2)} (the γ-cancellation between the Γ and ζ(s+1) poles).
  Certified amplitude 0.14040 against the closed form 0.14027492...
  (0.09%) at N = 1e5, converging from above at O((log n)²/n). Every
  estimator in the report is documented (binned slope, two-point slope,
  mean-of-logs and endpoint amplitudes, free least squares) per the
  audit redline.
- **`edge_mass(p, beta, eps)`** — the auditor's edge-softness law at the
  prime resonance (β+1) log p:
  mass(gap < ε) = [M^{-β} (log M/β + 1/β²) + O(log M · M^{-β-1})] /
  (-ζ′(β+1)), M = exp(log p / (e^ε - 1)). Per-level evaluation
  convention pinned (note v0.1.1) — the boundary-rounding artifact that
  put the two engines on opposite sides of 1 is closed.
- **`gates(quick=True)`** — the §7 green suite as a callable: 10
  recorded gates covering G1–G7 (the delivery's 13-line log splits
  G6/G7 further; G0's environment anchors are the pytest suite's job).
  Quick tier ~3 s (interval at 1e6, gas at 8 000, 3 refined zeros);
  full tier ~16 s reproduces the delivered numbers (interval width
  1.41e-12 at 1e7, cold slope 0.936235, amplitude 0.14040, ζ′-zero
  refinement residual 2.2e-56 at dps 50).
- Legend entries: `marked_gas`, `kms_flow`, `cold_gas_amplitude`,
  `spectral_edge_soft` — all Tier 5d, upstream `Z_D_closed_form`.
  `mtft.legend.trace("kms_flow")` walks the chain to the integers.

```python
>>> import mtft
>>> rep = mtft.cold_gas_report(100_000)
>>> rep["amplitude_mean_of_logs"]   # 0.14040 vs closed form 0.14027492...
>>> mtft.marked_gas_gates()["all_green"]
True
```

## Fixed — `envelope_slope` now follows the A.7 discipline (audit S.4 / T-E4)

- The old default was a hard-coded 91-point grid over 5.2 decades
  (~17.5 samples/decade against γ₁ = 5.18 periods/decade) with no
  minimum-bin guard: every bin underpopulated, the slope dominated by
  terminal-bin leverage — exactly the failure mode Tier 11's
  `estimator_standards` was written to prevent.
- Now: default density from `recommended_samples_per_decade(γ₁,
  per_period=10)` = **52 samples/decade** (272 points over the default
  range); binning and the terminal-bin guard via `binned_log_slope`
  (bins with < `min_bin` points dropped; `ValueError` if fewer than 3
  usable bins remain — the legacy 91-point grid now correctly discards
  all 11 of its bins).
- Validated: on-line slope 0.0000; off-line quadruplet
  β₀ = 0.6/0.75/0.9 → -0.069/-0.237/-0.392 against the 1/2 - β₀
  anchors; densification-stable (52 vs 208 samples/decade agree to
  ≤ 6e-4).
- Regression coverage: `tests/test_envelope_a7.py` (7 tests), including
  the γ₃ × 6-decade stride-resonance trap (54.994 cycles) from the A.7
  sessions.

## Changed

- `bakry_emery_curvature` / `tower_rigidity` docstrings: the
  "κ(y) ≥ 0 ⟺ RH" framing is now explicitly marked **legacy** — false
  unconditionally (κ^Λ(y) < 0; audit B4). The live criterion is the
  boundedness diagnostic of `corrected_rh_diagnostic()` (Th 1). No
  code-path change.

## Tests

- Suite: 401 → 440 (32 new in `tests/test_marked_gas.py`, 7 new in
  `tests/test_envelope_a7.py`).
