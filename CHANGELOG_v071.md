# v0.7.1 — Speiser–Hadamard Lab (July 18, 2026)

New tested section in `mtft/riemann.py` implementing the independently
verified ζ′ machinery behind the three-ensemble program's Dirichlet leg
(audit Addendum I; request: "anyone else can run these computations").

## Added

- **`ZETAPRIME_ZEROS`** — the certified census of all 19 nontrivial zeros
  of ζ′ with 0 < Im s ≤ 100 (17-digit, Newton-refined to |ζ′(ρ′)| < 1e-29;
  argument-principle certified complete; all Re > 1/2 — a numerical
  Speiser check consistent with RH to height 100).
- **`zetaprime_negative_zero(n)`** — the unique real zero of ζ′ in
  (−2n−2, −2n), solved *exactly* via the functional equation
  (H2(s) = log 2π + (π/2)cot(πs/2) − ψ(1−s) − (log ζ)′(1−s)); no ζ
  evaluations on the negative axis.
- **`zetaprime_refine(z0)`**, **`zetaprime_logcurvature(s)`**,
  **`zetaprime_zero_count_berndt(T)`** (N′(T) = (T/2π)log(T/4πe) + O(log T) —
  the ζ′ ensemble is sparser than ζ's by ≈ (T/2π)log 2).
- **`hadamard_zetaprime_check(s)`** — numerical evaluation of
  ∂²log(−ζ′(s)) = 2/(s−1)² − Σ_{ρ′}(s−ρ′)⁻² with both tails carried
  (negative-axis integral + high-γ Berndt density). Balances to
  |residual| ~ 1e-5 for s ∈ [3, 10]. Conditioning warning in docstring:
  the two tails nearly cancel at large s, and the LHS decays like (2/3)^s —
  use s ∈ [3, 10] for demonstrations.
- **`dirichlet_curvature(beta)`** — the decomposition lemma
  g_D(β) = ∂²log ζ(β) + ∂²log(−ζ′(β+1)) with component shares
  (ζ′ piece: 48.588864% at β = 3), plus **`von_mangoldt_curvature`**
  cross-check (∂²log ζ(β) = Σ Λ(n)(log n) n^{−β}).
- **`divisor_log_weights`**, **`weighted_theta`** — the Emergent weights
  w_n = Σ_{d|n} (log d)/d with W(s) = −ζ(s)ζ′(s+1).
- **`filtered_moment_identity(y, N)`** — the exact shift identity
  μ_N(y) = (1/4π²)[T″(y) − Re T″(y − i/N)] (verified to < 1e-12 relative;
  for N = 3 the (1 − cos(2πn/3)) factor is the SU(3) center projector up
  to 3/2).
- **`weighted_theta_cusp_fit`** — the modularity no-go: the double pole
  of W(s) at s = 0 forces Θ̃(y) = (−ζ′(2))/X − (1/4)ln²(1/X) + …, so the
  weighted theta has non-modular cusp asymptotics (fitted A = −0.2498
  vs predicted −1/4; C ≈ ζ′(0)).

## Verification provenance

All anchors recomputed independently (mtft audit, Addendum I): pole
coefficient 2 confirmed from raw ζ derivatives (no Hadamard input);
identity battery at s ∈ {3,…,30}; negative zeros bracket-certified to
n = 1000 and cross-checked against direct ζ′ sign changes; census
consistent with Berndt's count. 22 new regression tests
(`tests/test_zetaprime_hadamard.py`), **365/365 suite green**.

## Notes

- Version bump 0.7.0 → 0.7.1 (new API surface, fully backward compatible).
- Remaining v0.7.1-scope items from the earlier ledger are unchanged
  (trivial tower T(X), `predictions.lock`, pytest tier markers,
  paper-side y_c propagation).
