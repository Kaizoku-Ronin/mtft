# scripts/pari — PARI/GP provenance scripts

The computations behind the X₀(143) period/cusp analysis (Papers 32–33
era).  The `.txt` files are the verified session logs these scripts
produced (banner lines carrying local filesystem paths were trimmed;
the mathematical content is untouched).

| File | Computes |
|---|---|
| `mtft_period_matrix_v4.gp` | 50-digit modular-symbol periods on X₀(143): Petersson ⟨f₁,f₁⟩, mass-cycle {∞, 2/77} periods, trace/L1/L2 depth ratios vs the 1.52937 target, the analytic Q2/Q4 split, and the exactly-determined 2-parameter stiffness-period mass fit |
| `mtft_extended_analysis.gp` | Γ₀(143) cusp classes and widths, coefficient-field factorization (P irreducible ⇒ no algebraic Q2/Q4 split), the 544-cycle scan c ≤ 50, and 143a1 j-invariant / rank / torsion / L-values |

Run with:

```bash
gp -q mtft_period_matrix_v4.gp
gp -q mtft_extended_analysis.gp
```

Cross-checks against `mtft.x0_143` (performed at v0.7.0 integration):
field minimal polynomials = `FIELD_POLY_F2/F3`, H₃ roots =
`hecke_polynomial_f3_T2()` roots (max diff 5e−15), f₁ q-expansion =
`ORBIT_TRACE_F1`, cusp widths 143+13+11+1 = 168 = `INDEX`.  All match.

Re-validated at v0.7.0 on stock PARI/GP 2.15.4 (Linux):
`mtft_extended_analysis.gp` runs to completion with the same known
cosmetic warnings; the period script needs several minutes and a
≥ 2 GB PARI stack.  As part of the same pass, `HECKE_TRACES` (the
package's 200 trace-form entries) were re-verified exactly against
`mftraceform([143,2],0)`.

Known cosmetic issues in the logs (kept for honesty, flagged for a
cleanup pass): a `psi` variable name collides with a GP builtin, one
`ellj` call passes the wrong argument type, and `myphase` fails on
t_POL periods — none affect the quantities cross-checked above.  Also
flagged: PARI's `ellL1(E,1) = 0.94570` vs the package's
`L_VALUES["L'(f1, 1)"] = 0.791` (normalization to reconcile in v0.7.1),
and the log's "Cuspidal dim = 16" label actually reports dim M₂
(13 cusp + 3 Eisenstein).
