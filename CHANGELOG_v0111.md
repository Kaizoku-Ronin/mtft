# mtft v0.11.1 — documentation-correction wave (2026-08-04)

Scope: **zero computational changes.** Four labeling/documentation
corrections, one runtime guard, one new regression test, three-way
version bump. All engine numbers, certificates, and the wheel's import
surface are unchanged from v0.11.0.

## Closed — PET_F1 certification (supersedes the Add. BQ leg-3 open item)

The "+4.7%" gap between the corpus <f1,f1> and the modular-degree route
was a normalization mismatch, not a discrepancy: PARI's `mfpetersson`
normalizes by the index [PSL2(Z):Gamma_0(143)] = 168, while the BQ
comparison used the hyperbolic volume 56pi. The ratio 56pi/168 = pi/3
= 1.047198 is the entire effect. With the index normalization:

    deg(143.a1) = 4                                   (LMFDB)
    covol = Om_re * Im(om2)
          = 1.9699231645720704 * 2.0157723238708093
          = 3.9709165952963810                        (quadrature, 40-digit)
    <f1,f1> = 4 * covol / (4 pi^2 * 168) = 0.00239486886655019

matching the corpus 0.002394868866550 to 7.8e-14 — all published digits.
**PET_F1 upgraded Cert-corpus -> Cert** (three independent routes:
mfpetersson / Rankin-Selberg residue (BQ) / modular degree + period
lattice). Retro-calibration of BQ leg 3: the f2/f3 "2-3% uniform"
residuals are estimator bias (consistent with the +1.2% calibration);
PET_F2/F3 remain Cert-corpus, per-embedding order DIAGNOSTIC.

**Explicit retraction (corpus protocol):** the v0.11.0 statement
"<f,f>_raw = PET * V_N, V_N = 56pi" is superseded by
"<f,f>_raw = PET * 168". The v02 source comment is corrected and an
annotation is appended to CHANGELOG_v0110.md (the record is preserved,
not rewritten).

## Corrected — photon-ledger label (owned error: Claude; caught: Kimi BQ)

The v0.11.0 changelog recorded this relabel but the source never
received it. Now landed in `studies/x0143_particle_box_v03.py`
(docstring, inline comment, console line, figure title): the identity
f + M^T w = 0 is forced by the vanishing column sums of the master
generator — a construction certificate (catches indexing bugs at
1e-16), not an E2 pair. The honest numerical figure is the trajectory
quadrature residual (~7e-3, time-grid limited). The E2 label on v01's
charpoly-vs-oracle certificate is genuinely two independent routes and
is unchanged.

## Corrected — stale corpus flag in the v02 study

`studies/x0143_particle_box_v02.py` still printed "flag for corpus
correction" for `CURVE_143A1.j_invariant`, which v0.11.0 already fixed —
so the shipped study contradicted the shipped package. Replaced with a
runtime oracle guard: the study parses the oracle's j string and asserts
agreement with the exact a-invariant value (-262144/1859) and with the
lattice j, failing loudly on any regression. The Paper 32 tau flag
(0.039+0.980i / "nearly square") remains open on the paper side —
Delta = -1859 < 0 forces tau = 1/2 + it.

## Tests

- `test_pet_f1_modular_degree_certification` — locks the index
  normalization and the three-route agreement (LMFDB degree + lattice
  quadrature vs PET_F1; rel. tol 1e-8, scipy-only, ~0.1 s).

## Open items (carried)

- Capture ceiling = 0.500000 exactly: gamma0-independence is structural
  for uniform rate scaling and is NOT evidence of anything; the AG-D5
  test is a sink-exchange automorphism check plus V0 variation. Staged.
- Paper 32 tau / "nearly square" erratum (paper-side).
- f3 Petersson norm <-> embedding pairing (PARI mfpetersson diagonal
  matched to embedding order) — the gate on the orbit-Zeno interval.

Version: 0.11.1 three-way (pyproject / __init__.py / CITATION.cff).
Audit lineage: Addendum BQ (Kimi) -> this wave (Claude) -> Roger.
