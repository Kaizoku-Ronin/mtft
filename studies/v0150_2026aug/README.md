# v0.15.0 study bundle — 2026-08-12

Wave artifacts from four development chats (V6 period matrix, W1 Weil
compression, arithmetic area geometry, M8/M8b/M9 complex structures), staged
for integration, plus the auditor's correction and adjudication artifacts.

**Byte-preservation policy.** Every wave artifact listed under "Wave files" is
staged exactly as received — including statements the audit found to be wrong.
Corrections live in the append-only `CC-0X_*.md` files and in the auditor
artifacts; nothing upstream was rewritten.

## Wave files (byte-preserved)

| file | content |
|---|---|
| `INTEGRATION(1).md` | area-geometry wave integration note (v0.14.1 header; folded into v0.15.0) |
| `arithmetic_area_geometry.py` + `_ledger.json` | rigidity theorem, Brioschi route, gates A1–A9 |
| `RIEMANN_MATRIX_V6_REPORT.md` | V6 period-matrix wave: 19 PASS / 1 FAIL (C9b) |
| `period_matrix_manin_v6.txt`, `x0143_period_data_v6.json`, `run_v6.log` | V6 period data artifacts |
| `mtft_period_matrix_manin_v6.gp.txt` | V6 GP transcript |
| `W1_weil_compression_study.md` | W1 Weil-form compression study + CC-02 adjudication |
| `w1_study_driver.py` | clean reproducer of every W1 number |
| `Tano_Weights_Fact_Ledger.md` | weights fact ledger (CC-05 source, §12.1) |
| `m8_hodge_channel.py` + `_ledger.json` | M8 Hecke-commutant study (see CC-06) |
| `m8b_arithmetic_complex_structure.py` + `_ledger.json` | M8b complex structure (amplitude superseded — CC-06) |
| `m9_oldspace_complex_structures.py` + `_ledger.json` | M9 oldspace census (certified) |

## Auditor artifacts

| file | content |
|---|---|
| `CC-02_wseries_shift_chain.md` | w-series = −ζ(s)ζ′(s+1), three-route E2 |
| `CC-03_adaptive_precision.md` | finite_atom_curvature silent-precision fix |
| `CC-04_cold_amplitude.md` | cold amplitude closed form; retracts 16th digit |
| `CC-05_chapter5_weight_table.md` | Chapter 5 Table 1 four-value correction |
| `CC-06_m8_basis_mixup.md` | M8 mixed-basis correction (proposal) |
| `m8_verify_kimi.py` / `.json` | six-prime commutant reproduction: dim = 2 |
| `m8_deep_kimi.py` / `.json` | basis-mixing diagnosis; Z projector |
| `c9b_exact_symbol.py`, `c9b_exact_symbol.json`, `c9b_routeB.json` | C9b two-route E2 adjudication |

## Internal → repo correction-label mapping

The wave files carry their own internal "CC" numbers. The repo ledger
(`src/mtft/hodge_polarization.py` docstring) is the single source of truth:

| internal label (file) | repo label |
|---|---|
| CC-02 w-series shift chain (W1 study §1) | **CC-02** (unchanged) |
| CC-03 adaptive precision (INTEGRATION(1).md) | **CC-03** (unchanged) |
| CC-02 cold amplitude (INTEGRATION(1).md, area study docstring) | **CC-04** |
| "Candidate CC-04" Chapter 5 table (fact ledger §12.1) | **CC-05** |
| M8 basis mix-up (auditor, this bundle) | **CC-06** |
| CC-02 du03 metric (studies/hodge_2026aug) | **CC-07** (renumbered 2026-08-12; note prepended to the file) |

## Disclosures

1. **A8 historical-demonstration gate.** `arithmetic_area_geometry.py` gate A8
   ("A8_dyadic_shield_and_CC03") imports `finite_atom_curvature` from
   `mtft.curvature` and *demonstrates the pre-patch failure mode*: it asserts
   that `dps=90` at β=80 returns garbage far from 1/4. The staged
   `src/mtft/curvature.py` is the **patched** module (dps may only raise
   precision; adaptive doubling), so A8's footgun clause cannot pass against
   the installed package. The study is kept byte-preserved as the historical
   record of the CC-03 failure; do not re-gate it against the patched module.
   Its other gates (A1–A7, A9) are module-independent and stand.

2. **M8/M8b mixed basis (CC-06).** The M8 study reports
   dim commutant({Hecke, V}) = 1; the auditor's verbatim reproduction on six
   primes gives **2** (= span{I, Z}, Z/2 = (I − ι*)/2). The canonical M8b
   amplitude is **0.20610964892935077**, not the 0.3805 in the m8b ledger.
   Wave files byte-preserved; corrected values are the certified ones.

3. **C9b adjudicated — E2 certified.** The V6 wave's single FAIL (the
   [∞, 2/77] f1-period sign against the Paper 33 v2 archive) is resolved by
   two routes sharing no steps: Route A (exact rational Manin symbols) gives
   Re λ₁ = −1/2 exactly; Route B (q-expansion slash-integrals, own aₙ sieve,
   no PARI `mfsymboleval`) reproduces V6's per₂₇₇ to ~1e-7 and matches both
   engines' per₁₁ exactly, with {1/38, 2/77} = −per₁₁ to 4e-15 (Manin
   integrality). **V6 (PARI 2.15.4) sign CORRECT; the Paper 33 v2 archive
   (PARI 2.17) needs the sign-only correction** per the V6 report's
   pre-registration. Full budgets in `c9b_routeB.json` / `c9b_exact_symbol.json`.

4. **W1 Weil bundle certified.** Lemma 3.2 re-proved independently; W1-P1 E2
   3.472e-09 (pre-registered ≤ 1e-5); W1-P2 C/N = 0.73398 in the
   pre-registered [0.729, 0.739]. `tests/zeros_gamma_T100.npy` was regenerated
   by the auditor from the certified zero sweep (3994 zeros to γ = 4500.16):
   371 ordinates, max γ = 640.6947946688257 — the study's stated parameters
   (371 zeros, γ ≤ 640.7) exactly.

5. **M9 census certified; R7 open.** The M9 oldspace census
   (48 rational degeneracy primes for 11a1 below 20000; rationality iff
   q = a′² + k′²) reproduced. The R7 Eisenstein-congruence correlation
   hypothesis is filed explicitly **open** — no claim is made.
