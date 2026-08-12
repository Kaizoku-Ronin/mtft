# Study W1 — The compressed Weil form: port assessment, in-session certification, and one corpus correction

Date: 2026-08-10. Status: open study; module `weil.py` staged for v0.15.0, gated on the Kimi audit below. Source: "More than two thirds of the zeros of the Riemann zeta function lie on the critical line" (Claude / Anthropic, dated 2026-08-10; uploaded PDF).

## 0. Provenance and epistemic scope

The source paper is dated today and claims Lean 4 verification of its Theorems A–E. The cited repository exists and is public: `github.com/anthropics/zeta-23-lean` (HEAD `3635e74`, branches `main` and `rc2`, observed 2026-08-10 via `git ls-remote`). I did not build the Lean artifact in-session. Accordingly the asymptotic Theorems A–E carry the tag **[Ext]** throughout this study — external claims, plausible and formally attested elsewhere, but not corpus-certified. Nothing below depends on their truth. What this study certifies are **finite** identities and inequalities of the paper's machinery, each verifiable at machine precision at stated parameters, which is exactly the layer that can be brought home to mtft.

## 1. CC-02 (candidate): the Dirichlet series of w_n

Reading the paper's prime side (mean values over Λ(n)/√n and Λ(n)²/n) against the corpus surfaced an internal discrepancy in the corpus's own statements of the w_n generating function. Three printed values exist:

- Paper 1, Prop. 1.5, eq. (4): Σ wₙ n⁻ˢ = −ζ′(s+1).
- Arithmetica Generale, Pr 4.1.4: Σ wₙ n⁻ˢ = −ζ′(s)/ζ(s−1), "equivalently F(s+1)/ζ(s)".
- Forced by AG Pr 4.1.3 (the settled F(s) = −ζ(s−1)ζ′(s), re-confirmed here) together with wₙ = f(n)/n: Σ wₙ n⁻ˢ = F(s+1) = **−ζ(s)·ζ′(s+1)**.

Adjudication, three routes sharing no steps (E2 satisfied):

1. Direct sieve of wₙ = Σ_{d|n} (log d)/d, n ≤ 3·10⁵, versus the convolution sieve w = (Λ/id) ⋆ σ₋₁ (equivalently w = Λ₁ ⋆ 1, Λ₁(n) = (log n)/n): max pointwise difference **7.1·10⁻¹⁵**.
2. Numeric series at s = 3: partial sum 0.0828352629; −ζ(3)ζ′(4) = 0.0828352629, |diff| = **5.2·10⁻¹²**, exactly the size of the estimated tail. The printed alternatives miss by 1.39·10⁻² and 3.76·10⁻² respectively — excluded.
3. Consistency with the settled F: the same partial sum equals the F(4) partial sum identically (wₙ n⁻³ = f(n) n⁻⁴), and F(4) matches −ζ(3)ζ′(4) to the same precision.

**Corrected statement [Pr]: Σ_{n≥1} wₙ n⁻ˢ = F(s+1) = −ζ(s)·ζ′(s+1).**

Diagnosis of the error: −ζ′(s+1) is the Dirichlet series of the *summand* Λ₁(n) = (log n)/n. Since wₙ = (Λ₁ ⋆ 1)(n), the divisor sum multiplies by ζ(s); Paper 1's own proof line contains the ζ(s) factor and the final equality drops it (the same pathology as CC-01: a printed conclusion contradicting its own intermediate line). AG Pr 4.1.4's "−ζ′(s)/ζ(s−1)" is additionally inconsistent with its companion clause "F(s+1)/ζ(s)" (which itself is wrong — no division). Per protocol both statements receive append-only annotations; nothing is rewritten.

Downstream exposure: expected nil. Paper 18's explicit-formula pipeline and the Laplace-ensemble kernel use F(s) = −ζ(s−1)ζ′(s), which is correct and unaffected. A corpus grep for uses of "−ζ′(s+1)" as the w-series is assigned to Kimi below. Bridge note: the corrected series places ζ′ at the shifted argument s+1 inside the w-generating function — the object of the Speiser–Hadamard lab, and of the source paper's ξ′ extension (its Remark 7.3). See W3 below.

## 2. What the compressed Weil form is

Weil's explicit formula, read as a Hermitian form W(f,g) = Σ_ρ m_ρ ĥ_f(γ_ρ)·(pair-conjugate), is positive on all test functions iff RH. The paper compresses W to a d×d matrix G on a Gabor family φ(u)e^{−iτ_k u} at critical density (d ≈ λN test functions over [T, 2T]), and G has two independent computable expressions:

- **prime side**: G_kl = ∫ φ̂(τ−τ_k) φ̂(τ−τ_l) ν_X(τ) dτ, with ν_X = μ + Π_X + P_X built from Γ′/Γ, the pole, and prime powers n ≤ X = (T/2π)^λ — no zeros anywhere;
- **zero side**: G_kl = Σ_ρ m_ρ φ̂(γ_ρ−τ_k)·(conjugate structure), a sum over actual zeros.

Their equality is Weil's formula itself, and it is an exact, finitely checkable identity — a native E2 object. Everything downstream (traces, the ratio C = (tr G̃)²/tr G̃², rank–trace and Cauchy–Schwarz certificates, inertia bounds on synthetic configurations) is linear algebra on G.

## 3. In-session certifications (window I = [100, 200], λ = 1, η = 0.2 taper, unless stated)

| quantity | value | class |
|---|---|---|
| φ̂ quadrature vs mpmath quad (3 spot values) | max err 9.4·10⁻¹⁴ | Cert |
| Poisson/Gabor frame identity (Lemma 2.2), K = 4000 | rel err ≤ 6.2·10⁻⁵ | CERTIFIED(6.2e-5) |
| taper constants a, b vs paper §8 laws (1−0.603η, 1−0.688η) | 0.879332 / 0.862448 vs 0.8794 / 0.8624 | Cert |
| **E2: prime-side vs zero-side G** (371 zeros, γ ≤ 640.7, dps 15) | max rel discrepancy **3.386·10⁻⁶** | CERTIFIED(3.4e-6) |
| traces both routes | tr: 339.061206 vs 339.061905; tr G²: 3336.9952 vs 3337.0089 | Cert |
| N(I) exact vs Riemann–von Mangoldt main term | 50 vs 50.19 | Cert |
| measured C/N(I) vs finite-T law λ₁a²/(b+λ₁²J_T) | 0.6890 vs 0.6807 (both < F(1) = 0.75, per source Rmk 5.9) | DIAGNOSTIC |
| spectrum of G̃ | min 0.0021 > 0, n₊ = 44 = d (RH-verified range: PSD) | DIAGNOSTIC |
| rank–trace certificate 4trÂ − 2N(I′) − ‖Â‖²_F | +0.155·N(I) (truth s₁ = N(I′) = 60) | Cert (inequality exact on data) |
| Cauchy–Schwarz certificate | +0.178·N(I); ‖Â‖²_F/trÂ = 1.4615 (limit 4/3) | Cert |
| Lemma 3.2 audit (10⁵ random + structured instances) | 0 violations; equality configs to 1.4·10⁻¹⁴ | Cert on tested set |
| Montgomery–Taylor constant c*₁ | 0.753296067856; 2−1/c*₁ = 0.672500703679 | Cert (closed form) |
| MT window maximality (functional on grid) | flat 0.7499, cos(√2s) 0.7532, perturbed 0.7524 < MT | DIAGNOSTIC |
| synthetic inertia, p < d configs (corrected blocks) | 20 pairs @0.30: n₊=18≤20, n₋=3≤20; 10 pairs + 8 doubles: n₊=16≤18; cert ≤ s₁ in all | Cert on tested configs |

Notes. (i) At this height the Cauchy–Schwarz certificate exceeds the rank–trace one; the asymptotic ordering (rank–trace → 2/3 vs CS → 1/2) reverses only as ‖Â‖²/trÂ descends to 4/3. (ii) All classes are per-window statements; nothing here is a global bound.

## 4. Correction record within this study (preserved, per house rule)

Two errors of mine were caught and fixed during the session; both are kept in the record because each is a working lesson.

**W1-c1 (transcription).** My first pass at the paper's variational functional (7.3) read the OCR-ambiguous numerator as ∫v²; the functional then exceeded the closed-form maximum at the claimed maximizer. Eq. (7.2) disambiguates: the numerator is (∫v)² (it is a², not b). Corrected, the functional reproduces c*₁ to 1.2·10⁻⁴ and perturbations fall below it. Lesson: derive constants from the assembled formula, not from a single OCR'd display.

**W1-c2 (pair-block conjugation).** The first driver implemented the off-line pair contribution as 2m·Re(a·conj(a)ᵀ), which is positive semidefinite; the correct signature-(1,1) block, from the source's Prop. 4.1 proof (the conjugated factor is evaluated at γ̄_ρ), is **2m·Re(a·aᵀ) = xxᵀ − yyᵀ**. The bug was invisible in the first synthetic run because every configuration there had s₁+s₂+p ≥ d, so the dimension cap n₊ ≤ d made the bound vacuously true. The module's test suite, which includes a p < d configuration, exposed it immediately (n₊ = 44 > 30). Fixed; all inequalities now hold non-vacuously, and the negative index n₋ ≤ p is verified as well. Lesson: an inequality test whose bound exceeds the ambient dimension certifies nothing — pre-register the non-vacuous regime.

## 5. Portability map

**(a) `weil.py` — a fourth column for the ensembles program.** The Three Ensembles table (Laplace / Dirichlet / Critical) gains a natural fourth entry: kernel Φ(γ−γ′)² band-limited to |α| ≤ λ on zero pairs; localization = the bandwidth window [T, 2T]; criterion = inertia of the compressed form (n₋(G window) = 0 for all windows is the RH-side reading; the unconditional content is the pair of certificates). The prime-side density ν_X is a new computable object that sits beside the Dirichlet ensemble's resolvent data. The staged module provides: `Window` (taper, transforms, complex-argument φ̂), `gabor`, `nu_parts`, `G_prime`, `G_zero` (with synthetic depths and multiplicities), `certificates`, `rank_trace_gap`/`audit_rank_trace`, `mt_constant`, `w_series_check` (the CC-02 adjudicator), `selftest`. Sympy-free; scipy optional (complex digamma fast path with an mpmath coarse-grid fallback).

**(b) CC-02 → Speiser bridge (W3 proposal).** With Σwₙn⁻ˢ = −ζ(s)ζ′(s+1), the analytic structure of the w-generating function is governed jointly by the zeros of ζ (at s = ρ) and of ζ′ (at s = ρ′−1). The Speiser–Hadamard lab already certifies the ζ′ census to γ < 100; the source paper's Remark 7.3 runs the compression on ξ′ with stated constants 0.85838 / 0.92919 [Ext]. W3: compress the Weil form of ξ′ with the lab's census as the zero side — a finite E2 target entirely within existing corpus objects.

**(c) Theorem E degree-1 route (W2 proposal).** For a fixed primitive Dirichlet character the same engine runs with ν_{X,χ} (conductor term in μ_χ, χ(n)-twisted P, no pole term) [Ext for the asymptotics; the finite matrices are ours to compute]. The natural mtft instance is the character family mod 11, 13, 143 — the degree-1 harmonics of level 143, adjacent to the certified Eisenstein congruence structure in `eisenstein.py`. Deliverable: finite-T certificates per character; pre-register before running.

**(d) Falsify-engine additions.** Four pre-registrations, stated before any further computation:

- **W1-P1**: on I = [300, 600], λ = 1, η = 0.1, with all zeros γ ≤ 1000 at dps ≥ 15, the prime/zero max relative discrepancy is ≤ 1·10⁻⁵.
- **W1-P2** (external anchor): on I = [2000, 4000], λ = 1, η = 0.05, C_G̃/N ∈ [0.729, 0.739] (source §8(3) reports 0.734). ~3700 zeros; overnight local job on the laptop.
- **W1-P3**: for 100 random ρ↔1−ρ̄-symmetric synthetic configurations with depths U[10⁻⁴, 0.45], 20% doubled multiplicities, and **s₁+s₂+p < d enforced**, zero violations of n₊ ≤ s₁+s₂+p, n₋ ≤ p, and cert ≤ s₁.
- **W1-P4**: Lemma 3.2 on 10⁶ random plus ≥ 300 adversarially optimized Hermitian instances: zero violations; equality configurations to ≤ 10⁻¹².

## 6. WI-N1 — honest negative, filed before anyone burns a week

**No direct critical-line certificate exists for the X₀(143) newform L-functions (f₁, f₂, f₃ or any member of the three orbits) by this method, for any window.** Mechanism, not proximity: the Montgomery–Vaughan mean-value step caps the arithmetic length at X ≤ T^{1−ε}, so the normalised bandwidth cannot exceed Λ* = 1/m for degree m; at m = 2 the best achievable constant is c = (1/2)/(1 + 1/12) = 6/13, and the certified proportion 2 − 1/c = −1/6 < 0 — the certificate is empty whatever the window (source Rmk 7.2(ii), a structural statement of the same flavour as the M6/M7 reality obstruction: the wall is the degree, not the implementation). Escape hatches noted by the source and left open: family averaging (q → ∞, not our fixed level) and higher-correlation inputs of Hardy–Littlewood strength [Ext, "not carried out here"]. Corpus posture: the degree-1 route (W2) is the honest port; degree-2 attempts are pre-declared dead.

## 7. Kimi handoff (gate for v0.15.0 integration)

1. Independently re-derive Lemma 3.2 (von Neumann trace inequality + x² ≥ cx − c²/4; short) and its multiplicity-aware form; confirm the equality case.
2. Reproduce the E2 identity with an implementation sharing no code with `weil.py` (the source paper's §2 is a complete recipe); target W1-P1's window and tolerance.
3. Confirm CC-02 by an independent route (e.g., Abel summation on Σ (log d)/d^{s+1} · ζ(s) or direct mpmath evaluation at a complex s); then grep the corpus for any downstream use of "−ζ′(s+1)" as the w-series (expected: none; F-based results are unaffected).
4. Run the W1-P2 anchor locally (overnight); report C_G̃/N against the pre-registered interval.
5. Adversarial pass on `weil.py` itself, with W1-c2 as the worked example of the failure mode to hunt: any inequality test whose bound is ≥ d.

Files: `weil.py` (candidate module), `test_weil.py` (7 tests, all passing; E2 test uses the shipped zero cache), `w1_study_driver.py` (clean reproducer of every number in §3), `zeros_gamma_T100.npy` (371 ordinates, γ ≤ 640.69, mpmath dps 15).

Windows note for local runs: `py -m pytest`, `--break-system-packages` where needed; the zero sweep for W1-P2 should checkpoint to `.npy` every 200 zeros.
