# mtft v0.10.0 — The Spectral Reconstruction Toolkit (stages 1–5 of 5)

**FINAL — all five stages of the Integration Plan (v0.1) landed and
audited (Add. AZ → BP; stage 5 audited Add. BP, gauge certificate
independently reproduced, S5-1 dispositioned). Package version:
0.10.0 (pyproject.toml, `__init__.py`, CITATION.cff all bumped —
three-way guard satisfied).**

This release migrates the marked-primon-gas spectral program (the
rung-4/rung-5 corpus: construction note through PR-36, audit record
Add. U → BG, integration audits Add. AZ → BN) into the package, per the
Integration Plan §3 layout. Five new modules, all small — the machinery
is compact once de-duplicated (16 copies of `internal()` and 9 of
`gsq()` across the suite scripts collapse into one implementation each;
two copies of the A/C/C³ closed forms collapse into one, Add. BN §6).

Every module was built under the mutual-audit protocol: the builder's
code is independently audited before first use, findings are verified
against the auditor's own engines before acceptance, and errors are
owned in writing on both sides (audit trail: Add. AZ → BN).

## New — `mtft.ledger` (stage 1; K. K3; audited Add. AZ; extended Add. BN)

Every certified constant of the spectral program **as data**, with
value, uncertainty, class, and provenance — the fix for the
stale-constant problem (Plan §2: five superseded values hardcoded
across ten suite scripts). Superseded values are recorded, not
deleted; every entry names the note and addendum that certified it.

- `LEDGER`, `Entry`, families A and C; closed-form accessors `A_of`,
  `C_of`, `C3_of`, `u_ep_limit`, response coefficients.
- `verify()` — 46 checks recomputing every entry from its closed form
  (the ledger is self-auditing; nothing downstream may hardcode a
  number that lives here).
- `RELATIONS` — derived identities between observables, with
  derivation and provenance. First entry: `dev_centre_vs_dev_EP`
  (S3-1, resolved Add. BK/BL — see below).
- **Stage-4 registrations (Add. BN):** `mu_0`/`mu_1` (dressed diagonals
  at κ* = 5 — the BM proposal, accepted as amended; certification path
  is N-converged f64, and the entries carry the BN-E1 warning that
  low-nsite mp under-reads) and `tau_c_star`/`V_b_tree` (the coupled
  model's band-merging coupling and tree binding threshold, so
  `coupled.selftest` can assert through the `_L` guard — BN-F1).
- **Duplication resolved (Add. BN §6):** `A_of`/`C_of`/`C3_of`/
  `C3_channels` are now lazy delegates to `mtft.expansion` — one closed
  form, one home; the ledger keeps the data and the `verify()`.

## New — `mtft.chain` (stage 2; R. Tano; audited Add. BI; fixes verified Add. BN)

The one kernel pipeline: sites n = 2..N, β = 2, ρ(n) = ln n·n⁻³,
K = (min/max)^κ, T = DKD. `internal()` builds the gaps g_i and the
dressed hopping B = e^(−h/2)(VᵀKV)e^(−h/2) once, at either precision
tier — `backend='f64'` for exploration, `backend='mp'` (closed-form T)
for certification. The crossover is measured and enforced in the API:
`F64_EP_FLOOR = 2×10⁻¹²`, `F64_WINDING_FLOOR = 2×10⁻⁴` (PR-20 §2 —
inherited float64 g and B cap any measurement regardless of working
precision). `gap(i)` is continuous-i; `crossing_limit(i)` gives the
strong-coupling EP magnitudes that seed every extraction.

**BI-queue fixes landed and verified (Add. BN §4):**
- **BI-F2 fixed** — chain/expansion selftests now assert through the
  `_L`/`_LF` ledger guards, which raise `KeyError` on any unregistered
  name; a selftest may not certify against an unregistered number.
- **BI-F3** — `extract()` control-stencil default on; validation never
  skipped. **BI-F4** — completeness walk-counting uses full
  connectivity (the ±3 cutoff missed jump-4 walks). **BI-F5** —
  `scan_zeros`' blind spot documented.

## New — `mtft.expansion` (stage 2; R. Tano; audited Add. BI; fixes verified Add. BN)

Cauchy/Taylor monomial machinery: `extract()` with the `/4` stencil
correct and per-stencil positive controls; the closed forms A_i, C_i,
C³_i with channel decompositions and cancellation ratios — now the
**single implementation** the ledger delegates to (Add. BN §6).

**One required fix remains with the author (Add. BI, still open):**
- **BI.F1** — `richardson()`'s recurrence reintroduces c₁h₀²/15 for
  ≥ 3 points (18× measured degradation). **Blocked from the public
  surface in this integration** (not re-exported in `__init__`). The
  Neville-tableau fix (auditor-supplied, author-landed) is verified
  exact to 8.9×10⁻¹⁶ on the auditor's battery (Add. BN §4); the
  re-export decision is the author's.

## New — `mtft.ep` (stage 3; K. K3; audited Add. BK, dispositioned Add. BL)

Exceptional points: `gsq` (closest-pair selection only — the
fixed-index variant is not exported; an index is not an identity),
`newton` (explicit FD; no `mp.diff` on implicit observables),
`levels_of`, `ep_of`, `pair_winding`, argument-principle `census`,
`staircase`, `nearest` (requires `levels=` — census ≠ search).

Measured structure encoded, not intuited:

- **Finite-κ crossings are avoided** (min gaps 10⁻⁵–10⁻² at κ = 30).
  `ep_of` extracts the **diabatic centre** — the minimum of the
  closest-pair gap² — by Newton on gsq′, never on gsq (no zero exists).
- **S3-1 resolved (Add. BK derivation; Add. BL independent
  verification):** the sealed reduced-block hierarchy describes the
  complex EP branch point; the diabatic centre differs by exactly one
  copy of the within-pair repulsion constituent —
  `dev_centre(i;κ) = dev_EP(i;κ) − rep_C(i)·s²ᵏ`, from the 2×2
  dressed-block identity (EP: −2mv²/D³; centre: −4mv²/D³), with
  v² = e^(−(a+b))s²ᵏ(1 + O(10⁻⁶)). Full-model residuals ≤ 1.2×10⁻⁶
  relative at κ = 60, rungs 1–4 — now asserted in the suite, and
  registered as `ledger.RELATIONS['dev_centre_vs_dev_EP']`.
- **Staircase structure:** at rung i's EP the crossing pair is always
  sorted pair (0,1); the rung's identity lives in the **characters**
  (diabatic-basis labels of the exchanging branches), certified by
  `pair_winding`'s two witnesses.
- **Census counts EP pairs, not rungs** — complex EPs attach to every
  diabatic pair crossing; wide strips reproduce the closed-form
  diabatic count exactly, tight strips are γ-selective, and a zero
  count now carries a not-evidence-of-absence note (BK-F3).

Disposition of Add. BK (Add. BL, all findings reproduced, owned, and
fixed by the builder):

- **BK-F1** — at small κ the veers overlap and *no character-pure
  endpoint exists* (s_far non-monotone in half_width; measured
  0.16/0.64/0.002/0.12 at κ* = 5). `pair_winding` now auto-widens
  (×2, then ×1.25 refinement — pure doubling overshoots crowded
  windows) and its refusal distinguishes the two diagnoses by the
  measured trajectory: monotone → genuine half_width problem;
  non-monotone → regime limit, message says so. Second genuine
  instance found during disposition: rung 4 at κ = 30.
- **BK-F2** — frame-to-frame label tracking **anti-converges** (finer
  stepping follows the adiabatic branch through each veer; rung 3's
  characters flipped between 80 and 160 steps and stayed wrong at
  400/800 — a resolution-stability gate would certify the
  converged-wrong regime). Fixed at the root: labels are read against
  the fixed strong-coupling basis (`argmax |V[a,j]|`), parameter-free.
  The crowded-window physics this exposed (rung 3's clean window
  necessarily spans spectator diabatic crossings) redrew the
  certificate's spectator semantics: leak = pair↔spectator only,
  consistency = boundary crossing, spectator swaps reported as real
  pencil structure.

## New — `mtft.coupled` (stage 4; R. Tano; audited Add. BN)

One Hamiltonian, two measures: H(x) = h(κ) − τxB with x drawn from the
**Bloch** measure (ℤ: x = 2cos k, arcsine density, hard edges) or the
**Kesten** measure (coordination-p tree: Kesten–McKay density on
[−2√p, 2√p], soft edges). The vacuum-dressing rule is structural —
`H_of` takes `chain.Internal` and uses its dressed `.B`, never
re-derives it.

- `Measure` (frozen dataclass) with Gauss–Legendre `nodes()`; the tree
  θ-substitution cancels the √ endpoint singularity exactly.
- `band`/`bands`, `tau_c` (band-merging coupling, bisection with
  bracket check), `binding_threshold` (1 = V·G_dd(E); ℤ G diverges with
  exponent −1/2 — any attraction binds; tree G converges — V_b finite),
  `moments`, `_tree_walks` (closed-walk counts: m₂/m₄/m₆ = q,
  q(2q−1), q(5q²−6q+2)).
- Audited numbers (Add. BN, all independently re-derived): Kesten
  moments 3/15/87 at q = 3; **τ_c ratio tree/ℤ = 1/√2 to 12 digits**
  (PR-8 geometry-free scaling: H depends only on τx, so τ_c·x_max is
  measure-independent); **tree V_b = 0.0373264 finite** vs ℤ exponent
  −0.5004 (theory −1/2). Certified constants registered as
  `ledger.tau_c_star` / `ledger.V_b_tree`.

**BN findings closed by the author's patch (Add. BO; see "Closed this
release" below).**  [Audit note, Add. BP-F2: the author's FINAL draft
reverted this section to the pre-BO "open" wording while also adding
the "Closed" section — the contradiction was reconciled at landing.]

## Rules as code (enforced, not documented)

Closest-pair selection only · precision tiers with measured floors ·
no `mp.diff` on implicit observables · `nearest()` requires `levels=` ·
census is a counting theorem, not a scan · positive controls
mandatory · superseded values recorded · **tests assert against the
ledger, never literals** — now enforced by the `_L` guard, which
refuses unregistered names · the Legend is generated from the ledger
(`src/mtft/mtft_legend.md`, shipped as package data).

## Tests

`tests/test_spectral_stage1_3.py` wraps each module's certification
suite in the Plan §6 tiers: fast (`ledger.verify` 46 checks; chain,
expansion, and coupled selftests — 15 checks ≈ 2 s) and slow-marked
(`ep.selftest`, 62 checks ≈ 75 s, mp-tier winding certificates and the
resolved S3-1 relation). Full repo regression: 446 passed (fast tier),
spectral suite 5/5 including slow.

## New — `studies/` (stage 5; R. Tano; audited Add. BP — gauge certificate independently reproduced leg by leg)

The 33 suite scripts (OP3 → PR-34, rung studies, census drivers)
migrated from session artifacts into `studies/`, importing the package.

- **16 `internal()` copies → `mtft.chain.internal`** via three-line
  adapters preserving each study's signature and `(g, B)` return.
  Certificate: g to 9e-15; full H(u) spectrum to 2.5e-14 over
  u ∈ [0,8], κ ∈ {5,12,30,60}. The B representation differs at low κ
  (up to 1e-2) — pure gauge inside near-degenerate blocks, seeded by
  the session scripts' Z₂-normalised ρ vs the package's unnormalised
  form (a global scalar on T; operator unchanged). All observables
  engine-independent: gsq 7e-15, μ₀/μ₁ and band edges 1e-15, BG4
  survival P(t) 2.2e-15 (analytic gauge-invariance χ→χR, U→RᵀU
  confirmed empirically).
- **rung5_bloch's five-value return** rebuilt from `Internal` fields
  (`.K_raw`, `.V`) — nothing re-derived.
- **pr24/pr25 closed forms → `mtft.expansion.A/C/C3/channels_C3`**,
  equal to 30 dps at every tested index (completes Add. BN §6).
- **S5-1 (NEW finding; dispositioned Add. BP):** nine local
  gsq copies evaluate complex u in backends `mtft.ep` deliberately
  does not offer (f64 is real-u only; PR-20 floors, real-axis
  diabatic-centre design, Add. BK/BL). They are FROZEN verbatim with
  markers so each study's historical record stands. [Audit note,
  Add. BP-F1: the author's enumeration said eight; pr34's
  self-contained mp 3×3 minimal-block gsq is the ninth.] **Disposition
  (auditor's call, Add. BP): the frozen copies are the permanent
  record.** `ep.gsq` does not grow a complex-u path — at complex u
  the non-Hermitian gap is a different observable (sort-by-Re,
  min |Δλ|²), not a backend variant of the certified real-axis
  selection rule. Any future complex-u need gets a new function with
  its own definition and certification arc, not a gsq retrofit.
- Drivers' sibling file-loads anchored to `__file__`; studies run
  from any CWD. Excluded from pytest (they are measurement scripts,
  not tests; some run for hours by design).

## Closed this release (previously queued)

- **BN-F1 closed** — `coupled.selftest` asserts through
  `_L('tau_c_star')`/`_L('V_b_tree')`; verified 15/15 green.
- **BN-F2 closed** — quadrature cache is single-`eigh` per node;
  identical physics (V_b 0.0373263988, exponent −0.5004).

## Remaining author decision (not a release gate)

- **BI.F1** — `expansion.richardson` (Neville fix verified to
  8.9e-16) is still not re-exported from `__init__`; flip when ready.

## Error ownership (house rules, this release cycle)

- Auditor-owned: BJ-E1/E2/E3 (stage-3 test bugs — positional argument,
  index-as-identity assertion, near/far convention), BL-E1 (a
  gate-loosening misstep during BK disposition, caught pre-delivery),
  the BK-F1/F2 root causes (refusal message that misdiagnosed the
  regime; the anti-convergent tracker — both "encoded the intuition
  instead of the measured structure"), **BN-E1** (judged the author's
  mu_1 wrong by 1.9×10⁻⁶ when his value was the converged one —
  arithmetic precision is not truncation convergence; the entry now
  carries the warning), **BN-E2** (sign slip on the first independent
  G_dd leg, caught by the n-convergence cross-check).
- Author-owned: BI.F1 (open, above; BI-F2 fixed and verified Add. BN);
  BN-F1/F2 (closed Add. BO); BK stand-in adapter gaps (owned in
  Add. BK).  Add. BP: two documentation findings against the stage-5
  wave (BP-F1 frozen-census count/banner, BP-F2 changelog reversion)
  — reconciled at landing, no code impact.

## Resumes now

- PR-37 (pre-registered) resumes on the integrated package.
