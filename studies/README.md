# studies/ — the rung-4/rung-5 suite, re-pointed (stage 5)

The 33 study scripts of the marked-primon-gas spectral program
(OP3 through PR-34, the rung studies, and the census/diagnostic
drivers), migrated from session artifacts into the repository per
Integration Plan v0.1 §3, stage 5.  These are the **measurement
records** behind the certified constants in `mtft.ledger` and the
audit trail Add. U → BN.

## Re-pointing policy (what changed, what is frozen)

**`internal()` → `mtft.chain.internal`** (16 local copies removed).
Each study keeps its original signature and `(g, B)` return through a
three-line adapter.  Equivalence certificate: `g` identical to 9e-15;
the full H(u) = diag(g) − uB spectrum identical to 2.5e-14 over
u ∈ [0, 8] at κ ∈ {5, 12, 30, 60}.  The B *representation* differs at
up to 1e-2 at low κ — a pure gauge rotation inside near-degenerate
blocks (the session scripts normalised ρ by Z₂ = −ζ′(3), the package
does not; the scalar perturbs LAPACK rounding, not the operator).
Every observable checked is engine-independent: gsq to 7e-15, B
diagonals (μ₀, μ₁) and band edges to 1e-15, and the BG4 survival
probability P(t) to 2.2e-15.

**`rung5_bloch_coupling.internal()`** returns five values
(g, B, M̃, V, M); all five are reconstructed from the `Internal`
dataclass fields (`.K_raw`, `.V`) — nothing re-derived.

**Closed forms → `mtft.expansion`** (pr24 `A_of`/`C_of`, pr25
`A_of`/`C3_of`).  Verified equal to the package `A`/`C`/`C3`/
`channels_C3` to 30 dps at every tested index (Add. BN §6: one closed
form, one home).

**`gsq()` copies: FROZEN, not re-pointed (finding S5-1).**  Nine
studies (nearcheck, pr14, pr15, pr19, pr20, pr21, pr22 `gsq3`, pr24,
and pr34 — pr34's is a self-contained mp 3×3 minimal-block engine)
evaluate gsq at *complex* u in the f64/mp backends they shipped with.
`mtft.ep`'s f64 path is deliberately real-u only (PR-20 floors; the
audited design extracts the real-axis diabatic centre, Add. BK/BL).
Re-pointing would put a different engine under a historical
measurement, so these copies are retained verbatim with an S5-1
marker.  **Disposition (Add. BP, auditor's call):** the frozen copies
stand as the permanent record.  `ep.gsq` does not grow a complex-u
path — at complex u the non-Hermitian "gap" is a different observable
(sort-by-Re, min |Δλ|²), not a backend variant of the certified
real-axis selection rule.  If future work needs complex-u scans from
the package, it gets a NEW function with an explicit complex-gap
definition and its own certification arc, not a gsq retrofit.

**Drivers** (`census2`, `ks3`, `nearcheck`) load sibling studies by
file path; the paths are now anchored to `__file__` so they run from
any working directory.

## The X₀(143) particle box (v0.11.0)

`x0143_particle_box.py`, `x0143_particle_box_v02.py`,
`x0143_particle_box_v03.py` — the three-generation spectral atom on
X₀(143): tessellation + exact Manin-symbol engine (v01); float Hecke,
periods, aₚ extraction (v02); cusp nuclei + Lindblad capture +
Petersson-metric Zeno + falsifiability battery (v03).  Audited
end-to-end in Addendum BQ; gates in `tests/test_x0143_particle_box.py`
(fast tier runs in ~2 s).  Outputs (certificates, figures, period
cache) are written next to the scripts.

## Running

Studies are scripts, not tests: many run for minutes to hours by
design (nearcheck's Newton grid alone is ~10⁶ eigendecompositions).
They require the package importable (`pip install mtft` or
`pip install -e .` from the repo root), then:

    python studies/pr17.py

They are excluded from the pytest surface; the machinery they exercise
is certified by the package selftests
(`test_spectral_stage1_3.py`).  `rung5_v1_postmortem_ARCHIVE.py` is an
archived falsification record and is not meant to be re-run.

## Provenance

Original artifacts: rung-4/rung-5 sessions, Add. U → BN.  Builders:
R. Tano (Claude engine) and K. K3 under the mutual-audit protocol.
Stage-5 migration certified in CHANGELOG_v0100.md.
