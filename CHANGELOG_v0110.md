# mtft v0.11.0 — X₀(143) Particle Box (the three-generation spectral atom)

**The X₀(143) particle box (R. Tano, v01–v03 lineage) joins the package,
audited end-to-end before landing (Addendum BQ, independent auditor
K. K3: own Manin-symbol/Hecke chain, own aₚ engine, own Rankin–Selberg
pipeline, own quadrature engine, positive controls at level 11 and at
the 11.a2 Petersson norm). Package version: 0.11.0 across
pyproject.toml, `__init__.py`, CITATION.cff — three-way guard
satisfied.**

X₀(143) (level 143 = 11×13, index 168, genus 13) carries a weight-2
cuspidal spectrum splitting into blocks of dimensions 2/8/12 over the
three irreducible factors of the T₂ polynomial — three "generations"
(f1 = the elliptic curve 143.a1, f2 a quartet, f3 a sextet). The box
builds the physics on top: four cusps as nuclei (capture requires
dissipation, with an exact γ₀ = 0 control theorem and a capture ceiling
of exactly 1/2), orbit-dependent Zeno rates on the 26-dim cuspidal
homology via the Petersson metric, and an AG-D5 falsifiability battery
whose verdict ("ordering only") survives verbatim.

## New — `studies/x0143_particle_box{,_v02,_v03}.py`

- v01: tessellation engine (P¹(ℤ/143), Farey triangulation: 56
  triangles / 84 edges / 4 cusps of widths {1, 11, 13, 143}) and the
  exact Manin-symbol machinery (29-dim quotient, 26-dim cuspidal).
- v02: float Hecke engine, period functionals, aₚ extraction,
  falsifiability inputs.
- v03: stage G (cusp nuclei, Lindblad capture, emission spectrum,
  photon ledger) and stage H (Petersson metric on homology,
  orbit-dependent Zeno, falsifiability battery H4).
- Runtime outputs (certificates, figures, period cache) anchored to
  `__file__` per the stage-5 convention.

## Fixed — `mtft.x0_143` oracle (Add. BQ §5)

- `CURVE_143A1.j_invariant`: `'-1/15'` → `'-262144/1859'`. Certified
  exact from the a-invariants (c₄ = 64, Δ = −1859 = −11·13²) by the
  auditor; the previous value was a corpus error (corpus_flags[1]
  disposition). The second corpus flag (Paper32 τ) was likewise
  adjudicated against the corpus; the engine's τ is correct.

## Documented (Add. BQ §3, §8)

- The PET_F1 / PET_F2_DIAG / PET_F3_DIAG constants are **per-unit-volume**
  Petersson norms: ⟨f,f⟩_raw = PET · V_N, V_N = 56π. Comment added at
  the definitions in v02.
- Auditor certification of the norms (Add. BQ §4.4): PET_F1 sits 4.7%
  above the exact modular-degree value (deg(143.a1) = 4); the f2/f3
  diagonals match the auditor's independent Rankin–Selberg residues to
  a uniform 2–3% on all eleven embeddings. PET_F2/F3 remain DIAGNOSTIC
  class, unchanged.
- Photon-ledger figure label corrected (the 7.8e−16 is a telescoping
  construction identity; the quadrature residual is 7.2e−3).

## Tests — `tests/test_x0143_particle_box.py`

Fast tier (default): tessellation combinatorics; cuspidal T₂
characteristic polynomial (exact factor string); f1 Hecke eigenvalues
from the Manin engine vs the certified oracle table (p ≤ 47); oracle
j-invariant value; G1 cusp-incidence integer certificates; G4 unitary
control (γ₀ = 0 ⇒ generator identically zero).
Slow tier (−m slow): capture-ceiling smoke (full stage-G run).

## Regression

- New tier: 6 passed, 1 skipped. Oracle tier
  (`test_x0_143_verified.py`): 21 passed. Full repo regression:
  443 passed, 1 skipped, 0 failed.
- sin²θ_W = 3/13 = 0.230769 vs measured 0.23122 remains presented
  with the tension stated, per the packet's own framing (Add. BQ §9.6).

Audit trail: Addendum BQ (this wave, full independent audit, 8 legs,
conditional GO). Open items unchanged: BI.F1 (PARI mfpetersson diagonal
certification — now optional rather than motivated by suspicion);
PR-37 resumes on the integrated package.


---

### Corrections applied in v0.11.1 (annotation appended per corpus protocol)

- The "Documented" section above states the PET normalization as
  per-unit-volume (V_N = 56pi) and PET_F1 as "+4.7% above" the
  modular-degree value. Both are superseded: mfpetersson normalizes by
  the INDEX 168; the ratio 56pi/168 = pi/3 accounts for the whole 4.7%;
  the modular-degree route matches PET_F1 to 7.8e-14 (all published
  digits). See CHANGELOG_v0111.md.
- The photon-ledger relabel recorded above was documented but not landed
  in the v0.11.0 source; it is landed in v0.11.1.
