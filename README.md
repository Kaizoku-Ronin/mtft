# MTFT — Modular Time Field Theory

[![PyPI](https://img.shields.io/pypi/v/mtft.svg?color=blue)](https://pypi.org/project/mtft/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://img.shields.io/badge/tests-401%2F401-brightgreen.svg)](https://github.com/Kaizoku-Ronin/mtft/actions)

Fundamental constants from the integers. **Zero free parameters.**

![The Stiffness Landscape](https://raw.githubusercontent.com/Kaizoku-Ronin/mtft/main/viz/hero_stiffness.png)

<p align="center"><sub>
μ<sub>N</sub>(y) over gauge groups SU(2)…SU(16), computed entirely from
w<sub>n</sub> = Σ<sub>d|n</sub> (log d)/d. The wall at small y is confinement;
uniform positivity <em>is</em> the mass-gap statement.
<a href="https://raw.githack.com/Kaizoku-Ronin/mtft/main/viz/stiffness_navigator.html">
<b>Open the interactive 3D navigator →</b></a>
</sub></p>

MTFT proposes a modular time field **τ(x) = t_R / t_U** mapping spacetime to
the upper half-plane **ℍ**, with SL(2,ℤ) symmetry and the modular curve
**X₀(143)** (genus 13, level 11×13) as the arithmetic backbone. Coupling
constants, the gauge tower, dark-sector profiles, and a Riemann-Hypothesis
toolkit all emerge from one weight sequence — no fitted parameters anywhere
in the pipeline.

## The Three-Ensemble Program (v0.7 → v0.8)

One weight sequence, three assemblies, three classical RH criteria — all
installable, all cross-certified by independent audit:

| Ensemble | Object | Curvature output | RH criterion |
|---|---|---|---|
| **Laplace** | Σ wₙ e^(−2πyn) | μ_N(y) = (1/4π²)[T″(y) − Re T″(y−i/N)] — exact | Th 1: limsup\|Δκ·X^(−3/2)\| < ∞ |
| **Dirichlet** | Z_D(β) = −ζ(β)ζ′(β+1) — exact | g_D = ∂²log ζ(β) + ∂²log(−ζ′(β+1)) | Speiser (1935): ζ′ ≠ 0 in 0 < Re s < ½ |
| **Critical** | log ξ at s = 1 | Li coefficients λₙ, three independent methods | Li (1997): λₙ ≥ 0 for **all** n |

```python
from mtft import dirichlet_curvature, li_criterion_report, hadamard_zetaprime_check

dirichlet_curvature(3.0)     # exact split: ζ piece + ζ′ piece (48.59%)
li_criterion_report(12)      # λ₁..λ₁₂ with the Bombieri–Lagarias caveat attached
hadamard_zetaprime_check(3)  # ∂²log(−ζ′) = 2/(s−1)² − Σ(s−ρ′)⁻²; residual 4e−9
                             # at s=3, certified < 1e−5 across s ∈ [3,10]
```

## The Legend — `python -m mtft.legend`

mtft computes artifacts of the number line. The Legend is its map key:
every function tagged by **nature**, **Arithmetica Generale primitive
signature** (⟳ ÷ Σ ↑ ∂), **epistemic status** (Df/Pp/Pr/Conj/Heur/Cert ×
EXACT/CERTIFIED/DIAGNOSTIC/PHENO), and **derivation chain**. Zero free
parameters means every chain terminates in the integers — and you can watch:

```
$ python -m mtft.legend trace alpha_inverse
alpha_inverse  [Pr, CERTIFIED(3.5ppm)]
   └─ monster_order  [Pr, EXACT]  └─ integers  → It was always the integers.
   └─ genus_13       [Pr, EXACT]  └─ X0_143 └─ N_143 └─ integers  → ...
   └─ dim_S2_11      [Pr, EXACT]  └─ X0_143 └─ N_143 └─ integers  → ...

$ python -m mtft.legend status        # epistemic audit of the whole surface
$ python -m mtft.legend status EXACT  # just the bedrock you can build on
```

No other scientific package can ship `trace` honestly: everyone else's
chains end in a fitted constant.

## Machine Certificates (Tier 11)

The **Jacobian Conjecture counterexample** (July 2026) ships as a
self-certifying, dependency-free module — its own exact-arithmetic engine
re-derives the whole mechanism (constant Jacobian, S₃ monodromy, non-surjective
image) on every call. A counterexample cannot be faked:

```python
import mtft
cert = mtft.jc_verify_all(verbose=True)   # det DF = -2, degree 3, JC false n>=3
```

## Installation

```bash
pip install mtft            # fresh install
pip install --upgrade mtft  # already installed? pip won't upgrade unless asked
```

## Quick Start

```python
import mtft

# ── Modular forms on X₀(143) ──────────────────────────────
C = mtft.X0(143)                     # the modular curve, LMFDB-anchored
C.genus                              # 13
mtft.dedekind_eta(1j)                # 0.7682254223260566 = Γ(1/4)/(2π^(3/4))

# ── Gauge-Higgs via Hosotani holonomy ─────────────────────
h = mtft.HosotaniMTFT()
h.find_vacuum()                      # the Hosotani angle θ₀
h.gauge_masses()                     # m_W, m_Z, m_H from the holonomy

# ── Mass gap / stiffness ──────────────────────────────────
mtft.filtered_moment_identity(0.18174, N=3)   # exact to machine precision

# ── Falsifiability: the honest scorecard ──────────────────
mtft.falsify.honest_report()         # 23 pre-registered zero-parameter
                                     # predictions, ppm deviations, no cherry-picking

# ── Arithmetic computation (v0.7.0) ───────────────────────
from mtft.arithmetic_machine import decompose_turing_machine
decompose_turing_machine()           # a TM as a five-primitive AG object
```

CLI: `python -m mtft verify | report | tower | screen | info`

## Package Structure — 35 modules in 16 tiers

| Tier | Modules | What lives here |
|---|---|---|
| 0 | `constants` `arithmetic` | wₙ weights, X₀(143) structural constants |
| 1 | `modular` `modular_curve` `forms` `x0_143` | Modular geometry, LMFDB-validated newforms |
| 2 | `hosotani` `tano_metric` | Gauge-Higgs unification, the Tano metric |
| 3 | `particles` `koide` `decay` | Spectrum phenomenology |
| 4 | `lattice` `burning_ship` | Lattice gauge, Burning Ship fermions |
| 5 | `dimensional_bridge` `cosmology` `dark_sector` | Bridges & the dark sector |
| 5b | `falsify` | 23 pre-registered predictions, honest scorecard |
| 5c | `tower` | SU(N) landscape |
| 5d | `riemann` | Explicit formula, corrected RH diagnostic, **Speiser–Hadamard lab** |
| 5e | `lhcb_analysis` | ROOT-bridge confrontation (uproot) |
| 6 | `quantum` | Topological qudits |
| 7 | `crypto` `monster_hash` | SL(2,ℤ)-sponge hashing |
| 8 | `info_geometry` `jacobian` | Fisher–Rao curvature, Jacobian engine |
| 9 | `arithmetic_machine` `arithmetic_wick` `busy_beaver` `music` `viz` | Computation as arithmetic, Wick bridge, sonification |
| 10 | `critical_ensemble` | Li coefficients λₙ — three cross-certified methods |
| 11 | `jc_counterexample` `estimator_standards` `legend` | Machine certificates, A.7 estimator standards, the Legend |

## Key Identities

```
α⁻¹        = ln|M| + genus − 1/11 + O(10⁻⁴)          Monster ↔ fine structure
μ_N(y)     = (1/4π²)[T″(y) − Re T″(y − i/N)]          mass gap = twisted-curvature gap (exact)
Z_D(β)     = −ζ(β)·ζ′(β+1)                            Dirichlet ensemble closed form (exact)
∂²log(−ζ′) = 2/(s−1)² − Σ_ρ′ (s−ρ′)⁻²                Hadamard over ζ′ zeros (certified < 1e−5, s ∈ [3,10])
λ₁         = 1 + γ/2 − ½ln(4π)                        critical-ensemble anchor (exact)
```

## Visuals

`viz/` ships the gallery: `hero_stiffness.png` + `stiffness_navigator.html`
(regenerate with `python viz/make_hero.py`), plus React components —
`X0_143_BurningMandelbrot.jsx`, `MTFT_HyperbolicTiling.jsx`,
`MTFT_MonsterFingerprint.jsx`, `mtft_enneper.jsx`.

## Testing

```bash
pytest          # 401 tests: LMFDB anchors, exact identities, audit regressions
```

Every numerical claim above is a test. The suite gates PyPI publication —
nothing ships if anything fails.

## Citation

See [`CITATION.cff`](CITATION.cff). GitHub's "Cite this repository" button
does the formatting.

## License

MIT — see [LICENSE](LICENSE).
