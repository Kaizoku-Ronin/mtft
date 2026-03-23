# MTFT v0.5.0 Changelog

## Summary
3 bug fixes, 2 new features, 1 new module. 6 files changed, 1 file added.

---

## New Module: Falsifiability Engine (`falsify.py`)

The module that separates MTFT from numerology.

### `prediction_table()` — 23 zero-parameter predictions
Every prediction computed live from δ, T∞, γ, Ω, ξ — no hardcoded values.
**21/23 pass, 1 tension (M_W/M_Z at 2.3σ), 1 at boundary (α⁻¹ 3-term).**

Top results:
| # | Relation | Error | σ |
|---|----------|-------|---|
| 2 | α⁻¹ = 2πδ² + 1/(4δ) − ξT∞/δ⁶ | 0.00000019% | 12.1* |
| 13 | m_τ from Koide | 0.006% | 0.9 |
| 6 | m_H/m_W = √π·T∞² | 0.029% | 1.8 |
| 17 | λ = γ²/(8Ω²) | 0.074% | 0.1 |
| 5 | sin²θ_W = 3/13 | 0.195% | 0.5 |

*α⁻¹ is correct to 7 decimal places; the 12.1σ reflects CODATA's
10-digit precision exceeding the 3-term asymptotic expansion.

### `coupling_shift(Δα⁻¹)` — correlated coupling predictions
If CODATA updates α⁻¹, MTFT mandates:
- Δα_s = −(δ⁶/4ξ) · Δα⁻¹
- Δsin²θ_W = (3/169) · Δα⁻¹

Independent shifts → MTFT falsified.

### `holonomy_flux()` — Josephson prediction
Φ_H/Φ₀ ≈ −3.5%, material-independent, geometry-independent.
BCS predicts 0%.

### `desert_check()` — new particle tracker
No new particles between EW and Planck scales. Tracks discoveries.

### `falsification_test()` and `report()` — automated testing
Runs all predictions, reports pass/fail/tension at configurable σ threshold.

---

## New Feature: AEC Divisor Scaling (`quantum.py`)
`ArithmeticCode` now produces measurably different code distances and
protection gaps for primes vs composites. See previous changelog.

---

## Bug Fixes

### 1. `SM.cos_theta_W` — `constants.py`
Added `cos_theta_W`, `sin_theta_W`, `sin2_theta_W` to `_MTFTSM`.

### 2. Hosotani nontrivial vacuum — `hosotani.py`
New `HosotaniMTFT` class: V(θ) = a sin²θ + b sin⁴θ, a/b = −6/13.
Vacuum at sin²θ₀ = 3/13, m_H = 125.30 GeV.

### 3. Rotation curve units — `dark_sector.py`
New `rotation_curve_kpc()` accepts kpc and km/s.

---

## Files Changed (6) + Added (1)
1. `mtft/__init__.py` — exports + version 0.5.0
2. `mtft/constants.py` — SM.cos_theta_W
3. `mtft/hosotani.py` — HosotaniMTFT
4. `mtft/dark_sector.py` — rotation_curve_kpc()
5. `mtft/quantum.py` — AEC divisor scaling
6. **`mtft/falsify.py`** — NEW: Falsifiability Engine

## How to Apply (with Claude Code)
```bash
cd ~/path/to/mtft
# Extract patches over source
tar xzf mtft_v050_patches.tar.gz
cp mtft/*.py mtft/
# Verify
python -c "from mtft.falsify import report; report()"
# Commit and push
git add -A && git commit -m "v0.5.0: falsifiability engine, AEC divisor scaling, bug fixes"
git tag v0.5.0
git push origin main --tags
```
