# MTFT — Modular Time Field Theory

[!\[Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[!\[License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python library implementing the core data structures and computational tools for **Modular Time Field Theory (MTFT)**.

## Overview

MTFT proposes a modular time field **τ(x) = t\_R / t\_U** mapping spacetime to the upper half-plane **ℍ**, unifying dark matter, particle masses, and cosmic structure through:

* **Modular symmetry** — SL(2,ℤ) acting on τ
* **Gauge-Higgs unification** — Higgs = A\_τ holonomy via the Hosotani mechanism
* **τ-vortex dark matter** — logarithmic τ-profiles give flat rotation curves with no dark particles
* **Information geometry bridge** — Fisher-Rao curvature of the logistic map connects chaos → geometry → dynamics

### The Three-Layer Architecture

```
Counting   →  Logistic map orbits, figurate-number degeneracies
Geometry   →  Fisher-Rao metric, Ricci curvature R\_core
Dynamics   →  SM masses, dark matter halos, decay rates, cosmology
```

Connected by the **spectral determinant identity**:

```
det(1 − q^{1/m} P\_m) = η(τ)^{−1/m} θ₃(0,τ)^{1/m}
```

linking Fredholm determinants (chaos/RG) to CFT partition functions.

## Installation

```bash
pip install -e .              # minimal (numpy only)
pip install -e ".\[full]"      # with scipy + matplotlib
pip install -e ".\[dev]"       # with pytest + ruff
```

## Quick Start

```python
import mtft

# ── Modular forms ─────────────────────────────────────
tau = 0.1 + 1.5j
eta = mtft.dedekind\_eta(tau)
j   = mtft.forms.j\_invariant(tau)
print(f"η(τ) = {eta:.6f}")
print(f"j(τ) = {j:.2f}")

# Verify the spectral determinant identity
result = mtft.forms.verify\_spectral\_identity(tau, m=2)
print(f"Spectral identity relative error: {result\['relative\_error']:.2e}")

# ── Hosotani mechanism ────────────────────────────────
hp = mtft.HosotaniPotential(fermion\_fraction=0.4, kappa\_ew=0.05)
theta0 = hp.find\_vacuum()
masses = hp.gauge\_masses()
print(f"Vacuum θ₀ = {theta0:.4f}")
print(f"m\_W = {masses\['m\_W']:.2f} GeV  (PDG: 80.37)")
print(f"m\_Z = {masses\['m\_Z']:.2f} GeV  (PDG: 91.19)")

# ── Particle spectrum ─────────────────────────────────
sm = mtft.StandardModel()
top = sm.by\_name("Top")
print(f"Top quark: m = {top.mass\_GeV} GeV, κ = {top.kappa}")

# Full κ-hierarchy
for name, kappa in sm.kappa\_hierarchy():
    print(f"  {name:15s}  κ = {kappa:.2e}")

# ── Dark sector ───────────────────────────────────────
import numpy as np
halo = mtft.TauVortexHalo(A=1e20, r0=1e30)
r = np.logspace(31, 35, 100)
v = halo.v\_circular(r)
print(f"v\_∞ = {halo.v\_infinity:.4e} (flat rotation velocity)")

# ── Information geometry ──────────────────────────────
R = mtft.info\_geometry.R\_core()
print(f"R\_core (Feigenbaum) = {R:.4f}")

# ── Cosmology ─────────────────────────────────────────
cosmo = mtft.FriedmannMTFT(Omega\_tau=0.25)
hist = cosmo.expansion\_history()
```

## Package Structure

```
mtft/
├── constants.py       # PDG masses, SM parameters, physical constants
├── modular.py         # τ-field, SL(2,ℤ), hyperbolic geometry
├── forms.py           # Dedekind η, Jacobi θ₃, Eisenstein, spectral det
├── hosotani.py        # Effective potential, vacuum finder, EWSB
├── particles.py       # SM particle database with κ-couplings
├── dark\_sector.py     # τ-vortex halos, rotation curves, Tully-Fisher
├── info\_geometry.py   # Fisher-Rao metric, Ricci curvature, logistic map
└── cosmology.py       # Modified Friedmann, perturbations, G(t) oscillation
```

## Key Equations

|Component|Equation|
|-|-|
|Modular time field|τ(x) = t\_R(x)/t\_U(x) ∈ ℍ|
|Hyperbolic metric|ds² = (dx² + dy²)/y², K = −1|
|W mass|m\_W = κ\_EW \|sin θ\_H\| / (2R\_τ)|
|Fermion mass|m\_f = κ\_f \|sin θ₀\| / R\_τ|
|τ-vortex density|ρ\_τ(r) = A²/(2r²)|
|Flat rotation|v²\_∞ = 2πGA²|
|Tully-Fisher|v⁴ ∝ M\_BH ∝ M\_baryonic|
|R\_core bridge|m\_τ² = c\_m \|R\_core\||

## Testing

```bash
pytest tests/ -v
```

## License

MIT

## Citation

If you use this package in research, please cite:

```bibtex
@software{mtft2026,
  title  = {MTFT: Modular Time Field Theory Python Package},
  author = {Roger},
  year   = {2026},
  url    = {https://github.com/Kaizoku-Ronin/mtft}
}
```

