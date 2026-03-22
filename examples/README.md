# MTFT Examples

Install: `pip install mtft`

Then run any example:

```bash
python examples/01_quick_start.py
```

## Examples

| # | File | What it shows |
|---|------|---------------|
| 01 | `01_quick_start.py` | Core idea in 30 lines — gauge couplings from the integers |
| 02 | `02_arithmetic_confinement.py` | Weights w_n, torque convergence, SU(3) confinement lock |
| 03 | `03_lepton_masses.py` | Koide theorem, Tano mass formula, X_0(143) structure |
| 04 | `04_lattice_monte_carlo.py` | SU(3) lattice simulation with MTFT action, Wilson loops |
| 05 | `05_dark_sector.py` | Tau-vortex halos, flat rotation curves, Tully-Fisher |
| 06 | `06_modular_forms.py` | Dedekind eta, j-invariant, dimensional bridge formula |
| 07 | `07_burning_ship.py` | Fermion vacuum, anisotropic Feigenbaum, three generations |
| 08 | `08_quantum_computing.py` | Holonomy gates, topological qudits, error correction |
| 09 | `09_cryptography.py` | Arithmetic hash, PRNG, SL(2,Z) key exchange, lattice LWE |
| 10 | `10_full_scorecard.py` | All 22 predictions vs PDG with formulas and errors |
| 11 | `11_info_geometry_particles.py` | Fisher-Rao metric, Hosotani mechanism, SM catalog |
| 12 | `12_visualization.py` | Matplotlib plots (needs `pip install mtft[full]`) |
| 13 | `13_decay.py` | Nuclear decay in modular time, neutron lifetime anomaly |

## Quick verification

```python
from mtft.verify import print_report
print_report()
```

## Requirements

- `mtft` (core): `pip install mtft` — numpy only
- `mtft[full]`: adds scipy + matplotlib for plotting
- `mtft[viz]`: adds plotly for interactive 3D
