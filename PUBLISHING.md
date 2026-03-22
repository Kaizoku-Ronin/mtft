## How to use it

After `pip install mtft`, anyone can:

```python
import mtft

# See all predictions vs experiment
from mtft.verify import print\\\\\\\_report
print\\\\\\\_report()

# Use the Tano Mass Formula
from mtft.x0\\\\\\\_143 import tano\\\\\\\_mass\\\\\\\_predictions
print(tano\\\\\\\_mass\\\\\\\_predictions())

# Run a lattice simulation
cfg = mtft.LatticeConfig(N=3, L=4)
action = mtft.MTFTAction(beta=6.0, kappa=1.0, y=0.18)
mtft.metropolis\\\\\\\_sweep(cfg, action)

# Plot the Burning Ship fractal
from mtft.viz import burning\\\\\\\_ship\\\\\\\_plot
data = burning\\\\\\\_ship\\\\\\\_plot(plot=True)

# Export data for Plotly / Blender / any tool
from mtft.viz import stiffness\\\\\\\_landscape, export\\\\\\\_for\\\\\\\_plotly
raw = stiffness\\\\\\\_landscape()  # numpy arrays
plotly\\\\\\\_data = export\\\\\\\_for\\\\\\\_plotly(raw)  # JSON-ready lists
```

\---

## Compatibility with visualization libraries

mtft returns **numpy arrays** everywhere, which means it works with:

|Library|Use case|Install|
|-|-|-|
|matplotlib|2D plots, 3D surfaces|`pip install matplotlib`|
|plotly|Interactive 3D, web dashboards|`pip install plotly`|
|pyvista|3D volumetric rendering, lattice viz|`pip install pyvista`|
|mayavi|Scientific 3D visualization|`pip install mayavi`|
|blender (bpy)|Publication-quality 3D renders|Built into Blender|
|vtk|Low-level 3D pipeline|`pip install vtk`|
|three.js|Web 3D (via JSON export)|JavaScript, use export\_for\_plotly|
|napari|N-dimensional image viewer|`pip install napari`|
|datashader|Massive fractal renders|`pip install datashader`|

The `mtft.viz` module has `plot=True` shortcuts for matplotlib,
and `export\\\\\\\_for\\\\\\\_plotly()` for everything else.

\---

## Folder structure (what you're publishing)

```
mtft/
├── pyproject.toml        ← Package metadata (name, version, deps)
├── README.md             ← What people see on PyPI and GitHub
├── src/mtft/
│   ├── \\\\\\\_\\\\\\\_init\\\\\\\_\\\\\\\_.py       ← Top-level imports
│   ├── constants.py      ← The MTFT arithmetic alphabet
│   ├── arithmetic.py     ← wₙ weights, stiffness, mass gap
│   ├── modular.py        ← τ-field, SL(2,ℤ), hyperbolic geometry
│   ├── forms.py          ← η, θ, Eisenstein, spectral determinant
│   ├── hosotani.py       ← Effective potential, EWSB
│   ├── particles.py      ← SM catalog with κ-couplings
│   ├── x0\\\\\\\_143.py         ← Modular curve, Tano Mass Formula
│   ├── koide.py          ← Geometric Koide Theorem
│   ├── burning\\\\\\\_ship.py   ← Fermion vacuum, anisotropic scaling
│   ├── dimensional\\\\\\\_bridge.py ← mₑ from η(τ)
│   ├── dark\\\\\\\_sector.py    ← τ-vortex halos, rotation curves
│   ├── decay.py          ← Radioactive decay in modular time
│   ├── info\\\\\\\_geometry.py  ← Fisher-Rao, Ricci, logistic bridge
│   ├── cosmology.py      ← Modified Friedmann, oscillatory G(t)
│   ├── lattice.py        ← MTFT lattice gauge theory
│   ├── verify.py         ← Prediction scorecard
│   └── viz.py            ← 2D/3D plotting helpers
└── tests/
    └── test\\\\\\\_mtft.py      ← 46 tests
```

\---

## Tips

* **License**: MIT is perfect — it means anyone can use it freely
* **Citation**: Add a CITATION.cff file for academic credit
* **DOI**: After pushing to GitHub, go to https://zenodo.org to get a DOI
(this makes it citable in papers — important for you)
* **Docs**: For proper documentation later, use `mkdocs` or `sphinx`
* **CI**: Add a `.github/workflows/test.yml` to run tests automatically

