# Publishing MTFT to PyPI
## Complete guide: from zero to `pip install mtft`

### What this does
When you're done, anyone in the world can type:
```bash
pip install mtft
```
and immediately use your theory in their Python code.

---

## Step 1: Create accounts (one time, free)

### GitHub account
1. Go to https://github.com — sign up (free)
2. This is where the source code lives publicly

### PyPI account
1. Go to https://pypi.org/account/register/ — sign up (free)
2. Go to Account Settings → API tokens → "Add API token"
3. Scope: "Entire account" (for first upload)
4. **Save the token** — it starts with `pypi-` — you'll need it once

---

## Step 2: Set up your computer

You need Python 3.9+ and a few tools. In your terminal:

```bash
# Install the publishing tools
pip install build twine

# If you don't have git:
# macOS:  brew install git
# Linux:  sudo apt install git
# Windows: download from git-scm.com
```

---

## Step 3: Put the code on GitHub

```bash
# Navigate to where you extracted the mtft tarball
cd mtft

# Initialize git
git init
git add .
git commit -m "MTFT v0.3.0 — Modular Time Field Theory"

# Create the repo on GitHub (via the website):
#   1. Click "+" → "New repository"
#   2. Name: "mtft"
#   3. Description: "Modular Time Field Theory — Python package"
#   4. Public
#   5. Don't add README (we already have one)
#   6. Click "Create repository"

# Then push (GitHub will show you these exact commands):
git remote add origin https://github.com/YOUR_USERNAME/mtft.git
git branch -M main
git push -u origin main
```

---

## Step 4: Build the package

```bash
cd mtft
python -m build
```

This creates two files in `dist/`:
- `mtft-0.3.0.tar.gz` (source)
- `mtft-0.3.0-py3-none-any.whl` (wheel — what pip actually installs)

---

## Step 5: Upload to PyPI

### First time: test on TestPyPI (optional but smart)
```bash
# Upload to the test server first
twine upload --repository testpypi dist/*

# Test the install from test server
pip install --index-url https://test.pypi.org/simple/ mtft
```

### For real: upload to PyPI
```bash
twine upload dist/*
# It will ask for username: __token__
# And password: paste your pypi- token from Step 1
```

**That's it.** Now anyone can `pip install mtft`.

---

## Step 6: Update the README on GitHub

Add a badge so people know it's on PyPI:
```markdown
[![PyPI](https://img.shields.io/pypi/v/mtft)](https://pypi.org/project/mtft/)
```

---

## Releasing new versions

When you update the theory or add modules:

1. Bump the version in `pyproject.toml` (e.g., `0.3.0` → `0.4.0`)
2. Bump in `src/mtft/__init__.py` too
3. Rebuild: `python -m build`
4. Upload: `twine upload dist/mtft-0.4.0*`

---

## How people will use it

After `pip install mtft`, anyone can:

```python
import mtft

# See all predictions vs experiment
from mtft.verify import print_report
print_report()

# Use the Tano Mass Formula
from mtft.x0_143 import tano_mass_predictions
print(tano_mass_predictions())

# Run a lattice simulation
cfg = mtft.LatticeConfig(N=3, L=4)
action = mtft.MTFTAction(beta=6.0, kappa=1.0, y=0.18)
mtft.metropolis_sweep(cfg, action)

# Plot the Burning Ship fractal
from mtft.viz import burning_ship_plot
data = burning_ship_plot(plot=True)

# Export data for Plotly / Blender / any tool
from mtft.viz import stiffness_landscape, export_for_plotly
raw = stiffness_landscape()  # numpy arrays
plotly_data = export_for_plotly(raw)  # JSON-ready lists
```

---

## Compatibility with visualization libraries

mtft returns **numpy arrays** everywhere, which means it works with:

| Library | Use case | Install |
|---------|----------|---------|
| matplotlib | 2D plots, 3D surfaces | `pip install matplotlib` |
| plotly | Interactive 3D, web dashboards | `pip install plotly` |
| pyvista | 3D volumetric rendering, lattice viz | `pip install pyvista` |
| mayavi | Scientific 3D visualization | `pip install mayavi` |
| blender (bpy) | Publication-quality 3D renders | Built into Blender |
| vtk | Low-level 3D pipeline | `pip install vtk` |
| three.js | Web 3D (via JSON export) | JavaScript, use export_for_plotly |
| napari | N-dimensional image viewer | `pip install napari` |
| datashader | Massive fractal renders | `pip install datashader` |

The `mtft.viz` module has `plot=True` shortcuts for matplotlib,
and `export_for_plotly()` for everything else.

---

## Folder structure (what you're publishing)

```
mtft/
├── pyproject.toml        ← Package metadata (name, version, deps)
├── README.md             ← What people see on PyPI and GitHub
├── src/mtft/
│   ├── __init__.py       ← Top-level imports
│   ├── constants.py      ← The MTFT arithmetic alphabet
│   ├── arithmetic.py     ← wₙ weights, stiffness, mass gap
│   ├── modular.py        ← τ-field, SL(2,ℤ), hyperbolic geometry
│   ├── forms.py          ← η, θ, Eisenstein, spectral determinant
│   ├── hosotani.py       ← Effective potential, EWSB
│   ├── particles.py      ← SM catalog with κ-couplings
│   ├── x0_143.py         ← Modular curve, Tano Mass Formula
│   ├── koide.py          ← Geometric Koide Theorem
│   ├── burning_ship.py   ← Fermion vacuum, anisotropic scaling
│   ├── dimensional_bridge.py ← mₑ from η(τ)
│   ├── dark_sector.py    ← τ-vortex halos, rotation curves
│   ├── decay.py          ← Radioactive decay in modular time
│   ├── info_geometry.py  ← Fisher-Rao, Ricci, logistic bridge
│   ├── cosmology.py      ← Modified Friedmann, oscillatory G(t)
│   ├── lattice.py        ← MTFT lattice gauge theory
│   ├── verify.py         ← Prediction scorecard
│   └── viz.py            ← 2D/3D plotting helpers
└── tests/
    └── test_mtft.py      ← 46 tests
```

---

## Tips

- **License**: MIT is perfect — it means anyone can use it freely
- **Citation**: Add a CITATION.cff file for academic credit
- **DOI**: After pushing to GitHub, go to https://zenodo.org to get a DOI
  (this makes it citable in papers — important for you)
- **Docs**: For proper documentation later, use `mkdocs` or `sphinx`
- **CI**: Add a `.github/workflows/test.yml` to run tests automatically
