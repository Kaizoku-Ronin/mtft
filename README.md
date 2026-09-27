# MTFT — Modular Time Field Theory

[![PyPI](https://img.shields.io/pypi/v/mtft.svg?color=blue)](https://pypi.org/project/mtft/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

MTFT is a research toolkit for arithmetic weights, modular geometry, spectral
dynamics and tests of proposed connections to physics. Its central objects are

\[
w_n=\sum_{d\mid n}\frac{\log d}{d},\qquad X_0(143),\qquad 143=11\times13.
\]

The software has grown from stiffness and phenomenology tools into an exact
arithmetic and modular-surface laboratory: Hecke operators, canonical rings,
periods and Hodge structures, spin characteristics, field dynamics, Yukawa
tensors, anomaly tests and model-selection gates. It also includes experimental
cryptography, arithmetic computation, sonification and visualization tools.

**Research status:** integer-derived formulas, numerical certificates and
physical hypotheses have different evidence requirements. MTFT does not yet
supply a validated Standard Model or a general simulator of matter. The new
interaction atlas provides external SM targets against which future MTFT
derivations can be checked. Its reference parameters are not MTFT predictions.

**Start here:** [Capability index](docs/CAPABILITIES.md) ·
[SM interaction atlas](docs/SM/INTERACTION_ATLAS.md) ·
[Research and correction registers](docs/SM/) ·
[Release history](docs/changelog/README.md) · [Examples](examples/README.md)

## Install

```bash
python -m pip install --upgrade mtft
```

For development, including changes not yet released on PyPI:

```bash
git clone https://github.com/Kaizoku-Ronin/mtft.git
cd mtft
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell instead:
# .venv\Scripts\Activate.ps1
python -m pip install -e ".[dev,full]"
```

| Dependency set | Purpose |
|---|---|
| Base package | NumPy, SymPy and mpmath arithmetic/numerical engines |
| `mtft[full]` | SciPy solvers/quadrature and Matplotlib plots |
| `mtft[viz]` | Matplotlib and Plotly visualization helpers |
| `mtft[dev]` | pytest, coverage and Ruff |
| `mtft[lhc]` | uproot and Awkward for the optional ROOT/LHCb analysis bridge |
| PARI/GP, installed separately | Recompute GP-dependent modular-symbol and period data |
| MadGraph5_aMC@NLO, installed separately | Optional SM process-diagram/code generation |

The bundled SM atlas works offline with the base install. Many geometry engines
can use frozen data without PARI/GP; regenerating those records is separate work.

## Quick start

```python
import mtft
from mtft.interactions import load_catalog, select_vertices

C = mtft.X0(143)
print(C.genus)                         # 13
print(mtft.dedekind_eta(1j))
print(mtft.filtered_moment_identity(0.18174, N=3))

# Inspect the SM reference without computing an MTFT amplitude.
catalog = load_catalog()
print(len(select_vertices(catalog, sector="qcd")))  # 8 for this model

# Existing phenomenology comparison, with its own assumptions and inputs:
print(mtft.falsify.honest_report())
```

```bash
python -m mtft info
python -m mtft report
python -m mtft.legend search Hodge
python -m mtft.legend card surface_higgs_slope_level
python -m mtft.interactions render -o atlas-output/index.html
```

Open `atlas-output/index.html` in a browser. Search vertices, inspect diagrams
and symbolic tensors, filter sectors, follow related MTFT research and export
the reference or evidence ledger. The atlas contains **153 reference vertices**
from a pinned stock model and explicitly lists its omissions. It does not
represent every possible SM process diagram. See the
[atlas guide](docs/SM/INTERACTION_ATLAS.md) for evidence statuses, import limits,
scope, licenses and MadGraph integration.

A [prebuilt offline atlas](viz/sm_interaction_atlas.html) is also included:
download the HTML file and open it directly, without installing MTFT.

## What the package can do

The table groups the complete module surface by task. The
[generated module index](docs/CAPABILITIES.md) links every public module to its
source and entry points, covering both older tools and current research code.

| Work area | Capabilities and principal modules |
|---|---|
| Arithmetic foundations | Constants, divisor-log weights, filtered moments and combinatorial ancestry: `constants`, `arithmetic`, `combinatorial` |
| Modular functions and curves | SL(2,ℤ), eta/modular forms, X₀(N), newform data and X₀(143) anchors: `modular`, `modular_curve`, `forms`, `x0_143`, `levels`, `coset_reps` |
| Riemann and ensembles | Explicit formula, ζ′/Speiser diagnostics, Dirichlet curvature, Li coefficients and compressed Weil forms: `riemann`, `arithmetic_wick`, `critical_ensemble`, `weil` |
| Weight moments and information geometry | Closed-form moments, statistical metrics, curvature and stiffness landscapes: `moments`, `curvature`, `info_geometry`, `jacobian`, `tower` |
| Spectral dynamics | Marked primon gas, KMS diagnostics, internal/coupled chains, gap extraction and exceptional-point studies: `marked_gas`, `marked_gap`, `chain`, `coupled`, `ep`, `expansion`, `exception_spectrum` |
| Mellin and L-function channels | Bulk/skeleton peel, Dirichlet channels, GL(2) diagnostics and analytic benchmarks: `peel`, `lchannels`, `gl2_peel`, `hardy_ramanujan` |
| Hecke and integral arithmetic | Manin symbols, Hecke blocks, Eisenstein congruences, cuspidal torsion, codifferents and integer lattices: `hecke`, `eisenstein`, `cuspidal`, `codifferent`, `integral_lattice`, `quadratic_forms` |
| Canonical geometry | Canonical ideals, quadrics, Atkin–Lehner descent, mod-p and Petri gates: [`canonical`](src/mtft/canonical/) |
| Periods and Hodge geometry | Period matrices, symplectic frames, polarization, Bergman density, Hamiltonian channels and stability: [`periods`](src/mtft/periods/), `hodge_polarization`, `al_morphology` |
| Homology, spin and theta | Integral homology, affine spin actions, Arf parity, theta functions and finite Kakeya diagnostics: [`homology`](src/mtft/homology/), [`thetachar`](src/mtft/thetachar/), `thetafun`, `kakeya` |
| Discrete and boundary geometry | Origami/dimers, insertion calculus, boundary diagnostics and finite graph tools: [`origami`](src/mtft/origami/), [`boundary`](src/mtft/boundary/) |
| Modular-surface geometry | Manin meshes, cycle lattices, exact transport, Hodge frames, marked/oldform sectors, CRT dessins, spin circles and Hopf bookkeeping: [`surface`](src/mtft/surface/) |
| Fields on the surface | Yang–Mills partition functions, line operators, Ising models, Hamiltonian evolution, FEM Laplace/Bochner spectra, Petersson and HYM normalization: `surface.gauge`, `.ising`, `.dynamics`, `.spectral`, `.magnetic`, `.petersson`, `.hym` |
| Fluxes and Yukawa tensors | Spin/CM data, Riemann–Roch spaces, canonical multiplication tensors, flux scans and condensation: `surface.arithspin`, `.rrspace`, `.yukawa`, `.smflux`, `.condensation` |
| Parent-theory gates | Parent records, anomaly polynomials, charge lattices, discrete/bordism tests, tensor factorization, mode/chirality gates and compactification: [`research`](src/mtft/research/) |
| Product-surface research | Vector Higgs operators, X₀(143)×E slopes/indices, Künneth selection rules and torus theta factors: `research.vector_higgs`, `.product_surface`, `.yukawa_triangle`, `.theta_torus`, `.unified_parent`, `.involutions` |
| SM interaction reference | Pinned UFO vertices, tensor/parameter data, SVG/HTML diagrams and an evidence ledger: [`interactions`](docs/SM/INTERACTION_ATLAS.md) |
| Physics phenomenology | Particle data/embeddings, Hosotani mechanism, Koide relations, decay, dimensional bridges, cosmology/dark-sector models and comparison reports: `particles`, `hosotani`, `koide`, `decay`, `dimensional_bridge`, `cosmology`, `dark_sector`, `falsify`, `verify` |
| Lattice, algebra and materials | Lattice gauge experiments, Burning Ship models, Lie-algebra closure/SVD gates and material-property diagnostics: `lattice`, `burning_ship`, `liealg`, `tano_metric` |
| Quantum and experimental cryptography | Qudits/holonomy gates, arithmetic codes, SL(2,ℤ)-sponge experiments and Jacobian orders: `quantum`, `crypto`, `monster_hash`; research implementations, not a security certification |
| Arithmetic computation | Primitive signatures, Turing-machine decompositions, halting diagnostics and the package's exact-polynomial certificate: `arithmetic_machine`, `busy_beaver`, `jc_counterexample` |
| Provenance and usability | Epistemic metadata, constant ledgers, estimator standards, PARI runner, plots and sonification: `legend`, `ledger`, `ledger_peel`, `estimator_standards`, `gprun`, `viz`, `music` |
| Experimental data bridge | Optional ROOT-file/LHCb analysis via uproot and Awkward: `lhcb_analysis` |

The original three-ensemble route remains available:

| Assembly | Object | Tools |
|---|---|---|
| Laplace | Weighted theta sums and filtered curvature | `weighted_theta`, `filtered_moment_identity`, `peel` |
| Dirichlet | −ζ(β)ζ′(β+1) and its log curvature | `dirichlet_curvature`, `hadamard_zetaprime_check`, `marked_gas` |
| Critical | Li coefficients from log ξ | `li_criterion_report`, `critical_ensemble` |
| Weil extension | Compressed explicit-formula quadratic form | `weil` |

Finite scans and truncated tests remain diagnostics unless their documented
theorem or error bound supports a stronger conclusion.

## Read the evidence and corrections

The Legend is a curated map, with tags `Df`, `Pp`, `Pr`, `Conj`, `Heur`, `Cert`
and labels such as `EXACT`, `CERTIFIED`, `DIAGNOSTIC`, `PHENO` and `GIVEN`.
It does not yet contain a hand-written entry for every callable.

```bash
python -m mtft.legend status
python -m mtft.legend trace alpha_inverse
python -m mtft.legend trace sm_interaction_atlas
```

The atlas trace terminates at **external SM reference inputs**. Its data is
not presented as an arithmetic derivation. A passed software test verifies its
stated calculation, not an unrestricted claim of physical validity.

For the current v0.33.x research line (v0.33.1 adds the SM interaction atlas and the review corrections listed in [the v0.33.1 changelog](docs/changelog/CHANGELOG_v0331.md)), start with:

- [Compendium implementation notes](docs/SM/V0330_COMPENDIUM_NOTES.md): exact
  inputs, operator conventions, rank qualifications and independent routes.
- [CC-33 tensor coefficient correction](docs/SM/CC33_TENSOR_COEFFICIENT.md):
  the corrected six-dimensional gravitational-anomaly obstruction.
- [CC-26 normalization retraction](docs/SM/CC26_NORMALISATION_RETRACTION.md):
  an earlier Yukawa normalization claim and its corrected interpretation.
- [Wave C gates](docs/SM/WAVE_C_GATES.md),
  [chirality register](docs/SM/R2C01_CHIRALITY_REGISTER.md) and
  [product-surface register](docs/SM/R2C07_PRODUCT_SURFACE_REGISTER.md).

Historical results and frozen handoffs stay in `studies/`. They may contain
superseded hypotheses and release-specific scripts; their original bytes and
correction trail matter. New atlas mappings require explicit evidence.

## Visuals and long calculations

![Stiffness landscape](https://raw.githubusercontent.com/Kaizoku-Ronin/mtft/main/viz/hero_stiffness.png)

The [visualization gallery](viz/README.md) includes the stiffness navigator,
modular-surface and tiling components, arithmetic fingerprints and drawn-loop
tools. The graph-clock wavepacket animation is a model evolution on the
X₀(143) dual graph, not a simulation of experimentally established matter:

![Wavepacket on the X₀(143) graph](https://raw.githubusercontent.com/Kaizoku-Ronin/mtft/main/viz/wavepacket_X0143.gif)

For long PARI/GP work:

```bash
python -m mtft.gprun
```

The local browser runner locates `gp` (or the executable specified by `MTFT_GP`),
freezes the input script, stamps its hash and GP version, streams output to
disk, and records exit status. See [PARI scripts](scripts/pari/README.md) and
the [study guide](studies/README.md). Hours-long studies are separate from
routine tests; choose a study's documented parameters and preserve its logs.

## Test and maintain

```bash
# Atlas and its integration:
python -m pytest tests/test_interactions.py tests/test_tier11_and_legend.py -q
# Routine suite excluding tests marked slow (some unmarked tests still take minutes):
python -m pytest tests/ -m "not slow" -q
# The existing release gate:
python -m pytest tests/ -q
# Update or check the module index:
python scripts/build_capability_index.py
python scripts/build_capability_index.py --check
```

Optional dependencies and PARI/GP affect skips. Some research tests also use
their own opt-ins, such as `MTFT_SLOW=1` for the full theta census. Follow the
test/study's instructions rather than treating one command as every archived
experiment. No fixed test-count badge is maintained here.

Release notes are organized in [`docs/changelog/`](docs/changelog/README.md),
with a separate [unreleased record](docs/changelog/UNRELEASED.md). Publishing
uses the existing [release workflow](PUBLISHING.md); the package version,
`__init__` version and `CITATION.cff` must agree before a release.

## Citation and license

See [`CITATION.cff`](CITATION.cff) and [`LICENSE`](LICENSE) (MIT).
The external MadGraph reference retains its
[upstream license](src/mtft/interactions/_data/MADGRAPH_LICENSE.txt).
