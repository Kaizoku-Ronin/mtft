# MTFT boundary experiments v0.1.0

Boundary observations, sensor placement, and a reversed boundary bond on the
56-spin dual graph supplied by MTFT 0.26.0 at level 143. The chosen interior has
14 spins, the boundary 8, and the exterior 34.

Read `MTFT_Boundary_Experiment_Report.md` for the model, results, and limitations.
The three `Boundary_*.png` figures have matching vector SVG copies.

## Reproduce

The supplied MTFT source archive is unchanged. Python 3.12.13 was used. The
requirements file records the exact numerical library versions used here.
Node.js is optional and only used to cross-check the JavaScript engine.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install ./mtft-0.26.0.tar.gz
python boundary_model.py
python run_boundary_experiments.py
python export_browser_data.py
python verify_boundary.py
python make_boundary_figures.py
```

For an existing environment with these dependencies, unpacking the MTFT archive
and putting its `src` directory on `PYTHONPATH` also works. The model evaluations
do not use a network, external service, GPU, PARI/GP, random seed, or sampling.

`boundary_model.py` constructs exact integer coefficients, checks all 256
exterior fibers, and compares the reconstructed full density of states against
independent whole-graph elimination. `run_boundary_experiments.py` evaluates
48 temperatures × 256 sensor subsets × 2 bond models and runs an independent
full-interior check plus a small-graph separator/shortcut control. The scripts
raise an assertion on a failed gate. Runtime fields are machine dependent;
numerical values can vary in their last few digits across platforms.

## Query an observation

```sh
python boundary_cli.py --beta .64 --observe 6=+ 18=- 24=+
python boundary_cli.py --beta .64 --model one_reversed_bond --observe all=+
python boundary_cli.py --beta .64 --best-sensors 2
```

Boundary vertex labels are `[6,18,22,24,36,40,50,53]`. Unlisted vertices are
unobserved. These are passive observations from equilibrium, with the exterior
summed out. They are not interventions that choose boundary patterns uniformly.
The supported inverse-temperature range is `0 <= beta <= 3`, with coupling
magnitude and Boltzmann constant set to 1. There is no external magnetic field.

The first CLI form returns both the entropy for that particular outcome and
the information averaged over all outcomes for the chosen sensor set. Do not
interchange those quantities. The second bond model reverses edge 7, connecting
interior vertex 7 and boundary vertex 6. The `--best-sensors` option reports all
sets within 1e-9 bits of the numerically largest mutual information.

## Python and JavaScript interfaces

```python
from boundary_model import BoundaryModel
m = BoundaryModel()
state = m.evaluate(.64, 'ferro')
all_subsets = m.all_sensor_masks(state)
observation = m.observation(.64, mask=255, values=85, model='ferro')
```

Bit 0 corresponds to the first boundary vertex, and spin bit 0 means −1. The
value `85` therefore means +−+−+−+− in the declared boundary order.

`boundary_engine.js` exports `createBoundaryEngine(data)` in Node.js and can
also be included in a browser script. Supply `browser_coefficients.json` as
`data`; then call `evaluate(beta, model)`, `sensors(state, mask)`, and
`observe(state, mask, values)`. `observe` includes a suggested next sensor and
its expected conditional information, given the current outcome. This is a
one-step recommendation, not a globally optimal adaptive policy.

The in-conversation explorer uses this engine with embedded coefficients. Its
browser calculations passed 1,040 comparisons against the Python model plus
24 checks of conditional sensor advice. Local
browser previews were blocked by the workspace URL policy, so rendered layout,
theme, and interaction QA could not be completed. The three exported scientific
figures were rendered and visually inspected.

## Data map

| File | Meaning |
|---|---|
| `partition.json` | Fixed A/B/C partition, source-order graph edges, selected defect |
| `boundary_coefficients.json` | Exact interior and exterior integer coefficients, whole-graph counts, construction gates |
| `interior_enumeration.npz` | All 16,384 interior bit patterns, internal cuts, and contact patterns |
| `boundary_results.json` | Temperature summaries, numerical peaks, sensor optima, observation examples, controls |
| `sensor_mask_results.npz` | 48 × 256 arrays for each model and information/prediction metric |
| `verification.json` | Interface checks and explicit browser QA limitation |
| `graph143_drawing_positions.npy` | Numerical drawing coordinates only; no metric interpretation |
| `joint_143.json` | Previously verified whole-graph joint count data, retained for provenance |
| `ising_exact.py`, `ising_extensions.py` | Existing companion integer-elimination routines used by this experiment |
| `PROVENANCE.json`, `SHA256SUMS.txt` | Versions, source identity, and delivered file hashes |

The source MTFT license is retained in `MTFT_LICENSE`. The experiment is a
companion bundle; it does not modify or release a new version of MTFT.
