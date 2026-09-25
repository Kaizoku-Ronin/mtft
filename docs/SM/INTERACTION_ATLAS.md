# Standard Model interaction atlas

The atlas organizes **reference interaction vertices** and the evidence needed
to connect MTFT calculations to them. It supplies diagrams and symbolic tensor
data, not a new SM Lagrangian, amplitude calculator or matter simulator.

The finite list of vertices in a specified model generates indefinitely many
process diagrams at increasing external multiplicity and loop order. An atlas
of vertices therefore has a different scope from a repository of every diagram.

## Open the atlas

From a checkout with `pip install -e .`:

```bash
python -m mtft.interactions summary
python -m mtft.interactions render -o atlas-output/index.html
```

Open `atlas-output/index.html` in a browser. The generated file is self-contained:
no server, internet connection, Mathematica, MadGraph installation or JavaScript
package manager is required. It includes search, sector filters, diagrams,
source links, symbolic tensors, parameter definitions, research pointers and
JSON/SVG downloads. The source links need internet access when clicked.

```bash
python -m mtft.interactions list --sector qcd
python -m mtft.interactions list --particle e- --physical-only
python -m mtft.interactions show V_6
python -m mtft.interactions svg V_6 -o atlas-output/four-higgs.svg
python -m mtft.interactions export -o atlas-output/sm-catalog.json
```

Existing output files are protected; add `--force` to replace an export.
Filters accept an exact UFO particle symbol, name or **signed** PDG code.
`--physical-only` hides vertices with ghosts or Goldstones. It does not perform
a change of gauge or authorize omission of those fields in loop calculations.

## Reference coverage

The snapshot is normalized from MadGraph5_aMC@NLO's `models/sm` at commit
[`abfd2c92873b3cfe4baa65c661b580e5cc07d45f`](https://github.com/mg5amcnlo/mg5amcnlo/tree/abfd2c92873b3cfe4baa65c661b580e5cc07d45f/models/sm).
It contains **153 vertices, 43 particle records and 108 symbolic couplings**.
There are 71 vertices without ghosts or Goldstones.

| Sector | Vertices |
|---|---:|
| QCD | 8 |
| Electroweak fermion currents | 45 |
| Electroweak gauge self-interactions | 6 |
| Higgs Yukawas | 6 |
| Higgs–gauge interactions | 4 |
| Higgs self-interactions | 2 |
| Goldstone interactions | 58 |
| Ghost interactions | 24 |

This is the complete `vertices.py` inventory **of that stock model**. Its scope
is explicitly limited:

- The upstream model fixes the up, down and strange quark masses to zero; their
  Higgs Yukawa vertices are absent. The viewer lists these as three coverage
  gaps. Their associated Goldstone Yukawa terms are also affected by the same
  approximation. A fully massive SM export is a future reference expansion.
- Neutrinos are massless. CKM entries use the upstream truncated Wolfenstein
  parametrization; this is not an exact arbitrary unitary CKM matrix.
- No `restrict_default.dat` is applied. That card would remove additional
  masses, Yukawas and mixing. Intrinsic assumptions in the Python model remain
  even with no restriction card.
- Gauge fixing is Feynman gauge, with the ghost and Goldstone vertices present.
- This is tree-level data. It does not contain NLO counterterms, loop amplitudes
  or effective loop-induced Higgs–gluon/Higgs–photon interactions.
- The `QED` coupling-order label includes electroweak and Higgs powers; it does
  not label every such vertex as electromagnetic.
- Numerical parameter defaults are retained as source data. They are not a
  current measurement compilation or predictions of MTFT.

Particle order is preserved. Tensor expressions use one-based leg indices;
the coupling dictionary uses zero-based `(color_index, lorentz_index)` slots.
The sum over slots is the vertex structure. Individual diagrams and slots are
not generally separately gauge-invariant observables.

The normalized JSON includes parameters, function expressions, Lorentz
structures, coupling orders, source line numbers and SHA-256 hashes for all
seven parsed files. Expressions remain strings and are never evaluated.
The upstream license is shipped in
[`MADGRAPH_LICENSE.txt`](../../src/mtft/interactions/_data/MADGRAPH_LICENSE.txt)
and embedded in JSON/HTML exports so standalone copies retain attribution.

## Record MTFT evidence

All initial rows are `unmapped`. Related-code links are research suggestions,
not evidence of a match. In particular, the surface Yang–Mills toolkit is not
automatically a derivation of the four-dimensional QCD vertices.

```bash
python -m mtft.interactions init-ledger -o atlas-output/mtft-sm-ledger.json
# Edit this JSON as the research produces evidence.
python -m mtft.interactions validate-ledger atlas-output/mtft-sm-ledger.json
python -m mtft.interactions render --ledger atlas-output/mtft-sm-ledger.json -o atlas-output/mapped.html
```

| Status | Meaning | Required record |
|---|---|---|
| `unmapped` | No supported correspondence recorded | A row for the vertex |
| `symmetry-compatible` | Candidate fields/representations can support the interaction | Evidence references with assumptions |
| `vertex-derived` | An explicit normalized interaction has been derived | Derivation references, including color, Lorentz/chiral and normalization conventions |
| `amplitude-tested` | A process-level comparison has passed the recorded numerical gate | Derivation evidence and at least one complete benchmark |

A ledger is bound to the **content hash of the entire catalog**. Changing a
reference, restriction, coupling or convention requires an explicit migration;
the same vertex number in two exports need not mean the same interaction.
Every reference vertex must occur exactly once.

An amplitude benchmark records:

```json
{
  "process": "declared external states",
  "energy_gev": 200.0,
  "perturbative_order": "tree-level",
  "observable": "spin/color averaged squared matrix element at specified kinematics",
  "unit": "declared units",
  "conventions": "gauge, masses, coupling scheme, phase-space point, basis and normalization",
  "reference_artifact": "path/to/reference-result.json",
  "mtft_artifact": "path/to/mtft-result.json",
  "mtft_revision": "commit used for the calculation",
  "reference_value": 1.0,
  "mtft_value": 1.0,
  "absolute_tolerance": 1e-10,
  "relative_tolerance": 1e-8
}
```

The numbers above illustrate the schema; they are **not computed results**.
For `amplitude-tested`, the validator checks finite numbers and
`abs(mtft - reference) <= absolute_tolerance + relative_tolerance * abs(reference)`.
Record tolerances before comparisons. Complex amplitudes should be compared
through separately specified real observables or components.

Validation checks the ledger's structure and recorded comparisons. It does not
read the referenced papers, verify the derivation or rerun either amplitude.
Physical review must still examine those artifacts and all contributing diagrams.

## Import and reproduce

The importer accepts a conservative **declarative tree-level SM UFO subset**:
literal constructor assignments, particle `.anti()` declarations and imports
(which are not executed). Unsupported dynamic statements, unresolved tensor
references and counterterm files fail explicitly. It is not a general UFO
interpreter. A trusted general UFO implementation should be consumed by a
dedicated tool such as MadGraph, outside this static reader.

```bash
python -m mtft.interactions import-ufo /path/to/model \
  --source-url https://github.com/mg5amcnlo/mg5amcnlo \
  --revision abfd2c92873b3cfe4baa65c661b580e5cc07d45f \
  --model-path models/sm \
  --gauge 'Feynman gauge' \
  --restrictions 'No restriction card; stock intrinsic mass assumptions' \
  -o atlas-output/imported.json
python -m mtft.interactions --catalog atlas-output/imported.json summary
```

To reproduce the bundled curated snapshot byte-for-byte (requires network):

```bash
python scripts/interactions/update_reference.py --check
# Or supply already downloaded model files:
python scripts/interactions/update_reference.py --local /path/to/models/sm --check
```

The bundled scope notes describe the pinned model. For any new revision, review
the masses, mixing conventions, gauge, license and exclusions before updating
the snapshot. Rebuild dependent ledgers deliberately. Never silently replace
the pin with a floating branch.

## MadGraph process diagrams

MadGraph is optional and installed separately. The supplied
[`reference_processes.mg5`](../../scripts/interactions/reference_processes.mg5)
generates reference process code and diagrams; it does **not** run integration
or claim an MTFT prediction. In a MadGraph installation:

```bash
./bin/mg5_aMC /absolute/path/to/mtft/scripts/interactions/reference_processes.mg5
```

It selects `sm-full` to avoid the default restriction card and explicitly sets
Feynman gauge. The stock model's intrinsic light-quark assumptions still apply.
The output is a MadGraph reference project in the working directory. Parameter
cards, kinematics and compiler requirements belong to that external workflow.
Deriving an MTFT UFO and comparing amplitudes is subsequent physics work.

References: [FeynRules SM implementation](https://cp3.irmp.ucl.ac.be/projects/feynrules/wiki/StandardModel),
[UFO 2.0 specification](https://arxiv.org/abs/2304.09883),
[MadGraph5_aMC@NLO](https://github.com/mg5amcnlo/mg5amcnlo),
[FeynArts](https://feynarts.de/).
