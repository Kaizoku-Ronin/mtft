# Exact Ising experiment for MTFT v0.26.0

Start with `MTFT_Ising_Exact_Report.md`. The two PNG figures preview the results;
SVG copies are provided for editing and publication. The original MTFT release
is unchanged. These scripts are a separate experimental companion.

## Reproduce

Python 3.9+ with MTFT, NumPy, SciPy and Matplotlib is required for the full
workflow. The supplied archive pins the MTFT source used here:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install source/mtft-0.26.0.tar.gz scipy matplotlib
python ising_exact.py --levels 6 11 15 35 55 77 105 143 --out results.json
```

This regenerates the exact density of states, two elimination orders,
small-level direct spin histograms, cut/cycle transform, and floating checks.
On a slow machine the 24-spin direct enumeration takes longer than the
polynomial computation. The exact counter itself needs only NumPy, with MTFT
needed to construct a Manin graph.

For the fully independent N=143 cycle enumeration, use a GCC- or Clang-compatible
C++17 compiler:

```bash
g++ -O3 -std=c++17 cycle_counts.cpp -o cycle_counts
./cycle_counts < graph143.txt > cycle_counts143.json 2> cycle_enumeration.log
python analyze_results.py
```

The checked edge list `graph143.txt` is included. `analyze_results.py` asserts
coefficient-by-coefficient equality of the C++ histogram and the cut/cycle
transform before generating the figures. C++ is needed only for this independent
verification; evaluating the thermal curves from the supplied counts needs no
compiler or PARI/GP.

## Use the counts at a new temperature

```python
import json
from ising_exact import thermodynamics

data = json.load(open("results.json"))
row = next(r for r in data["levels"] if r["N"] == 143)
observables = thermodynamics(
    row["density_of_states"], row["spins"], row["edges_count"], beta=0.65
)
print(observables)
```

## Use the exact counter on a different graph

```python
from ising_exact import packed_dos

# Four-vertex cycle; multiple edges and self-loops are supported too.
counts, certificate = packed_dos(4, [(0,1), (1,2), (2,3), (3,0)])
print(counts)  # [2, 0, 12, 0, 2]
```

Every vertex must be numbered from zero through n−1. A width guard prevents
inadvertently creating huge intermediate tables. The measured order width is
an upper bound, not an optimal-treewidth claim.

Counts are exact integers. Thermal evaluations and located peaks are numerical
finite-graph diagnostics. The total partition function does not distinguish
the individual spin-structure contributions of the Pfaffian formulation.

## Files

- `results.json`: eight levels; graphs, full integer counts, checks, orders and timings.
- `analysis.json`: thermal peak estimates, N=6/11 graph correspondence, max-cut witness.
- `density_of_states_143.csv`: exact coefficients and energy levels. The last column
  counts even subgraphs by **occupied-edge count** equal to the row index; it is
  a different count from spin configurations with that disagreement count.
- `thermodynamics.csv`: numerical curves evaluated from exact counts.
- `cycle_counts143.json`, `graph143.txt`: independent cycle counts and their input graph.
- `release_checks.log`, `experiment.log`, `cycle_enumeration.log`: run records.
- `provenance.json`, `SHA256SUMS.txt`: environment, source and artifact provenance.
- `source/mtft-0.26.0.tar.gz`: the supplied release used for the experiment.

The custom Python and C++ companion code is supplied under the MIT license in
`LICENSE`. The MTFT archive retains its own included license and attribution.
