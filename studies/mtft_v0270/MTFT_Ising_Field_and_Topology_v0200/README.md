# MTFT Ising companion v0.2.0

A reproducible field-and-topology playground built against the supplied **mtft v0.26.0** source archive. This companion adds joint energy/magnetization counts, flat surface twists, bond defects, small bilayers, complex zeros, and exact-weight equilibrium sampling. The bundled MTFT release is unchanged.

Start with `MTFT_Field_and_Topology_Report.md`. Five PNG/SVG figure pairs and `MTFT_Equilibrium_Draws.gif` show the results. The JSON files contain exact integer tables, witnesses, and verification records.

## Try the playground

Use Python 3.10 or later; the recorded run used Python 3.12.13. The precomputed-data playground needs NumPy. Full reconstruction uses the dependencies below and MTFT.

```bash
python playground.py field --beta .64 --eta .03
python playground.py twist --basis 25 --beta .64
python playground.py defect --edge 4
python playground.py bilayer --N 35 --beta .64 --kappa .1
python playground.py sample --a 4 --d 1 --count 5 --seed 143
python playground.py zeros
```

Every command writes JSON to standard output. For example, append `> my_field_result.json` to keep a result. `python playground.py --help` lists the commands.

- `field` evaluates the exact N=143 table at numerical natural parameters β and η. Its Fisher metric is the full covariance matrix of interaction statistic S and magnetization M; it reports Gaussian curvature K.
- `twist` selects one of 26 recorded tree/cotree basis representatives. Indices are zero based. The output includes the changed bond signs and the numerical free-energy cost.
- `defect` reports the complete ground-state result for one AF-to-ferro bond reversal, including a closest witness. Edge 4 requires at least nine flips.
- `bilayer` supports N=6, 11, 15, and 35. κ is the explicitly introduced inter-layer coupling.
- `sample` uses positive integer edge weights `a,d` and site weights `w0,w1`. It has β=log(a/d)/2 and η=log(w1/w0)/2. Bits 0 and 1 mean spins −1 and +1. Weights remain exact integers; displayed β and η are numerical.
- `zeros` summarizes the exact field-zero certificates and numerical temperature-zero calculations. Full coefficients and rational isolating intervals are in `lee_yang.json` and `fisher_zeros.json`.

Numerical evaluation can lose precision when parameters make the ensemble nearly deterministic. The displayed geometry was checked over β∈[0.05,1], η∈[−0.08,0.08]. A null curvature indicates that the numerical determinant guard was reached. The underlying integer counts retain their exactness.

## Reconstruct the experiments

Install into an environment of your choice:

```bash
python -m pip install ./mtft-0.26.0.tar.gz -r requirements.txt
```

`requirements-tested.txt` records the versions used here; it is an environment record, not a promise that those exact versions support every platform. Our run imported the supplied release from its extracted `src` directory and used already available dependencies. Installing the source archive is a standard alternative.

Run from this directory, without Python's `-O` option because the computational gates use assertions:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python ising_exact.py --levels 6 11 15 35 55 77 105 143 --out baseline_results.json
python run_extensions.py
c++ -O3 -std=c++17 twisted_cycle_counts.cpp -o twisted_cycle_counts
./twisted_cycle_counts < twist_cycle_input143.txt > twisted_cycle_counts143.json 2> twisted_cycle_enumeration.log
python analyze_geometry.py
python analyze_zeros.py
python make_figures.py
python verify_artifacts.py
```

The first command regenerates the earlier zero-field data. `run_extensions.py --stages joint twists defects bilayer sampling` can run selected stages; `joint` must precede sampling, and twists generate the C++ input. The C++ verifier enumerates all 2^29 even subgraphs of the N=143 dual graph. It caps cycle rank at 29 and supports at most 128 edges.

`verify_artifacts.py` checks the stored C++ output against the exact transform; it does not silently rerun the enumeration. It also independently checks the rational Lee–Yang intervals by exact sign changes, tests the shipped CLI examples, and writes `verification.json`. The delivered run passed 165 checks: 142 direct exact checks, ten recorded exact construction checks, seven numerical diagnostics, and six CLI integration checks.

## Files and API

| File | Purpose |
|---|---|
| `ising_extensions.py` | `signed_counts`, `IntegerSampler`, `signed_minimum`, `bilayer_counts`, `joint_observables` |
| `ising_exact.py` | Zero-field counting kernel, min-fill order, high-temperature transform, thermal evaluation |
| `playground.py` | Six user-facing commands and reusable `bilayer_observables` |
| `run_extensions.py` | Exact experiments and construction checks |
| `analyze_geometry.py` | Full Fisher metric, three curvature routes, numerical field grid |
| `analyze_zeros.py` | Exact reciprocal reduction, real-root isolation, complex-root diagnostics |
| `twisted_cycle_counts.cpp` | Independent cycle-space enumeration with two homology parities |
| `verify_artifacts.py` | Stored-artifact certificate checks and CLI checks |
| `make_figures.py` | All figures, animation, and flat CSV exports |
| `joint_*.json` | D[k,m]; m is the number of up spins, not magnetization |
| `surface_twists.json` | Flat cocycles, signed counts, gauge controls, two-coordinate projection |
| `defect_census.json` | All 84 reversals, exact ground counts, closest witnesses |
| `bilayer.json` | B[k,j]; k is total intra-layer disagreement, j is layer mismatch |
| `equilibrium_samples.npz` | First 256 configurations and full 20,000-draw histograms per coupling |
| `field_geometry_grid.npz` | β, η, and numerical response/curvature grids |
| `graph143_drawing_positions.npy` | Deterministic drawing positions; no physical metric is assigned |
| `provenance.json`, `SHA256SUMS` | Runtime, source identity, scope, and integrity records |

Python's JSON reader preserves the large integer coefficients. Other consumers should use arbitrary-precision integers; conversion through IEEE-754 numbers can discard exactness. The original source archive and the prior targeted release-test log are included. No compiled verifier binary or dependency directory is distributed.

The counting routines are research instruments for the declared finite models. The table-size guard limits spin-scope entries, not the total byte size of arbitrarily large integer polynomials. Scaling to larger graphs should explicitly budget both elimination width and integer size.
