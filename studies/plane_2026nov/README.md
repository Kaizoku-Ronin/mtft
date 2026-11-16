# studies/plane_2026nov — `mtft.plane` against the Modular Plane lab's own core

`mtft.plane` (added 16 November 2026) is a port: its kernel and its quantity lattice come from the core of the Modular Plane lab,
version 8.3, which is written in JavaScript. This folder holds the comparison of the port with the original.

| file | what it is |
|---|---|
| `make_inputs.py` | writes `inputs.json`: 52 matrices (12 named, 40 random, seeded) and 49 sites of the quantity lattice (9 named, 40 random, seeded) |
| `probe_core.cjs` | runs the lab's core on those inputs with Node and prints what it computes |
| `lab_core_outputs.json` | that output, frozen on 16 November 2026 (Node 22.23.4), with the core's SHA-256 |

The core used: `f1ac4b53368bff4d0ac32da43d3f735c4438a46c6e9b3e72dcd69e30c927b809` (SHA-256). It is the lab's file, with the lab's
instruments in it, and it is not shipped here; the lab keeps it. With a copy:

    python make_inputs.py
    node probe_core.cjs /path/to/the/core.js inputs.json > lab_core_outputs.json

## What was compared, and the result

- **The kernel** (`tests/test_plane_kernel.py::test_the_port_equals_the_lab_core`). For each matrix: its kind, its region in the plane
  of trace and determinant, its decomposition into dilation, two strains and a rotation, the shape and the numbers of its numerical
  range (centre, both semi-axes, numerical radius, spectral norm), the matrix exponential at t = 0.7, and where they are defined its
  length, Rademacher's Ψ, its order and its Jacobi monodromy. The 52 matrices meet all three kinds, all four shapes of the range
  (ellipse, disk, segment, point) and nine regions. The package agrees with the lab's core to 3.0 · 10⁻¹⁶ relative on every number,
  and exactly on every label and integer. Both of the package's routes to the range were compared: its closed forms and its copy of
  the lab's own method (`kernel.lab_numerical_range`).
- **The prospector** (`tests/test_plane_quantities.py::test_the_prospector_equals_the_lab_core`). For each site: the catalogued
  quantities that live there, the pairs that build it, the number of mates, its mass dimension and charge parity, its dimension
  string and its Planck coefficients. All 49 sites agree name for name (7 of them occupied); the index of the Planck sublattice is 4
  on both sides.

A wider run of the same comparison, on 312 matrices and 209 sites, was made once on the same day and gave the same result
(worst relative difference 1.2 · 10⁻¹⁵, no site differing); it is recorded in the release's build record and its output is not kept.

These tests compare the port with the original. That the original is right is tested elsewhere in `tests/test_plane_kernel.py`, by a
route that uses no formula for the numerical range (numpy, the top eigenvalue of the Hermitian part of e^{−iφ}M).
