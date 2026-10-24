# Sphereland check scripts

Each script certifies the numbers shown by one or more instruments of
[Sphereland](../index.html). It prints PASS or FAIL for every check and ends
with `ALL CHECKS PASSED`, or exits non-zero. The scripts are research
certificates kept beside the page so that every claim on it can be re-run.
They are not part of the `mtft` package and are not shipped in its sdist or
wheel. Their docstrings keep the notes on where each question came from.

| Script | Instruments | What it certifies |
|---|---|---|
| `siegel_slices_20261018.py` | 8, 9 | bridge presentations of the trefoil, the genus of a slice's branched cover, the slice's modular parameter by the AGM, Clebsch's invariants |
| `modular_knots_20261018.py` | 10 | the trefoil complement as SL(2,ℝ)/SL(2,ℤ), modular knots and Rademacher's linking function |
| `hopf_qubit_20261018.py` | 11 | the Hopf fibration as a qubit, the Berry phase, the quaternionic Hopf map and two-qubit entanglement |
| `ym_gap_20261019.py` | 12 | Derrick scaling, the two-loop Λ, strong-coupling SU(2) observables by Weyl integration |
| `kakeya_20261019.py` | 13 | Besicovitch's Perron tree and finite-field Kakeya sets |
| `period_clock_20261019.py` | 14 | repeating expansions, repunits, cyclotomic polynomials, the trefoil's Alexander polynomial, cyclic covers |
| `bessel_knots_20261019.py` | 15 | knotted nodal lines of Bessel beams and the zeta function they carry |
| `uniform_gap_20261020.py` | 16 | spectral gaps under refinement: Cayley graphs, Ramanujan graphs, lattice gauge theory |
| `casimir_knot_20261020.py` | 17 | vacuum energies as constant terms of modular forms, Dedekind sums, linking numbers |
| `gravity_knot_20261021.py` | 18 | Kepler's cusp, its link, and the gravitational-bridge formulas |
| `kronecker_thermal_20261022.py` | 19 | Kronecker's limit formula, the Eisenstein Laurent expansion, the balanced source |
| `poisson_mass_20261022.py` | 20 | Rayleigh–Lamb modes, Michell's thresholds, Fierz–Pauli, Willmore integrals, dislocation energies |
| `aether_torus_20261024.py` | 21 | four aethers, the torus as a modular medium, the inverse-square law and kissing numbers, Kepler's coil |
| `modular_flip_20261024.py` | 22, 23, 24 | the modular flip and the Farey tiling, AdS₅ glueballs and quasinormal modes, Burgers' strips, Maass forms by Hejhal's method, the scattering phase |

Instrument 16 also quotes a Monte Carlo of the repaired SU(2) action. That
script is a tool, not a check, and it already lives with CC-69 as
[`scripts/cc69/bulk_mc_su2.py`](../../../scripts/cc69/bulk_mc_su2.py):
`python scripts/cc69/bulk_mc_su2.py L kappa out.json [sweeps] [betas...]`.

## Running

```bash
python viz/sphereland/checks/modular_flip_20261024.py
```

Requirements: `numpy`, `scipy`, `sympy` and `mpmath`, and `shapely` for
`kakeya_20261019.py` (tested with Python 3.13).
`python-flint` and `gmpy2` are optional and only make the high-precision parts
faster. Run times range from seconds to several minutes; `modular_flip_20261024.py`
takes about six, mostly for the Burgers poles at up to 270 digits.

`python modular_flip_20261024.py --json data.json` also writes the data that
Instruments 22–24 embed (the soliton glueballs, the quasinormal modes, the
Burgers poles and the scattering phase).
