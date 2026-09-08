# MTFT Hecke atlas and stability, v0.1.0

Two bounded experiments using MTFT 0.26.0's X0(143) data: exact prime splitting in the quartic/sextic coefficient fields, and conservative perturbations of the four rational Hecke sectors. Start with `Hecke_Atlas_and_Stability_Report.md` and the two PNG/SVG figures.

The bundle includes frozen input matrices and compact labels from the preceding irregular-prime experiment. No MTFT or PARI installation is needed to reproduce these computations. Python 3.12 with the pinned dependencies was used. The frozen data retain MTFT's source provenance; this experiment is an external companion and does not modify the package.

Run these commands in this directory, in order:

```bash
python -m pip install -r requirements.txt
python build_atlas.py
python build_sectors.py
python simulate_stability.py
python validate_results.py
python render_report.py
```

For a quick inspection:

```bash
python inspect_hecke.py --prime 3511
python inspect_hecke.py --interval 9
python inspect_hecke.py --sectors
```

`prime_atlas.json` contains every prime p < 10201, both modular factorizations, square-interval position, irregularity index, Wieferich status, ramification flags, and a separate level-prime flag. Polynomial coefficients are descending and nonnegative modulo p. At ramified primes, the multiplicity-expanded degree list is recorded for inspection but is not a Frobenius cycle type.

`splitting_types.json` includes all 55 joint types, even those not seen here, exact limiting densities as rational strings, finite counts, and regularity overlays. Infinitude labels refer to the splitting class itself; intersections with irregular primes are not assigned an infinitude theorem.

`field_certificate.json` contains exact discriminants, Galois computations, independent reduction witnesses, and proofs of the group and compositum conclusions. `sector_certificate.json` gives rational CRT projectors, exact Atkin–Lehner intersection dimensions, selection rules, and the matrix iteration lemma.

`stability_arrays.npz` stores transfer with axes `[seed, class, epsilon, time, destination, source]`. Columns are population fractions for a uniform source sector in Hodge-orthonormal coordinates. Other arrays hold perturbations, 13 dimensionless positive frequencies, seeds, and grids. Perturbation statistic columns are off-block Frobenius norm, T2 commutator norm, W11 commutator norm, W13 commutator norm, operator norm, Frobenius norm, and trace. The class order is in `common.py` and `stability_summary.json`.

The plan predates counts/flows. Its arithmetic typo is corrected transparently in `PLAN_ERRATUM.md`; the original plan and hash are retained. Output scripts overwrite only their generated outputs when rerun. `SHA256SUMS.txt` is a delivery manifest, so reruns that change environment-dependent serialization may require rebuilding it.

The algebraic certificates establish finite identities and support the accompanying mathematical arguments. Simulated leakage depends on the specified energy model, perturbations, and time window. No statistical claim about irregular-prime independence, Legendre proof, physical mass prediction, or new mathematical theorem is asserted.
