# Reproduction and verification

## Verify the delivered bundle first

From the outer handoff directory:

```bash
python3 checks/verify_bundle.py
```

This uses only Python's standard library. It verifies the combined manifest, the original study manifests, original ZIP integrity, and the expected MTFT source hash. It does not rerun mathematical experiments or certify their conclusions. Verify before rerunning programs, because reproduction writes new output files and may change their hashes.

## Environment

The recorded Hodge transport environment was Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, and Matplotlib 3.10.8. `requirements-reproduction.txt` pins these three Python packages. The standalone study programs do not require installation of MTFT. Work in a copy if preserving the delivered results is important.

An optional isolated environment can be prepared with:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements-reproduction.txt
```

On Windows, activate the virtual environment using the command appropriate for your shell. The dependency installation may require network access. The programs themselves use the included data.

## Study 1: graded geometry

Run from `studies/MTFT_Graded_Moonshine_Investigation/`:

```bash
python3 verify_graded.py
python3 audit_operators.py
python3 audit_bracket.py
```

The operator audit is standard-library-only; the other two need NumPy. They regenerate `graded_results.json`, `operator_results.json`, and `bracket_results.json`. The exact lower-bound/upper-bound argument and frozen-input assumptions remain necessary when interpreting rank results.

## Study 2: triangle, Feigenbaum, and Hodge tower

Run from `studies/MTFT_Triangle_Feigenbaum_Hodge_Investigation/`:

```bash
python3 scaling/compute_scaling.py
python3 geometry/verify_geometry.py
python3 plot_results.py
```

The first two need only the standard library. Plot regeneration needs Matplotlib. The stored calculations reach the principal superstable period 1024 with 60/90-digit comparisons. Numerical root searches and precision agreement are not interval certification.

## Study 3: Hodge transport

Run from `studies/MTFT_Hodge_Dynamics_Transport/` in the order shown:

```bash
python3 root_checks/elliptic_exact.py
python3 geometry/compute_transport.py
python3 dynamics/verify_fricke.py
python3 root_checks/interval_transport.py
python3 monodromy/compute_monodromy.py
python3 monodromy/new_component_monodromy.py
python3 plot_results.py
```

The first three use only the standard library. ODE transport needs NumPy/SciPy; plotting needs Matplotlib. The order allows later programs to reuse the generated exact connection data. Numerical solver variation may alter last digits and output checksums without changing an exact theorem.

## Oloid follow-up

`notes/OLOID_RETURN.md` gives the source equation, the rational conjugacy derivation, the physical parameter range, and the limitations. No rolling simulation, new ODE experiment, or configuration-space/Hodge lift accompanies that note. Its rational identities can be checked directly by substitution.

## Historical source snapshot

`sources/mtft-0.32.0.tar.gz` is the exact distribution used in the earlier source audit. The original study reports did not carry this tarball inside their individual bundles; it is now included in the outer handoff for convenience. It is retained as an archive to avoid confusing a historical snapshot with an active checkout. The full MTFT suite was not rerun in these studies or during consolidation.

## Validation actually performed while packaging

- All 47 payload hashes in the three existing study manifests were checked.
- Original ZIP CRC checks passed; standalone reports matched the reports inside the ZIPs.
- The included source archive matched the SHA-256 recorded in all three studies.
- The final handoff was checked against its manifest and the resulting ZIP was tested for corruption.

The detailed results are in `provenance/PACKAGING_REPORT.json`. No additional mathematical test pass is implied by these packaging checks.
