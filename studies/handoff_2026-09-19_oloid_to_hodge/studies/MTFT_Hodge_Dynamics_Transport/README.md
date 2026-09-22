# Hodge dynamics and explicit period transport

Read `REPORT.md` first. This continuation investigates the gauge-geometric meaning
of the previously constructed quadratic-map curve tower. Exact algebraic proofs
and numerical diagnostics are kept distinct.

Run from the bundle root:

```bash
python3 root_checks/elliptic_exact.py
python3 geometry/compute_transport.py
python3 dynamics/verify_fricke.py
python3 root_checks/interval_transport.py
python3 monodromy/compute_monodromy.py
python3 monodromy/new_component_monodromy.py
python3 plot_results.py
```

The first three commands need only Python's standard library. Numerical transport
uses NumPy and SciPy; plotting uses Matplotlib. Results are written next to their
programs, and the plot is `transport_comparison.png` at the root.

## Contents

- `root_checks/`: exact elliptic reduction and deliberately imposed chaotic
  interval forcing with bounded flat transport.
- `geometry/`: rational genus-two connection, direct genus-three verification,
  and flat inherited/new projectors.
- `dynamics/`: Fricke isogeny lift, connection compatibility, normalization laws,
  and proofs of the dynamical obstructions.
- `monodromy/`: complex parameter loops, monodromy matrices, gauge checks, and
  bounded/parabolic/hyperbolic word growth.
- `MANIFEST.json`: dependency versions and checksums of the distributed files.

The standalone calculations do not install or modify MTFT. The earlier source
audit used `mtft-0.32.0.tar.gz`, SHA-256
`46ffc563a0cd81a7a82d975497bf64d82f8471eb424409e3d154532b5d8d6069`.
No new source release, full-package test gate, interval certification, or
physical-constant prediction is claimed.
