# MTFT T2 short-horizon study v0.1.0

Read MTFT_T2_Short_Horizon_Report.md for findings and SHORT_HORIZON_IDENTITIES.md
for the proof. The complete numerical geometry is diagnostic; the bound is
exact for the explicitly frozen dyadic Hermitian matrices.

Dependency-free certificate replay:

```sh
python verify_certificate.py
```

To regenerate the numerical diagnostics and exact certificate:

```sh
python -m pip install -r requirements.txt
python run_study.py
python render_results.py
```

Regeneration overwrites generated outputs in this directory. Use a working
copy to preserve the frozen evidence. Input blocks and the selected subspace
are copied unchanged from the T2 Memory Approximation parent. Frozen model
construction is explicitly recorded in the report; FROZEN_MODELS.json allows
exact replay without NumPy, SciPy, the source checkout, or a network connection.

The module dyadic_bounds.py is a reusable research-tool candidate. The supplied
v0.27.0 handoff includes its validation evidence and does not add it to a public
MTFT API automatically.
