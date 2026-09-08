# MTFT active-module experiment

Read `MTFT_Active_Module_Report.md` for findings and evidence limits.
`LOCALIZATION_RANK_LEMMA.md` proves the support-rank bound. The module is
eight COMPLEX dimensions (16 real), with five complex inactive directions.

The release source is unchanged. The original uploaded archive is included
as `inputs/mtft-0.26.2.tar.gz`, SHA256
`35cfc00c32877bd6c12b16464ddc8190bcca9e2f054474f4ab0a6c605fe9da38`.

## Reproduce

Use Python 3.12 and the package dependencies, plus the versions recorded in
`requirements.txt`. Extract the archive and run:

```bash
tar -xzf inputs/mtft-0.26.2.tar.gz
python run_study.py --source mtft-0.26.2
```

The runner sets the source import path and one OpenBLAS thread. It executes
five scripts in order, records `run_log.txt`, and writes `VALIDATION.json`.
No Sage, GP executable, network lookup, theta sum, or large simulation is
required. Study validation preserves ambiguous package Lie-closure outcomes;
it does not convert them to successful closure certificates.

Scripts and recorded outputs:

- `frame_check.py`: independently verifies operator transport, Hodge adjoints,
  and the distinction between linear and antilinear channels.
- `kernel_check.py`: exact rational support-rank certificate; primitive,
  Gram, and localized-support projector constructions; STAR checks.
- `coupling_study.py`: frozen primary/replication operators, four-block
  coupling norms, coordinate/precision comparisons, and invariant hulls.
- `independent_coupling_check.py`: independent localized projector, direct
  T2 span enlargement, exact commuting pairs, and exploratory compression.
- `render_coupling.py`: standalone PNG/SVG from the recorded data.

`PROTOCOL.md` was fixed before the coupling values were read. Its Amendment A
clarifies channel types. `RUN_NOTES.md` preserves the initial execution stop
and exploratory follow-up. NPZ files retain the projectors and operator
matrices. No object arrays or pickle loading are used.

`SHA256SUMS.txt` covers the delivered snapshot. Reruns may change runtime paths,
timings, float bytes, plot metadata, and NPZ bytes. They do not regenerate the
written report or manifest. Independent routes share the package's frozen
period/arithmetic inputs; independence concerns computation, not source data.
