# MTFT Hecke-block coupling study

Read `MTFT_Hecke_Block_Report.md` and `W143_COUPLING_LEMMA.md`.
The figure is `Hecke_Block_Map.png` (also SVG). The report keeps the exact
arithmetic decomposition separate from the numerical local active/fixed split.

## Reproduce

Use Python 3.12 with `requirements.txt`. The supplied release archive is
included under `inputs/mtft-0.26.2.tar.gz`. From this folder:

```bash
tar -xzf inputs/mtft-0.26.2.tar.gz
python run_study.py --source mtft-0.26.2
```

The default prior-input directory is `inputs/previous`. It contains the
two frozen numerical matrix files from MTFT_Active_Module_v0.1.0 required
for the primary and independent routes. They are input data, not results
re-estimated in this stage. Their provenance and original manifest are
included there. To use a freshly replayed prior study instead, pass
`--prior /path/to/MTFT_Active_Module_v0.1.0`.

The runner uses one OpenBLAS thread and invokes:

1. `exact_blocks.py`: rational CRT/basis projector comparison, ten exact
   operator-preservation checks, exact sign data and W143 rank explanation.
2. `block_coupling_study.py`: direct-kernel local projector versus exact CRT
   arithmetic projectors; three coordinate seeds; angles, ranks, commutator
   attribution and signed Gram tables.
3. `independent_block_geometry.py`: localized-support projector versus exact
   arithmetic bases with independent QR; comparisons with all primary runs.
4. `render_block_map.py`: scientific plot from recorded results.

`run_log.txt` preserves output. `VALIDATION.json` contains the study checks.
No Sage, GP executable, external lookup, theta sum, or Lie-closure run is
required. The original source archive is included unchanged; no patch is
applied. Its SHA256 is
`35cfc00c32877bd6c12b16464ddc8190bcca9e2f054474f4ab0a6c605fe9da38`.

`PROTOCOL.md` was frozen before the new observables were read.
`RUN_NOTES.md` identifies explanatory follow-ups. The exact projector JSON
uses rational strings; NPZ files use ordinary numeric arrays with no pickle.
Matrix magnitudes and saturated ranks remain double-precision diagnostics.

The delivered `SHA256SUMS.txt` is a snapshot manifest. Reruns may change
timings, paths, binary plot/NPZ bytes and floating-point bytes; they do not
regenerate the written report or manifest. The independent routes share
upstream period and arithmetic data.
