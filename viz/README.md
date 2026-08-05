# viz/ — React visualizations of the MTFT geometry

Five self-contained React components (hooks only, no external state, no
build config beyond React 18).  Render any of them in a React sandbox,
a claude.ai artifact, or your own app.

| Component | Shows |
|---|---|
| `mtft_enneper.jsx` | Interactive Weierstrass–Enneper minimal-surface explorer (Paper 32 — the modular Enneper surface) |
| `MTFT_HyperbolicTiling.jsx` | Poincaré-disk SL(2,ℤ) tessellation colored by the three X₀(143) newform orbits |
| `MTFT_MonsterFingerprint.jsx` | Deterministic avalanche "fingerprint" built from a 13-round mixer (a JS *lookalike* of MonsterHash, not the real algorithm) |
| `X0_143_BurningMandelbrot.jsx` | Burning-Ship/Mandelbrot hybrid modulated by the 200 Hecke traces (first 50 = the verified `TRACE_TOTALS_50`) |
| `MTFT_DrawnLoop.jsx` | The drawn-loop stage: sketch closed walks on the 56-triangle dual of X₀(143), read off the 29 exact intersection numbers with the Manin arcs, the cuspidal/Eisenstein split, and Zeno survival curves (exact rates 35/16, 10/3); interface data exported by `studies/du02_cycle_space_map.py` |

These live at the repo root rather than inside the package on purpose:
a `mtft/viz/` package directory would shadow the Python module
`mtft/viz.py`, and the wheel has no business shipping JSX.
