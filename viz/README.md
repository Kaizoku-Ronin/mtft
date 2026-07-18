# viz/ — React visualizations of the MTFT geometry

Four self-contained React components (hooks only, no external state, no
build config beyond React 18).  Render any of them in a React sandbox,
a claude.ai artifact, or your own app.

| Component | Shows |
|---|---|
| `mtft_enneper.jsx` | Interactive Weierstrass–Enneper minimal-surface explorer (Paper 32 — the modular Enneper surface) |
| `MTFT_HyperbolicTiling.jsx` | Poincaré-disk SL(2,ℤ) tessellation colored by the three X₀(143) newform orbits |
| `MTFT_MonsterFingerprint.jsx` | Deterministic avalanche "fingerprint" built from a 13-round mixer (a JS *lookalike* of MonsterHash, not the real algorithm) |
| `X0_143_BurningMandelbrot.jsx` | Burning-Ship/Mandelbrot hybrid modulated by the 200 Hecke traces (first 50 = the verified `TRACE_TOTALS_50`) |

These live at the repo root rather than inside the package on purpose:
a `mtft/viz/` package directory would shadow the Python module
`mtft/viz.py`, and the wheel has no business shipping JSX.
