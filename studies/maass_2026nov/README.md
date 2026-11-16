# studies/maass_2026nov — the Maass spectrum of Γ₀(143): provenance of `mtft.maass`

The package's subpackage `mtft.maass` (v1.1.0) was ported from the Modular Plane lab's records MP14 to MP20, 12 to 14 November 2026.
This folder keeps what that port rests on.

- `build_frozen_spectrum.py` builds the package's frozen table (`src/mtft/maass/_data/`): for each of the lab's 83 forms it starts the
  package's own solver within 10⁻⁵ of the lab's eigenvalue and keeps what the solver finds. About six minutes. Run from the repository
  root. Its last run: 83 forms; largest change from the lab's value 5.0 × 10⁻¹²; two horocycles agreeing to 1.2 × 10⁻¹¹; Hecke relations
  to 2.3 × 10⁻¹¹.
- `lab/` holds the lab's scripts and number files as they were, byte for byte. They are records, not package code: they import each
  other by relative path, write their outputs beside themselves, and some need scipy or cypari2.

| Folder | Record | What it did |
|---|---|---|
| `lab/mp14` | MP14 (T-G10) | the scattering matrix of the four cusps; the Eisenstein series; Φ(½) = −I |
| `lab/mp15` | MP15 | Hejhal's solver on Γ₀*(143); the first lists below R = 1.6; controls at levels 133 and 145 |
| `lab/mp16` | MP16 | the trace formula twisted by the Atkin–Lehner involutions and the mirror, for odd forms (needs cypari2) |
| `lab/mp17` | MP17 | the untwisted odd trace: the odd lists checked sector by sector |
| `lab/mp18` | MP18 | the η-twist: 42 forms below R = 1.8; the pole-free detector; the twisted trace formula |
| `lab/mp19` | MP19 | the even sectors with the continuum; the untwisted spectrum complete below R = 1.98 |
| `lab/mp20` | MP20 | the Ford domain's quadrature; participation ratios; a Monte Carlo second route |

- `build_trace_weights.py` recomputes the geometric weights of the trace formulas for lengths up to 17 with cypari2 and writes
  `src/mtft/maass/_data/trace_weights.npz`. Given the folder of the lab's stored tables it compares: they agree exactly. About ten seconds.

What is in the package and what is not. The solver, the reduction, the frozen table, the scattering matrix, the Ford domain and the trace
formulas are in `mtft.maass`, with their gates in `tests/test_maass_*.py`. The trace formulas' geometric side needs class numbers and
regulators of real quadratic orders, which PARI supplies (conditional on GRH in general; the package checks a sample against Dirichlet's
finite sum). The package freezes the weights to length 17, which a test function of width 0.43 to 0.5 needs, and certifies the lists
complete below R = 1.959 (untwisted) and 1.607 (twisted, together). The lab's sharper radii, 1.98 and 1.8, came from a width of 0.35
and lengths to 23.5: 126,751 traces per coset, megabytes of weights that are not stored here. They are a recorded result of these scripts.
