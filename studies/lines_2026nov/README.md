# studies/lines_2026nov — the closed lines of Γ₀(143)\PSL(2, ℝ): provenance of `mtft.lines`

The package's subpackage `mtft.lines` (v1.1.0) was ported from the Modular Plane lab's records MP12, MP13 and MP14, all of
12 November 2026. This folder keeps what that port rests on.

- `build_frozen_lines.py` builds the package's frozen record, `src/mtft/lines/_data/lines_143.json`, from the lab's three number files,
  and writes `PROVENANCE.txt` beside it. It copies values and computes nothing. Run from the repository root; it takes a second.
- `lab/` holds the lab's scripts and number files as they were, byte for byte. They are records, not package code: they read each
  other's number files by path, and Route A of MP12 needs sympy.

| Folder | Record | What it did |
|---|---|---|
| `lab/mp12` | MP12 (test T-G4) | what replaces the Rademacher number at level 143: the thirty integers of 28 closed lines. `mp12_tg4_verify.py` is Route A; `mp12_tg4_mtft.py` is Route B, which ran on `mtft.hecke` of mtft 1.0.0 |
| `lab/mp13` | MP13 (test T-G9) | which unions of closed lines bound, and their linking numbers: every crossing of 24 geodesics by exact cutting (`mp13_tg9_verify.py`), and the Lorenz template (`mp13_lorenz.py`) |
| `../maass_2026nov/lab/mp14` | MP14 (test T-G10) | why the Pell lines meet nothing algebraically: three more motions, 38 geodesics, the mirror rule. Kept once, in the other study folder, because its second half is about the Eisenstein series |

The pages' own JavaScript checks of these records (`test_core_tg4.js`, `test_core_tg9.js`) ran on the lab page's core, which is teaching
material and is not in this repository.

What is in the package and what is not. The graph, the classes, the windings, the lines, the periods by both routes, the involutions,
the cutting, the crossings and the Lorenz template are in `mtft.lines`, each computed by the package's own code and gated against the
frozen record in `tests/test_lines_*.py`. One part of the record is used and not only compared: the lab's integer coordinates of the 29
generators on the four Hecke blocks, a choice of bases that keeps its table of periods reproducible; `mtft.lines.periods.compare_bases`
checks that choice against the package's own Manin symbols. Not reproduced: MP14's shadows on X₀(13), a level with elliptic points,
which `mtft.lines` refuses; and the lab's numerical check that the three η-quotients are functions on the group (160 digits), for which
the package has two exact routes instead.
