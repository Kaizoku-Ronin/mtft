# MTFT AXG-04 reproduction

Start with `MTFT_AXG04_Report.md`. This is an additive research bundle, not a
new mtft package release. C3X is a chosen six-dimensional EFT candidate with
new family multiplicity, gauge global form, flux and coupling inputs.

## Run

Python 3.12 was used:

```bash
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_all.py
```

The scripts regenerate their JSON results and a combined summary. No network
is needed after dependency installation. No sibling AXG directory or full
mtft installation is required.

| Audit | Content |
| --- | --- |
| `chiral_parent_audit.py` | Original-product controls, independent Cartan reconstruction, C3X subgroup factorization |
| `charged_action_audit.py` | Physical-representation reconstruction, explicit action data, tensor lattice, BF/discrete checks, classical AdS control |
| `higgs_action_audit.py` | Yukawa contractions, exact divisor/evaluation proofs, zero modes and flavor limits |
| `gravity_background_audit.py` | Einstein equations, flux energy, radius solution, radion and shape limitations |
| `larger_parent_audit.py` | U8 Dirac-adjoint alternative, mirrors, enhancement and Yang–Mills Morse index |
| `global_bordism_audit.py` | Ordinary spin-BG bordism calculation using explicit mod-two matrices |

The arithmetic inputs S0²=K, H0(S0)=span(1,u) and the four CM-point values are
inherited from the supplied v0.31.4 arithmetic-spin construction. Their source
is bundled as `input/arithspin_v0314.py`; this investigation proves the new
degree-one divisor consequence and does not claim a new independent proof of
all earlier arithmetic inputs. `input/smflux_v0314.py` provides the original
M1 comparison data. The source license and preceding audit summary are included.

`higgs_action_notes.md` and `larger_parent_notes.md` give further derivations.
`manifest.json` hashes delivered files except itself; Python bytecode and the
local virtual environment are omitted. Different Python patch versions can
change the environment line of the regenerated summary.

The exact checks do not supply a globally defined GS partition function, a UV
completion, observed flavor parameters, a Minkowski/de Sitter vacuum, or a
mechanism selecting the arithmetic complex structure.
