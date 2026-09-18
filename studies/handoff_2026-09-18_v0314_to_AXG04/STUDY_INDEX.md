# Study index and chronology

The sequence below follows the investigations in this conversation, starting with the supplied v0.31.4 candidate. Original reports and scripts retain their contents. Folder names have been organized for this handoff; `provenance.json` maps each extracted file to its original ZIP member.

| Order / date | Study and report | Recorded verification | Main result / physical interpretation |
|---|---|---|---|
| 0 — 16 Sep | [v0.31.4 spin and vacuum audit](studies/00_v0314_audit/MTFT_v0314_Spin_and_Vacuum_Audit.md) | Symbolic audit PASS; seven targeted shipped tests | Theta-compatible Q8 correction; M2 reconstruction diagnostic; SUSY-flatness claim unestablished; restricted 6D radius runaway |
| 1 — 16 Sep | [SC7-01 spin circle](studies/01_SC7_01/MTFT_SC7_01_Spin_Circle_Investigation.md) | Implemented exact checks PASS; 36 primary local-system configurations and inverse-character duality checks | Spin-circle topology is valid; required smooth 7D bulk spectrum FAILS; restricted radius gate FAILS |
| 2 — 16 Sep | [Arithmetic geometry and spacetime study edition](teaching/MTFT_Arithmetic_Geometry_and_Spacetime_Study_Edition.docx) | 187 exact hand checks; original 69-page document QA | Teaching synthesis, explicit matrix models, glossary, exercises and solutions through SC7-01 |
| 3 — 17 Sep | [HOPF-02](studies/02_HOPF_02/MTFT_HOPF02_Report.md) | 85 exact checks, zero failures | Arithmetic Hopf reduction; Hecke quaternionic obstruction/control; CRT; Hodge sphere join |
| 4 — 17 Sep | [AXG-01](studies/03_AXG_01/MTFT_AXG01_Report.md) | 159 exact checks, zero failures | Full local 4D two-axion anomaly factorization, vector masses and operator dressings; no classical radius stabilization |
| 5 — 17 Sep | [AXG-02](studies/04_AXG_02/MTFT_AXG02_Report.md) | 84 exact checks, zero failures | 6D anomaly obstruction; zero-index chirality sensitivity; native-axion flux constraint; conditional two-field AdS control |
| 6 — 17 Sep | [AXG-03](studies/05_AXG_03/MTFT_AXG03_Report.md) | 149 exact checks, zero failures | Tensor signature obstruction; integral factor of three; candidate Z3 fermion anomaly; elementary-Higgs and 6D Yukawa obstructions |
| 7 — 18 Sep | [AXG-04 / C3X](studies/06_AXG_04/MTFT_AXG04_Report.md) | 322 exact checks, zero failures | Explicit charged parent/action, pure unit-flux modes, integral local GS factorization, bordism calculation, AdS control; global completion and realistic vacuum open |

These counts are **historical records in the included ledgers**, not a claim that all studies were rerun during packaging. Do not sum them into a claim of independent theorems, physical evidence, or package-release coverage. The audit and SC7 scripts do not use the same counted-assertion ledger as the later studies; no count has been invented for them.

## Source and result locations

| Study | Entry point(s) | Primary ledger(s) |
|---|---|---|
| v0.31.4 audit | `studies/00_v0314_audit/audit_v0314.py` | `results.json`, `TEST_RESULTS.md`, `run.log` |
| SC7-01 | `studies/01_SC7_01/investigate_spin_circle.py` | `results.json`, `run.log`, `PREREGISTRATION.md` |
| Teaching | `teaching/paper_work/verify_hand_calculations.py` | `hand_calculation_checks.json`, `final_qa.json` |
| HOPF-02 | `studies/02_HOPF_02/investigate.py` | `results.json`, `PREREGISTRATION.md` |
| AXG-01 | `studies/03_AXG_01/run_all.py` | `results_summary.json`, anomaly/operator/vacuum results |
| AXG-02 | `studies/04_AXG_02/run_all.py` | `results_summary.json`, parent-anomaly and flux/scalar results |
| AXG-03 | `studies/05_AXG_03/run_all.py` | `results_summary.json`, four audit ledgers |
| AXG-04 | `studies/06_AXG_04/run_all.py` | `results_summary.json`, six audit ledgers, two derivation notes |

## Original archives

- [v0.31.4 audit bundle](original_bundles/MTFT_v0314_Audit_Bundle.zip)
- [SC7-01 bundle](original_bundles/MTFT_SC7_01_Reproducible_Bundle.zip)
- [Teaching reproduction bundle](original_bundles/MTFT_Study_Edition_Reproduction.zip)
- [HOPF-02 bundle](original_bundles/MTFT_HOPF02_Reproducible_Bundle.zip)
- [AXG-01 bundle](original_bundles/MTFT_AXG01_Reproducible_Bundle.zip)
- [AXG-02 bundle](original_bundles/MTFT_AXG02_Reproducible_Bundle.zip)
- [AXG-03 bundle](original_bundles/MTFT_AXG03_Reproducible_Bundle.zip)
- [AXG-04 bundle](original_bundles/MTFT_AXG04_Reproducible_Bundle.zip)

## Teaching and visual coverage

The study edition contains number-space maps; rotations, boosts and shears; V4/D8/Q8 matrices; Lorentz and Clifford matrices; modular generators; cusp permutations and character transforms; Hodge/Hecke blocks; triality controls; the Fano plane and octonion multiplication; the Gray graph; index and anomaly calculations; Yukawa normalization; and the spin-circle/radius tests. See its matrix/theorem locator, glossary and worked solutions.

The original 13 figure files and their deterministic generator are in `teaching/paper_work/figures/` and `make_figures.py`. Large exact matrices and permutations that would be cumbersome on paper are retained in the study JSON files, particularly HOPF-02 and AXG-04's bordism ledger.

## Inputs and scope boundary

`baseline/` contains the v0.31.4 candidate actually supplied, not a fetched PyPI or GitHub substitute. The two earlier Atlas/Arithmetic Geometry documents are included under `inputs/earlier_study_documents/` as context. Earlier TRI-01 results are discussed and reproduced in part by the teaching edition, but the separate pre-v0.31.4 investigation is outside this handoff's start boundary.

The supplied September 16 conversation is in `inputs/conversation/`. The wiki PDFs and three supplied images are reference inputs, not evidence that their suggested physical identifications were established. HOPF-02 is the label of the delivered Hopf investigation; no separate HOPF-01 deliverable is asserted here. Unrelated figurate-zeta and other research branches have not been mixed into this window.
