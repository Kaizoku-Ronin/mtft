# mtft v0.31.3 — KK-TOWER-02: tower certification, FEM Yukawa cross-check of CC-26, KK decay overlaps (2026-09-14)

`surface.magnetic`: `MagneticMesh.eigenmodes` / `divisor_modes` (eigenpairs in the unitary frame, kinetic-orthonormal),
`m1_fem_yukawa` (Yukawa tensor from FEM zero modes, no Gram/Cholesky, plus KK -> f + H overlaps), `ratio_distribution`.
Certified on three meshes: family-tower first excitation Delta lambda_1 = 0.023 (generic bundle 0.066), Higgs gap 0.55;
the KK-TOWER-01 "19th Higgs mode" retracted as an eigensolver/mesh artefact.  Cross-pipeline gate (slow test): the FEM
tensor reproduces the CC-26-corrected M1 ratio distribution to 1–3% (medians 0.61/0.26), excluding the retracted values
by 6–8 orders.  KK1 couples to fermion + Higgs at Yukawa strength (overlap ratio 0.92).  Register:
outputs/SM/KKT02_CERTIFICATION_AND_FEM_YUKAWA.md.  Base: the 0.31.2 candidate (CC-26 + SM-15 + KK-TOWER-01).  Pin four-way.
