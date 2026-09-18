# mtft v0.32.0 — Standard Model campaign corrections, KK towers, handoff integration, research gates (2026-09-18)

Consolidates the 0.31.2–0.31.4 candidates and the 0.32.0.dev1–dev3 waves on top of PyPI 0.31.1.

## Corrections (register: docs/SM)
- **CC-26** (Astra audit V0311-C01): `hym` family normalisation used (L^-1)^T; fixed, `hym.normalise_yukawa` single implementation,
  `normalisation_is_basis_invariant` gate (tol 1e-6 over measured 1e-8).  Retracted: all SM-08…14 mass-ratio/hierarchy/CKM claims and the
  7.8e-5 "spin-structure invariant" as physics.  Corrected leading order is anarchic (M1 up medians 0.26/0.62), independently confirmed by
  kinetic-orthonormal FEM eigenmodes (KK-TOWER-02, 1–3%).  Frozen M2 file ships raw tensors and all Gram matrices.
- **CC-27** (handoff audit §3): theta-compatible Atkin–Lehner lift on H^0(S0) is Q8 (`arithspin.theta_compatible_al_lift`), D8 section lift retained, Sym^2 intertwiner test.
- **CC-28**: CW-01 register written; SUSY protection of Higgs directions is conditional on an undeclared supersymmetric parent and F-flatness.
- SM-15: degeneracy theorem — a W13-symmetric vacuum implies m_c = m_t; the grading cannot source the hierarchy.
- Review V0311 adopted: relative rank tolerance, rank <= 2 wording, exact cubic/mixed anomaly functions, Majorana obstruction (nu^c has B-L = +1), registers shipped.

## New machinery
- `surface.magnetic`: P1 Peierls FEM Bochner spectra of O(D) in the compact metric (coexact connection from the dual-graph solve with invisible divisor
  lumps); Landau levels realise h^0 (degree 30: 18 modes +0.6%; degree 15: 3); M1 towers (family first excitation 0.023 on two meshes; Higgs gap 0.55);
  FEM Yukawa tensor and KK -> f + H overlaps.  Slow gates included.
- `surface.spin_circle`, `hopf_geometry`, `crt_dessin`, `hodge_blocks` (Wave B, exact; independent-construction tests).
- `mtft.research`: `anomalies`, `charge_lattices`, `mode_operators`, `discrete_anomalies`, `tensor_gs`, `compactification`, `bordism`, `parents`, `pipeline` (Waves B/C).
- `arithspin.al_lift_on_S0` (D8), `smflux.hypercharge_normalisation` (sin^2 theta_W = 3/8), `hym.w13_graded_normalisation` / `w13_epsilon_scan`, `rrspace.cubic_higgs_space` / `mixed_yukawa`.
- Frozen data: `x0143_cubic_higgs.npz`, regenerated `x0143_m2_tensors_h02.npz`.
- Frozen handoff studies (Astra): SC7-01, HOPF-02, AXG-01…04, v0.31.4 audit, with `SHA256SUMS`/`provenance.json`; AXG-04 rerun on import (322/0).
- Legend: 11 new entries (60 registered).  Registers: 24 files in `docs/SM/`.

## Decisions recorded
Gate 4 route 2: one parent, three internal modes; C3X's consistency machinery is the gate battery.  Teaching materials excluded from the package.

## Verification (this build)
Full fast suite 820 passed, 4 skipped, 0 failed (NumPy 2.4 / SciPy 1.17); slow gates run explicitly; clean tree; fresh-venv install; four-way pin 0.32.0.
