# Unreleased — PV-01: parent actions and vacua (2026-09-29)

Candidate for v0.34.0 once the release gate is run on Roger's machine (the ledger's release shape: `parent` gates plus one B1 toolchain).

- `research/anomaly8d.py`: the ten-form anomaly of chiral 8D content, exactly (lambda-ring Chern characters; tr_{Lambda^2} F^5 = (N - 16) tr F^5
  + 5 tr F tr F^4 + 10 tr F^2 tr F^3; SU(N) quintic-free condition n_F + (N - 16) n_A + (N + 16) n_S = 0 with the remainder s_3 X_4 — one
  Green–Schwarz 2-form, an axion in addition for U(N); no exact cancellation without GS); `classify`, `gs_decomposition`, `specialise`
  (explicit Chern roots, SU(N) tracelessness).  Split-stack parity lemma and the B1 finite search (`b1_content_blocks`, exact rational
  `hypercharge_assignment` — the ledger's row-3 tool — `sm_net_content`, `b1_scan`, `b1_trilemma`).  Result: U(8) = U(3) x U(2) x U(1)^3 with
  Lambda^2(8) + 8 x 8 gives exactly three Standard-Model families on split backgrounds (60 in box 3, 156 in box 4; none with one or two singlet
  stacks), but never a Yukawa, a slope-free Higgs and a colour-matched slope at one ratio (best: r_H = 2 with r_c = 5/2).
- `research/vacuum_energy.py`: the flux energy of a split background as an exact function of the shape r = A_X/A_E, its AM–GM floor per line
  bundle, the Cauchy–Schwarz excess, the Hermitian–Yang–Mills floor of the topological class, and the Einstein-frame scaling (no volume
  stationary point; shape stabilised at fixed volume).  S1 and S2 share E = 4 pi^2 (19/r + 11 r); floors 32/3 at r = 1 and r = 7.
- `mtft.parent`: `ParentCandidate`, `gates`, `ledger_table`, a catalogue (three-stack S1/S2, three B1 U(8) exemplars, M3, C3X) and the CLI
  subcommand `python -m mtft parent [ID]`; the three-stack column reproduces the ledger's status column.
- Docs: `docs/SM/PV01_PARENT_VACUUM_REGISTER.md`; `PARENT_ACTION_REQUIREMENTS.md` section E (append-only); Legend entries
  `anomaly_8d_ten_form`, `b1_split_stack_search`, `flux_energy_landscape`, `parent_ledger_gates`; `tests/test_pv01.py` (13 fast, 2 slow scans).
- Version unchanged (0.33.1) until the release is chosen.  No historical numerical results or frozen artifacts change.

## PV-02 (2026-09-29): supersymmetry named in the code, and the blow-up (B3)
- `surface/blowup.py`: NS lattice of Bl_p(X0(143) x 143a1) (U (+) <-8> (+) <-1>, modular degree 4), K~ = 24 f1 + E, not spin; the supersymmetric
  topological twist (net = e - 24 b, a coboundary phi_i - phi_j for adjoint blocks), ampleness 0 < eps < min(A_X, A_E), slopes with the exceptional
  modulus, `polystable_locus`, the exact-SM phi conditions (n_A = n_B + 1, minimal U(8)), `d_flat_sm_search` (116 D-flat exact-SM assignments in
  the box, fewest extras 31) and the closed-form `collinear_family` (traceless SU(8) flux along v = (1, 0, -3), D-flat on eps = A_E/3 for every
  r > 1/3, three families per unit charge).
- `tests/test_pv02.py`: the lattice, chi(O) = 0 and RR against the product, the SUSY-twist chirality -24 b, the twin-prime uniqueness of 143
  (g - 1 = (p + q)/2 only for k = 2), the D-term identity E_split - E_HYM = slope variance / r, the five means of (11, 13), the phi conditions,
  the search and the polystable locus.
- Docs: `docs/SM/PV02_BLOWUP_SUSY_REGISTER.md`; ledger section F; Legend `blowup_susy_twist`; `viz/arithmetic_susy_picture.html`.

## CC-35 (2026-09-29): reflection positivity of the MTFT lattice action retracted from PROVEN to OPEN
- `research/reflection_positivity.py`: power sums as virtual characters (Murnaghan–Nakayama `hook_expansion`, `power_sum_characters`), exact
  first-order kernel coefficients, Weyl-integration `kernel_coefficients` for SU(2)/SU(3), the site-reflection pairing, positivity brackets, the
  link-reflection note, and the character repair (`weights="sym"`).
- Finding: the lemma "c >= 0 makes exp(c Re Tr P^n) positive-definite" (Paper 13 §4.2, Paper 24 Thm 2.4) is false for n >= 2; with w_1 = 0 nothing
  compensates; at the package defaults (kappa = 1, y = 0.18174) the Polyakov kernel has negative 3bar (-0.004698), adjoint (-0.003969) and SU(2)
  spin-1/2 (-0.005944) coefficients, so the papers' site-reflection positivity fails at strong coupling.  Positive-definite only for kappa > ~666
  (SU(3)) / ~59.2 (SU(2)) at y_c over the tested representations.  Link reflections and the gauge-invariant sector: open.
- Wording corrected: mu_N(y) is the classical stiffness of the Polyakov potential, not a spectral gap (`lattice.py`, `arithmetic.py`, `tower.py`,
  `quantum.py`, Legend `mu_stiffness`, `filtered_moment_identity`); the Polyakov term is recorded as breaking O(4)/hypercubic invariance.
  Function names unchanged.
- Docs: `docs/SM/CC35_REFLECTION_POSITIVITY.md` with `figs/cc35_polyakov_kernel.png`; Legend `cc35_reflection_positivity`; `tests/test_cc35.py` (7 tests).
- Papers flagged for the next revision: 13 §4.2; 24 abstract, Thm 2.4, Steps 3–4, Thm 2.5 status; 5 Thm 7.8 wording; Dictionary line for Paper 24.
