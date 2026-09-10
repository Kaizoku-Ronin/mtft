# mtft v0.27.2 — periods-frame integer map: three frames, one lattice (2026-09-10)

The frozen v6 period basis (`periods/_data/X0_143_period_basis_v6.json`, 26 vectors in PARI's
msinit Manin-generator coordinates) and the surface cycle basis span the same sublattice:
C V = K_period exactly with V unimodular.  X = S^-1 V^-1 (S = `periods.symplectic_change()`)
is frozen as `X_periods` in `surface/_data` with five load-time gates: X conjugates the frozen
cycle-frame W11/W13 exactly onto `periods_frame_ops`, X^-T Jint X^-1 = -J_std, and
X J_true X^-1 = -hodge_complex_structure to 2.5e-15.  Together with Pi (v0.27.0) the cycle,
canonical (mtft.homology) and periods (mtft.periods) frames are one lattice with exact maps
(28 frozen gates).  The full 26 x 26 Hodge structure derived from the tree/cotree cycles by
mfsymbol agrees with the tau0 pipeline — the direct E2 the j-invariants only sampled.
Accessor: `intertwiner.periods_frame_map()`.  Tests: `tests/test_surface_periods_frame.py` (2).
Built on the live 0.27.1 tarball (pull-before-branching).  Release pin bumped four-way.
