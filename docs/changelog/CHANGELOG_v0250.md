# mtft v0.25.0 — the Modular Surface Laboratory as `mtft.surface` (2026-09-02)

Import on demand: `import mtft.surface` (not loaded by `import mtft`).

## New subpackage `mtft.surface`
- `manin`  EXACT: CRT canonical P^1(Z/N) (replaces the O(N^2 phi(N)) unit scan),
  Manin S/R/T with R T S = -I, cellular complex, non-definitional gates
  (independent Euler reconstruction, width census vs divisor formula).
  Elliptic levels refused (CC-MSL-02); nu3 uses Kronecker (-3/2) = -1 (CC-MSL-01).
- `cycles` EXACT: tree/cotree integral H_1 basis, Bareiss unimodular gate.
- `hodge`  CERTIFIED/DIAGNOSTIC: intersection form from period-dual Whitney forms;
  Whitney masses, quotient-aware refinement, polar star transported to cycles,
  exact/coexact branch split (the non-vacuous spectral convergence gate).
- `transport` EXACT [gp]: Hecke T_p and AL W_Q on the cycle lattice (route A:
  mspathlog + coset matrices) gated against mshecke/msatkinlehner (route B);
  integrality, charpoly, involutions, commutation, intersection preservation.
- `hodge_structure` CERTIFIED [gp]: J_true from mfsymbol periods over the same
  cycles; Riemann bilinear relations and Hecke/AL commutators at 1e-15; Siegel
  distance of the Whitney family (0.3764 -> 0.2463 -> 0.1755 at N=143, converging,
  rate not certified); elliptic-block E2: j(143a1) recovered to 3.5e-14.
- `gauge`  EXACT/OVERLAY: U(1) flux sectors and theta-series flux sum, flat
  torus U(1)^{2g}, 2D Yang-Mills (U(1), SU(2), SU(3)) closed and with cusp
  holonomies, same-genus control 105 vs 143, Riemann-Roch index, AST
  line-operator census with the theorem gate (squarefree N: cyclic Lagrangians
  = P^1(Z/N), Witten shift = Manin T, theta-orbit sizes = cusp widths).

## Conformal fact recorded
The equilateral Whitney metric on the Manin triangulation is in the conformal
class of X0(N) (Voevodsky-Shabat); the 1-form star is conformally invariant, so
Hodge-theoretic outputs of the refinement family converge to the true ones while
Laplace spectra remain DIAGNOSTIC for the hyperbolic metric.

## Retired wording
"Cusps as stable particles" (printable kit) -> cusps are theta-orbit /
line-operator sectors.  sigma(N) counts line-operator variants, not gauge groups.

## Cross-frame agreement
`mtft.periods.sector_census()` (1,6,5,1) equals the transport census halved.
An explicit Sp(26,Z) intertwiner to the `mtft.homology` canonical frame is the
open item for v0.26.

## Tests
`tests/test_surface.py`: 10 tests (1 GP-gated). Suite: 731+27 prior unchanged.
