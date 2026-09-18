# SM-03 — the up-type Yukawa tensor of M1 on X0(143) (2026-09-13)

Model M1 (SM-01): Q in (c, L-bar) and u^c in (c-bar, a) with internal bundles S0 (x) O(±(P1+P2+P3)); the Higgs
H_u in the (L, a-bar) sector as a 4D scalar from the internal gauge field, internal zero modes in
H^1(O(−2 SigmaP)) = H^0(K(2 SigmaP))^*.  The coupling Q_i u^c_j H_k is the multiplication
H^0(S0(SigmaP)) (x) H^0(S0(SigmaP)) -> H^0(K(2 SigmaP)), i.e. a 3 x 3 x 18 tensor.

## Construction (all point conditions at explicitly known CM points; DIAGNOSTIC 1e-10 precision, exact structure)
- Three-family sections phi_i = G_i/(f1 f_K): G_i in H^0(K^2) with cusp orders >= 8 at 0, 1/11 (exact) and simple
  vanishing at the 20 W143 points and P4: rank 33 of 36, dimension 3 (v0.30.1).
- Higgs target psi_k = G'_k/f1^2: G'_k in H^0(K^3) vanishing to ORDER TWO at the same 21 points (42 conditions,
  derivatives by the chain rule through the composite reduction matrix, checked against finite differences to
  4e-8): rank 42 of 60, dimension 18 = h^0(K(2 SigmaP)) (Riemann–Roch: 24 + 6 + 1 − 13).
- f_K = g1^2 w with w = eta(tau)^2/eta(13 tau)^2 (pole of order 1 at oo, offset −1): f_K = q − 2q^2 − q^3 + …,
  cuspidal divisor 2A.
- Y: G_i G_j = sum_k Y_ij^k f_K G'_k solved on 131 q-coefficients; maximal relative residual 3.2e-12.
  This residual is the E2 of the wave: the products of the 3-space land in f_K times the 18-space, as the
  divisor bookkeeping requires, although the two spaces were built independently.

## Results
- Y is symmetric in (i, j); for every Higgs direction tested (random and the 18 basis directions) the 3 x 3 mass
  matrix has rank 3: all three families become massive, no texture zero in these bases.
- Not yet: the W13 grading of the three sections and of the Higgs directions (a Z2 flavour symmetry that
  should force a texture), the compact-metric normalisation of the sections (kinetic terms) needed to turn Y
  into mass ratios, the down-type and lepton tensors (same construction with the (c-bar, b), (L-bar, d),
  (b, d-bar) sectors), and the choice of Higgs direction (potential).  The tensor is basis-dependent until the
  normalisations are in; its rank properties are not.
Tool: `surface.rrspace` (`higgs_target_space`, `up_yukawa_M1`, `f_K_series`, `eta_quotient_series`,
`eigen_vals_and_derivs`), tests in `test_surface_rrspace.py`.  Frozen: `x0143_weight1_g1.json`.
