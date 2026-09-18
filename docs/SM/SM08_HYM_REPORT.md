# SM-08 — the compact metric, Hermitian–Yang–Mills normalisations, and the first mass ratios (2026-09-13)

## The compact uniformising metric of X0(143) (new; gated)
Liouville's equation ∇²u = e^{2u} − 1 on the cusped surface with cusp asymptotics imposed as the Neumann flux
∂_n u = −2 pi Y0 + 1 on every horocycle (all cusps truncated at Y_c = Y0·width) is solved by Newton on the
SPEC mesh in 5 iterations (residual 1e-11).  GAUSS–BONNET GATE: ∫ e^{2u} dA_hyp = 151.37 vs 48 pi = 150.80
(0.38%; the mesh's own hyperbolic-area error is 0.33%), against 56 pi for the cusped metric — so the solution
is the smooth constant-curvature metric of the compact curve that KK-04 characterised and left uncomputed.
The Neumann flux accounts exactly for the curvature deficit −8 pi + 4/Y0.  e^{2u} ≈ 0.74 in the interior and
e^{−18.7} deep in the cusps.

## Green's functions and HYM metrics
Δ_c G = δ_p − 1/A_c solved for the four CM points (point loads, flux-compatible to 1e-16) and for the cusps 0 and
1/11 (unit Neumann flux through the horocycle), constant mode fixed by ∫ G e^{2u} = 0.  The HYM metric on O(D) is
h_D = exp(4 pi Σ m_p G_p); for the family divisor A + ΣP3 it tames the order-6 cusp poles and the simple CM
poles; the Higgs sector K(2ΣP3) carries y² e^{−2u} for K and e^{8 pi ΣG_P}.  Gram matrices by quadrature
(467,520 points) with the compact area element: Hermitian to 1e-17, positive, conditions 1e4 (family) and 3e6
(Higgs) in the construction bases.

## First mass ratios (DIAGNOSTIC — read the caveats)
With family and Higgs spaces orthonormalised, the up-type mass matrix for a unit Higgs direction v is
y(v) = Σ_k Ỹ_ij^k v_k.  Over 4000 random directions the spectrum is always dominated by one family:
   m2/m3: median 2e-4, 95th percentile 1.4e-3;   m1/m3 below 1e-5 (numerically zero).
The same holds on the exact W13 textures (even directions: rank-2 exactly; odd: 3e-5 median) and under 20%
random perturbations of the Gram matrices (m2/m3 stays below 1e-3).  In the orthonormal family basis a single
family pair carries the coupling; the family Gram's eigenvalues span a factor 1e4, i.e. the three sections have
HYM norms differing by 10^4, and the hierarchy is driven by that spread.
For orientation only: observed m_c/m_t ≈ 2–4e-3 and m_u/m_t ≈ 6e-6 (scale-dependent).  The one-heavy-family
pattern with the lightest family effectively massless is the structure observed; the numerical value of m2/m3
is within an order of magnitude and NOT certified.

## Caveats (each is a concrete next step)
1. Flattened tensor rank: Sym²(3) → 18 shows four singular values above the 18-space precision floor (8.7e-5);
   a genuine rank 4 is geometrically impossible for three base-point-free sections, so the two smallest are at
   the noise floor.  The smallest normalised couplings — and hence m2/m3 — need the 18-space rebuilt at higher
   precision (finer mesh/longer expansions) before any digit is trusted.
2. The point-source Green's functions carry the unresolved logarithmic singularity at the source node; the
   metric factor near the CM points is therefore approximate (a local log correction is straightforward).
3. No mesh-refinement run of the Liouville–Green–Gram pipeline yet; the Gauss–Bonnet gate certifies the metric
   to 0.4%, not the normalisations.
4. The Higgs direction is a modulus: the hierarchy is a statement about the whole moduli space (direction-
   independent), the specific ratios are not predictions.
Tool: `surface.hym` (assemble, solve_liouville, greens_functions, quad_with_nodes, gauss_bonnet_gate); test
`test_surface_hym.py` (Gauss–Bonnet within 1%, Green's function properties).  Down/lepton normalisations need
the same pipeline on the twisted sections (vectorised Gamma_0-only reduction at quadrature points): next.
