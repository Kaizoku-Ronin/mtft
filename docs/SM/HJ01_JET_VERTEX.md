# HJ-01 — the family vertex of the jet classes (2026-10-04)

Module `mtft.matter.jet_vertex`; tests `tests/test_hj01_jet_vertex.py`.

1. **Selection rule (EXACT, local model).** Near the blown-up point the geometry is rotation-symmetric in the curve direction. The jet classes come from holomorphic
   jets through an equivariant coboundary, so family q carries curve charge 2 - q. A messenger whose curve factor is a holomorphic differential — the curve copies of
   PA-01/PB-01, whose sum is the Bergman kernel — carries beta + 1 >= 1. The vertex needs total charge zero, and (4 - n) + (beta + 1) is always positive: the curve
   copies do not couple at leading order. Messengers that couple must depend antiholomorphically on the curve, i.e. be massive along it.
2. **Locality (STRUCTURAL).** Both vertices sit at the same point, so the coupling messengers' propagator is dominated by the blow-up scale, where the symmetry holds:
   local and charge-conserving at leading order. X_0(143)'s global data enter through the light Higgs's profile — its Taylor (q-expansion) coefficients at the cusp.
3. **The round blow-up gives no hierarchy (EXACT).** With the full SU(2), a coherent condensate supplies only the m = 0 spin-2 component along its axis; the texture has
   singular values 2 : 1 : 1 for every tilt. The hierarchy has to come from the curve/torus anisotropy and the Higgs's Taylor data.
4. **Leading texture (MODEL).** Y_qq' ~ nu_q nu_q' h_(4-n) t_n. Illustrative curve log-convexities h4 h2/h3^2: 0.25 (143a1 newform, a5 a3/a4^2), 0.43 (profile peaked at the cusp).

Corrections: CC-42 retracts XQ-4; XQ-1 becomes a caveat (it describes the forbidden curve-copy channel).
Open: the light Higgs's curve and torus profiles in the ratio-3 background; the torus Taylor data; the order offset; the condensate tilt's contribution.
