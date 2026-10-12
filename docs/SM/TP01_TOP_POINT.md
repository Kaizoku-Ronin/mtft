# TP-01 — the flux family fails, the constant-jet top survives (2026-10-05)

Module `mtft.matter.top_point`; tests `tests/test_tp01_top_point.py`.

1. **Flux families come in 24s (EXACT).** chi(V*) - chi(V) = 24 b on the surface, so one unit of torus flux gives 2g - 2 = 24 chiral families (one A-type, 23 Phi-type).
   TF-3 is retracted (CC-45). Chiral families come only from exceptional charges.
2. **Two units bought globally (ESTIMATE).** The condensate's global components carry Taylor coefficients of every order; the conj of the second supplies character -2 at
   (r/L)^2 to (2r/L)^2: 0.11-0.42 at eps = A_E/3, 0.70 at 0.55 A_E. A degree-0 top at a point of exceptional charge 1 is marginal, not excluded (CC-46 caveats TF-2, LB-4).
3. **Two-point design (STRUCTURAL).** (cusp, O) with charge 2: two degree-1 jet families; (cusp, P), P a torsion point on the same fibre: one constant-jet family. The
   coupling Higgs has exceptional charges (+4, +2) and stays a torus-constant Phi-type pullback. The jets need 3 units, the top 2 (paid globally).
4. **A unique Higgs (EXACT).** H^0(C, K_C (x) L_H) with deg 24 - (a_Q + a_u): thirteen cusp forms at 0 (trivial L_H; twelve for a generic flat one), exactly one, constant on the curve, when a_Q + a_u = 24 and
   L_Q (x) L_u = K_C. With b = 0 the curve degrees are free; D-flatness of the Q block then ties a_Q A_E to 2 eps_1 + eps_2.

Open: jet–top mixing across the torus (V_cb ~ 0.04), the second family's three units, PA-01 stability with two points, the condensate's global fraction.
