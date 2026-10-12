# DD-01 — a parity selection rule for the CKM matrix (2026-10-05)

Module `mtft.matter.double_duty`; tests `tests/test_dd01_double_duty.py`.

1. **Parities (EXACT).** The torus involution w -> -w fixes the jet point and both of the top's homes (2-torsion points are their own negatives). The constant-jet top and the
   curve-direction jet x are even; the torus-direction jet w is odd.
2. **A decoupled family (EXACT).** If the background respects the involution, every entry between opposite parities vanishes in both sectors: the w-family mixes with nothing,
   and the CKM matrix has a decoupled row and column. Excluded by data, so the involution must be broken.
3. **The hierarchy (STRUCTURAL).** Breaking of size eps_P with the odd family first: V_us ~ eps_P, V_cb from the top's tail (geometry), V_ub ~ c eps_P V_cb. The data need
   eps_P ~ 0.225 and c = 0.40, a natural O(1): Wolfenstein's lambda, lambda^2, lambda^3 read as parity breaking, geometry and their product.
4. **The source (STRUCTURAL).** The gamma-flux is odd under the involution, so gamma-charged matter blocks (TH-3's design) supply eps_P: the Cabibbo angle measures their
   gamma-charge at the jet point.
5. **SP-2 refined (ESTIMATE).** Through the tail's gradient, V_cb ~ 3 (r1/d)(r2/d)^3: separation 0.95-1.09, curve area over a 2-torsion point 18-28 A_E.

Open: the O(1) coefficient c from the overlaps; eps_P from the gamma-flux at the jet point.
