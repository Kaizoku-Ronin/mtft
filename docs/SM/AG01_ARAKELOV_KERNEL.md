# AG-01 — the canonical Green's-function kernel of X0(143) (2026-10-07)

Module `mtft.matter.arakelov_kernel` (frozen numbers); provenance `scripts/arakelov_kernel_x0143.py` (needs cypari2); tests `tests/test_ag01_arakelov_kernel.py`.

1. **The kernel is conformal (EXACT).** Green's functions for different volume forms differ by f(P) + f(Q) + C, so every mixed second derivative coincides:
   d_P d_Q g = B(P,Q) - pi w(P)^T Y^-1 w(Q). The canonical metric enters only through the local normalisation sqrt(mu(P) mu(Q)).
2. **The second-kind differential (EXACT construction).** No weight-2 eta quotient has a lone double pole at the cusp; the smallest modular function with poles only
   there has order 60. It is eta(t)eta(11t)/(eta(13t)eta(143t)) (orders -6, 6, 6, -6) times a weight-2 cusp form vanishing to order 5 at infinity and 7 at 1/13.
   Residue zero to 3e-14. Holomorphic part fixed by Petersson orthogonality over the Ford domain of the Atkin–Lehner quotient (area 14 pi to 3e-6; reproduces
   Im Omega to 8e-6). Eta factors are evaluated as products, since the q-series of a weakly holomorphic form is numerically unusable at low height.
3. **The kernel (NUMERICAL).** K from the blown-up cusp (area-1 units): Fricke point 74.55, tau/2 preimage 22.87, cusp 0 13.02, cusps 1/11 and 1/13 1.73 and 1.27;
   construction-independent to 0.8%. Against the flat 1/l^2: 0.85, 0.73, 0.60.
4. **The pre-registered test (NUMERICAL).** r = 2 pi r1 r2^3 K / V_cb, V_cb in [0.030, 0.055], against the Cabibbo window 91–896:

| Home | r range | central r | Outcome |
|---|---|---|---|
| Fricke point (CM) | 167–306 | 224 | pass |
| tau/2 preimage | 51–94 | 69 | edge only |
| cusp 0 | 29–53 | 39 | fail |
| cusp 1/11 | 3.9–7.1 | 5.2 | fail |
| cusp 1/13 | 2.8–5.2 | 3.8 | fail |

5. **Reading (STRUCTURAL).** One pass in four arithmetic homes, at a point CA-3 had flagged post hoc: modest evidence. It fixes the top's home at the Fricke point
   and r near 170–300, consistent with V_us (CA-1) and NU-02.

CC-56 caveats FA-2.
