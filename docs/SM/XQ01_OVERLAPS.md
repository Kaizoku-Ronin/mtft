# XQ-01 — what the overlaps X_qq' are made of (2026-10-03)

Module `mtft.matter.overlaps`; tests `tests/test_xq01_overlaps.py`.

1. **Factorisation (EXACT in the PB-01 messenger model).** Messengers are curve copies: a holomorphic differential on X_0(143) times a massive torus mode. Summing the
   curve factors in the propagator gives the Bergman kernel of X_0(143), so X_qq' = (family normalisation) x (Bergman jets at the cusp) x (torus-propagator jets) x
   (condensate components). The four cusps are Atkin–Lehner images of infinity and Atkin–Lehner preserves the Bergman kernel, so the blown-up cusp does not matter.
2. **Family normalisation (EXACT; CC-41).** The jets x^2, x w, w^2 have norms 2 : 1 : 2 on the exceptional curve, so normalised families carry nu = (1, sqrt 2, 1). Coherent
   and parity-even condensates then give m1 m3/m2^2 = 1/4 rather than 1; the data (0.57, 2.55, 0.081) still exclude them, and the rho targets become 0.568 (up),
   0.684 (down), 0.408 (leptons). CP-12's targets assumed unnormalised families and are downgraded.
3. **The curve factor (NUMERICAL).** B_nm = a_n^dagger (Im Omega)^-1 a_m, n, m = 1..6, relative to B_11: diagonal 1, 2.094, 2.926, 4.160, 4.922, 5.208; every off-diagonal
   correlation below 0.32 (B_42: 0.0018). On a disc the Bergman kernel is diagonal in powers of q; the off-diagonal entries measure how far the whole curve breaks the
   cusp's rotation, so X_0(143) nearly conserves that charge.
4. **Charge-conserving pairing (MODEL, conditional).** If each vertex conserves the rotation charge, the charge-n entry carries B_{5-n,n+1} = 0.316, 0.0054, 2.926. The
   charge-1 channel is suppressed about 180-fold, the texture acts parity-even, and the data exclude it. So either the family vertex pairs differently — decided by
   the harmonic representatives of the jet classes — or the torus factor supplies t1^2/(t0 t2) ~ 3e4, or the charge-1 entries arise elsewhere.

Open: harmonic jet representatives; torus-propagator jets; the condensate vertex beyond leading order.
