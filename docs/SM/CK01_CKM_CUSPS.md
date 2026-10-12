# CK-01 — quark mixing and the cusps of X0(143) (2026-10-02)

Module `mtft.matter.ckm_cusps`; tests `tests/test_ck01_ckm_cusps.py`; correction CC-37 (below).

1. **Inter-cusp distances (EXACT).** Between the height-1 horoballs: d(inf, 0) = log 143 = 4.963, d(inf, 1/11) = log(143^2/13) = 7.361, d(inf, 1/13) = log(143^2/11) = 7.528;
   the other pairs follow by Atkin–Lehner, so any three cusps form the same triangle. Any smooth metric agreeing with the hyperbolic one below height 1 can only
   lengthen them.
2. **The bound.** Families at different cusps mix through the tails of their harmonic jet representatives, which decay at least like e^{-d} (the twist mass only
   steepens it): |V| <= 1/143 = 0.0070 for the closest pair, 13/143^2 and 11/143^2 for the others. |V_us| = 0.225 and |V_cb| = 0.041 both exceed 1/143, so every
   split assignment (1+1+1, 2+1, 1+2) is excluded; all three families share one blown-up point — the Matter-04 background.
3. **The family triplet (EXACT structure).** With exceptional charge 3 the jets split 1 + 2 + 3 by normal degree against 2 + 1 + 0 in the conjugate orientation: the
   three net families are the degree-2 jets Sym^2 of the tangent plane, the monomials x^2, x w, w^2 (x along X0(143), w along 143a1). The torus rotation at the point
   gives them charges 0, 1, 2, broken only by the lattice of 143a1. Torus-constant heavy modes (the ones the D-term mixes with the Higgs) couple only the x^2 family:
   the leading Yukawa matrix is rank one, as in local F-theory Yukawas; torus-dependent components feed the lighter families one power of eps per unit charge.
4. **The texture (STRUCTURAL).** Y_{qq'} = X_{qq'} eps^{q+q'}: masses ~ (eps^4, eps^2, 1), mixing ~ (eps, eps, eps^2). With eps_u = 0.060, eps_d = 0.139 from the second-to-third
   ratios and independent order-one coefficients: m_u/m_t 2.3e-5 [1.0e-5, 5.5e-5] (data 7.4e-6), m_d/m_b 6.3e-4 [2.7e-4, 1.4e-3] (9.4e-4), |V_us| 0.15 [0.07, 0.31] (0.225),
   |V_cb| 0.15 [0.08, 0.26] (0.041), |V_ub| 0.021 [0.011, 0.039] (0.0037). Right pattern; the two small angles come out 4-6 times too large, the known weakness of
   symmetric textures (left- and right-handed jets carry the same charges).
5. **Consequence.** Mixing is local, so the size ratio r is free again: the large r the neutrino sector wants (NU-02) is compatible with the CKM matrix.

**CC-37.** The one-family-per-cusp background of AL-01 is excluded by quark mixing (|V_us|, |V_cb| > 1/143). Its index-level statements stand; it is retired for
phenomenology. The canonical background is again the single blown-up point of Matter-04 (288 fermion zero modes).

**Qualification (EPS-01).** Rank one at leading order holds when the D-term is torus-constant at the point (product-type or single-band gamma-localised condensate); epsilon is the condensate's torus-charge misalignment there.

OPEN: the Taylor coefficients of the actual heavy modes (the X_{qq'}), the source of the up/down asymmetry needed for |V_cb|, the Georgi–Jarlskog problem.
