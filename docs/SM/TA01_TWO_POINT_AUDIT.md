# TA-01 — auditing the two-point design (2026-10-05)

Module `mtft.matter.two_point_audit`; tests `tests/test_ta01_two_point_audit.py`.

1. **Fibre geometry (NUMERICAL).** Unit-area torus of 143a1: deep hole 0.626 from O; farthest torsion points (3-torsion) 0.598; 2-torsion 0.494, 0.563.
2. **Constraints (EXACT) and target (ESTIMATE).** D-flatness 2 eps1 + eps2 = a_Q A_E; on a shared fibre eps1 + eps2 < A_E. The top needs eps2 ~ 0.55 A_E (radius 0.418), giving
   eps1 = 0.225 (radius 0.268) for a_Q = 1; at the low end eps2 = A_E/3 both radii are 0.326.
3. **Mixing (ESTIMATE).** Degree-k jet tails fall like (r/d)^(3+k). At separation 0.598 the top's tail at the jets is 0.34 (or 0.16): V_cb ~ 0.16-0.34 against 0.041, and
   m_c/m_t ~ its square against 0.0036. The same-fibre design fails by 4-8x; TP-3 becomes a caveat.
4. **Where the top must sit (STRUCTURAL).** V_cb = 0.041 needs a separation of about 2.9 top-radii (0.95-1.21), beyond the fibre. The top's point must be on another curve
   fibre, at a non-cusp point (another cusp is too far, CC-37). Candidates: special points of X_0(143), e.g. preimages of torsion points of E. The jets' three units are
   odd, so with torus parity at least one must be curve-direction.

Open: a natural non-cusp point at the right distance; PA-01 stability with two points; the condensate's global fraction.
