# YA-01 — a character audit of the Yukawa amplitude for jet families (2026-10-04)

Module `mtft.matter.yukawa_audit`; tests `tests/test_ya01_yukawa_audit.py`.

1. **The overall rotation (EXACT, local model).** The overall U(1) (x, w) -> e^{i alpha}(x, w) is a symmetry near the blown-up point in every background we use (the
   gamma-flux has character (1, -1)). Characters as forms: each family jet +2; a Phi-type Higgs at Taylor order m, m + 2; the condensate density from jets of degree
   <= 2, between -2 and 2. Any amplitude with the Higgs entering holomorphically — tree level (YK-01) or PB-01's D-term mixing, int D conj(phi_X) phi_H — has
   character 4 + (m + 2) + d >= 4 and vanishes. The minimal gamma-neutral background, whose Higgs is Phi-type, is covered.
2. **Scales (EXACT given the backgrounds).** The exceptional area is A_E/3 (linear size 0.58 of the torus): the torus rotation is broken by 143a1's lattice at the
   blow-up scale. The curve is larger and nearly rotation-symmetric (XQ-3).
3. **Curve charge alone (EXACT under 2).** Only n = 3, 4 survive: texture [[0,0,0],[0,0,Y3],[0,Y3,Y4]], rank 2. Y3 comes from the Higgs's value at the cusp times the
   condensate's c2 c0* component, Y4 from its value and slope times the tilted components. An untilted condensate gives nothing; a coherent tilt t gives m2/m3 ~ t^2.
4. **Consequence (STRUCTURAL).** Jet-family Yukawas are set by the breaking of the point's rotations by the global geometry — mostly the torus lattice at the blow-up
   scale — and need the condensate's tilt. CC-44 caveats PB1-2, EPS1-2, HJ-2 and HJ-4.

Open: the torus lattice's breaking at the blow-up scale, and whether it gives an O(1) top Yukawa.
