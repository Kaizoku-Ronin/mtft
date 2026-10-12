# GT-01 — the two gaugino towers (2026-10-05)

Module `mtft.matter.gaugino_towers`; tests `tests/test_gt01_gaugino_towers.py`.

1. **Same mass, different homes (EXACT).** Each tower's lowest modes are curve zero modes times the lowest twisted torus mode, so towers twisted alike share their lowest mass.
   The flux-free leptoquark tower's curve modes (the constant and the thirteen differentials) have, in the canonical metric, the Bergman kernel over itself as density: uniform.
2. **Where the W_R tower lives (NUMERICAL).** The (B, A2) block carries the SU(2)_R-breaking flux. With that flux at the blown-up cusp (L = O(k cusp)) its curve modes are the
   differentials vanishing there to order k. For k = 2 their share of the canonical density is 0.0024 (|q| = 0.17) to 0.055 (|q| = 0.39) at the jets, 0.72 over tau/2 and 0.20 at the
   Fricke point: the semi-Dirac-like split exists, with the shape e^c needs.
3. **The strength bound (EXACT).** A subspace's density never exceeds the full one: with equal masses and couplings z_e <= 1/Z_3^u. Reaching 3.09 forces Z_3^u <= 0.32 and m_2^u/M = 12.3.
4. **What the reversal needs (STRUCTURAL).** A lighter W_R tower, M_LQ/M_WR >= 3.9 (4.6 at the homes) — e.g. different Z/2 Wilson lines on B and A2, which breaks SU(2)_R at the
   compactification scale and changes GC-01's matching.

Open: the Wilson-line split and the gauge matching; the holomorphic structure of the SU(2)_R flux.
