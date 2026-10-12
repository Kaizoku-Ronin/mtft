# LB-01 — the lattice of 143a1 at the blow-up scale, and the top (2026-10-04)

Module `mtft.matter.lattice_breaking`; tests `tests/test_lb01_lattice_breaking.py`.

1. **Parity (EXACT).** j(143a1) = -262144/1859, so the lattice's only rotations are +-1 and torus characters are conserved mod 2. Of YA-3's channels, the torus
   mismatch is 6, 6 and 7; the odd one dies, the Higgs's slope drops out, and Y3 ~ h0 c2 c0*, Y4 ~ h0 (c1 c0* + c2 c1*) both need six units of breaking.
2. **Multipoles (NUMERICAL).** Unit-area lattice: G4 = 1.4016, G6 = 3.8595. The first shell (angles 0, +-64, +-116, 180 degrees) has six-fold sum 5.66 against a
   hexagonal 6. The torus Green's function carries Re(G_2k w^2k)/(2k), so the six-unit breaking across a disc of radius r is |G6| r^6/6 at the edge, at most |G6| (2r)^6/6.
3. **Size (ESTIMATE).** Exceptional area A_E/3 gives radius 0.326: six-unit breaking 7.7e-4 to 0.049 (0.17 at A_E/2, 0.39 at 2A_E/3). With the condensate tilt the
   heaviest jet entry is at most 0.0095 for t = 0.14.
4. **Verdict (ESTIMATE).** y_t at the b–tau scale (1.1e6 GeV) is 0.70, about 73 times the jet bound; y_b there is 0.0104, at the bound. The first two families fit jet
   couplings; the top cannot be a jet in backgrounds of this class.

Open: a background whose third family comes from an unsuppressed, non-jet channel.
