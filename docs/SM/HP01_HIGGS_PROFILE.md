# HP-01 — the light Higgs of the ratio-3 background (2026-10-04)

Module `mtft.matter.higgs_profile`; tests `tests/test_hp01_higgs_profile.py`.

1. **Where (EXACT).** Ratio-3's Higgs entries carry gamma-charges 1 (up, (L, A2)) and 3 (down, (L, B)). Over a curve point x an entry restricts to the torus as the degree-0
   bundle P_{Delta phi(x)}, with cohomology only where Delta phi(x) = O. The modes form a skyscraper over phi^-1(E[Delta]): 4 Delta^2 of them (deg phi = 4). H_u: the four
   cusps, one of them the blown-up cusp. H_d: 36 points over the nine 3-torsion points, four at the cusps. Near the cusp the H_u mode follows the graph w = phi(q),
   phi(q) = q - q^3/3 - q^4/2 - q^5/5 - 2q^7/7 + ... (143a1: a_2 = 0).
2. **Type (EXACT).** The cohomology is H^1 only (chi = -4 Delta^2 - e(e-1)/2, checked against the surface intersection form): A-type modes, no Phi-type global Higgs.
3. **No tree-level Yukawa (EXACT).** The families are A-type jets; Tr Phi ^ [A, A] needs exactly one Phi-type field, so Q u^c H vanishes whatever the profiles, and PB-01's
   mixing also needs a Phi-type Higgs. FS-4 and RT-1 are superseded (CC-43); RT-2 and FS-3's eps-ratio reading become caveats.
4. **Where it works (EXACT).** The coupling Higgs carries exceptional charge +6, so with a gamma-neutral bidoublet its Phi-type modes are pullbacks: torus-constant times
   holomorphic sections of K_C (x) L_a on X_0(143). Their Taylor data at the cusp are q-expansion coefficients — the HJ-01 structure. FS-01's minimal family qualifies.

Open: in a gamma-neutral background, which Phi-type Higgs combination is light and takes the vev, and its q-expansion at the cusp.
