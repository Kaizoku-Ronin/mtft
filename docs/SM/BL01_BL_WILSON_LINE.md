# BL-01 — a B–L Wilson line on 143a1 (2026-10-06)

Module `mtft.matter.bl_wilson_line`; tests `tests/test_bl01_bl_wilson_line.py`.

Setup: a flat Z/2 holonomy on the colour stack twists quark bifundamentals by a character chi of the torus lattice and leaves lepton ones alone. Local entries are
untouched; each tail is multiplied by the image-sum ratio R. CI-3 needs lepton tails 1/R = 1.758 times quark tails. Tails: the curve-gradient of DD-01's cubic tail,
summed over torus images, with P2 over a 2-torsion point (tau = 1/2 + 1.02327i, unit area).

1. **Two characters of three kill the quark tails (EXACT).** Images pair at equal distance; a character that is -1 on the doubled offset cancels every pair. One
   character per 2-torsion point survives (R from 0.92 over 1/2, 0.79 over tau/2, falling with the curve separation). Twisting leptons instead only shrinks their tails.
2. **The size is wrong (NUMERICAL).** R = 1/1.758 at s_X = 0.716 (over 1/2) or 0.499 (over tau/2): twisted quark V_cb 0.133 or 0.184 against [0.030, 0.055]
   (factor 2.4–6.1). Where V_cb = 0.055–0.030 the muon factor is 6.6–18 instead of 3.09. Lepton pattern or V_cb, never both.
3. **The vacuum is a saddle (EXACT at one loop, matched spectra).** Every Z/2 point with a character distinct from the Scherk–Schwarz one lies above theta = 0 with a
   negative Hessian eigenvalue (all six pairings; weights |w|^-4, -6, -8). Equal characters make the leptoquark gauginos massless.
4. **By-product (NUMERICAL).** Over a 2-torsion point the untwisted tail is 3.1–3.5 times DD-01's single-image estimate (s_X = 0.85): the V_cb separation moves out ~1.3.
5. **Verdict (STRUCTURAL).** Dressing, insertion and the B–L Wilson line are all closed: the lepton pattern is an open failure of the current matter sector.

CC-54 caveats CI-5. Open: an audit of the flavour failures together (lepton pattern, TO-01's V_cb, GC-01's b–tau normalisation and scale).
