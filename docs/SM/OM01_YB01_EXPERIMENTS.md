# OM-01 and YB-01 — the two experiments (2026-10-08)

Modules `mtft.matter.odd_modes`, `mtft.matter.bridged_top`; tests `tests/test_om01_yb01_experiments.py`.

## OM-01: forced odd modes (the Artebani–Quispe idea, transplanted)

1. **The count (EXACT).** Under w -> -w, O(N O) has 1 + floor(N/2) even and floor((N-1)/2) odd sections; the Lefschetz sum over the four 2-torsion points gives
   #even - #odd = ((-1)^N + 3)/2. Odd modes are forced from N = 3 up (N = 3: the odd one is wp').
2. **No torus degree anywhere (EXACT).** One unit of torus flux makes 24 chiral families (TP-01), so every block has torus degree 0; gamma charges are Poincare-type,
   also degree 0. Zero modes are constant along the torus — even. The only forced odd modes are the w-jet families: states, not sources.
3. **Verdict (EXACT).** No down-sector field that could carry an odd vev or coupling has an odd mode. An odd B-stack mode would cost 72 chiral families. PA-1's
   verdict is topological.

## YB-01: the bridged top (pre-registered threshold y_t >= 0.5)

1. **The bridge (ESTIMATE).** Bridged/local = 4 (r1/r2)^2 t(s) = 0.41–0.45 at the Fricke separation (tail t = 0.25–0.28).
2. **The test.** y_t <= 0.29–0.32 (flat tail, marginal local 0.70), 0.17–0.27 with curvature: **fail** by 1.6 or more.
3. **Verdict (STRUCTURAL).** Both redesigns are closed: this two-point design does not carry the quark flavour structure. The arithmetic results stand.

CC-60 caveats PA-4.
