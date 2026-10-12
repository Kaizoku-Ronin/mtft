# TF-01 — can any jet configuration host the top? (2026-10-05)

Module `mtft.matter.top_family`; tests `tests/test_tf01_top_family.py`.

1. **Minimal-breaking theorem (EXACT, local characters).** Net jet families of a block with exceptional charge e are the degree-(|e| - 1) jets, character |e| - 1. The
   condensate lives in u^c's block, so its density reaches character -(|e_u| - 1) at best, and a quark Yukawa needs u = |e_Q| + 1 >= 2 units of rotation breaking
   whatever the degrees: 4 at the standard point, 2 for a constant-jet top at an |e| = 1 point. No jet configuration is unsuppressed.
2. **Two units (NUMERICAL).** E2-hat(tau) = 0.105 for 143a1 (zero for square and hexagonal lattices): the two-unit breaking at the blow-up scale is 0.019-0.075, so the
   best jet configuration gives y <~ 0.08 against y_t = 0.70 — nine times short. With LB-01 this closes every jet realisation of the top. CK1-2 becomes a caveat.
3. **The two-plus-one design (STRUCTURAL).** Two families as degree-1 jets (|e| = 2, the spin-1/2 doublet x, w) and a third flux family: global chiral index a b, so
   (a, b) = (1, 1) on Q and (-1, -1) on u^c; then b_H = 0 keeps the Higgs torus-constant and Phi-type, and the flux family's Yukawa is a global magnetised-torus
   overlap with no character obstruction. The jet families need 3 units, the right ballpark for the second family.

Open: build and audit the two-plus-one background (anomalies, D-flatness with torus flux, PA-01 stability, jet–flux mixing).
