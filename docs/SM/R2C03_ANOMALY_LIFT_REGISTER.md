# R2C-03 — kappa = 1; tensor integrality for the p2 term; the colour-cubic obstruction of M1's declared content (2026-09-19; EXACT)

## kappa = 1
6D -> 4D reduction with D = d - iA, kinetic -(1/2 g6^2) tr F^2, canonical complex Higgs h = 2 a_zbar/g6 (from |a5|^2 + |a6|^2 = 2(|a_z|^2 + |a_zbar|^2)),
vertex a_zbar (sigma5 - i sigma6) = 2 a_zbar sigma^- with unit matrix element between S^+ and S^- zero-mode components: y = g6 I = g4 sqrt(A) I exactly.
Hence y_t/g_4 = s_max(v) = 0.42 / 0.63 / 0.87 (5/50/95% over Higgs directions), observed ~0.7 at ~1e16 GeV.  `vector_higgs.reduction_constants`.

## p2 cancellation by chiral tensors
[p2]: Weyl -1/1440 per signed dimension, self-dual tensor +1/360, gravitino -49/288.  Need n_T - n_T' = n_grav/4: of the four ledger-preserving
assignments, (eps_cd, eps_ab) = (+,+) with n_grav = 24 (6 net tensors) and (-,-) with 16 (4) survive; (+,-) 22 and (-,+) 18 are excluded by
integrality; a gravitino cannot be balanced by tensors (269/4).  `gravitational_anomaly.m1_tensor_survivors`.

## The colour-cubic obstruction (2-form x 6-form)
dI8/dD3 = -f_L - f_a/2 - f_b/2 + (5/2) f_c - f_d/2 for the family signs.  The f_L term is sourced by the Q sector alone — the only coloured field with
U(1)_L charge among Hom-type bifundamentals (i, j-bar) — and no additional (c, d-bar) content changes it.  As a 2-form x 6-form it is not a product of
4-forms, so no Green–Schwarz term cancels it.  Consequence: M1's declared stack group U(3) x U(2) x U(1)^3 with flux (0,-3,3,3,0) has NO anomaly-free
6D completion preserving the three families within Hom-type content.  This sharpens AXG-02's finding to a closed statement with its two escape routes:
(a) non-Hom sectors (N_i, N_j) with eps = +1, each cancellation adding three vector-like coloured pairs; (b) an L-stack without its own U(1) and the
flux vector (3,0,0,0,-3), which keeps all six family indices at +-3 and turns both Higgs sectors into index-0 sectors (h^0 = h^1 = 2 for the
trivial twist).  Either is a new model record for the route-2 search, not an extension of M1.  Gate `six_d_anomaly_lift` updated (EXACT witness).
