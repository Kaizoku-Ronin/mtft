# mtft v0.29.1 — KK05 audit corrections; W11 parity theorem recorded (2026-09-12)

- Wilson lines: lifting condition corrected to "the charged line F^q is nontrivial" (order-q lines act
  trivially on charge q); test with an order-3 line on Gamma_0(11).
- `condensation_energy` now returns split / balanced / drop in the stated convention E = (1/2g²)∫tr B²:
  (0,72) → 216 pi, 108 pi, drop 108 pi (per g² R²); the earlier scalar was the split energy.
- `petersson_gate`: the cross-sector check now inspects the computed Gram (max correlation ~2e-3) instead
  of the prediction's own mask; docstring states the cusped-metric scope and the e^{2(1−k)u} scaling for
  compact-metric norms of K^k sections.
- Recorded (exact, from frozen AL traces): W11 has 0 fixed points, so every W11-equivariant bundle on
  X0(143) is a pullback from the genus-7 quotient and has even degree; net chirality chi(S⊗E) = deg E is
  therefore even in any W11-symmetric flux background — three families require breaking W11 (W13 with 4
  fixed points and W143 with 20 admit odd-degree equivariant bundles).  This is the structural change
  KK08 asked for.  No SM derivation is claimed; KK06/KK07/KK08 results are recorded as Astra reported them.
Base: live 0.29.0 (SHA-256 a5a0d5e0…).  Pin four-way.
