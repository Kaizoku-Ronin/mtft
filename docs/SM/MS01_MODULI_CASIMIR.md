# MS-01 — can the Casimir energy fix the curve's size? (2026-10-05)

Module `mtft.matter.moduli_casimir`; tests `tests/test_ms01_moduli_casimir.py`.

1. **Structure (EXACT).** With the twist on the torus only, bosons and fermions differ on odd windings, where every degree of freedom counts. Heat-kernel expansion on the curve
   and Poisson resummation on the torus give V_J = -K[(384/pi) N A_X A_E^-3 S0 + 128 a_1 A_E^-2 S1]; in the Einstein frame V_E = alpha/A_X + beta/A_X^2, alpha < 0, beta of sign -a_1.
2. **Condition (EXACT).** A minimum in A_X exists iff a_1 < 0, at A_X/A_E = (2 pi/3)(|a_1|/N)(S1/S0). chi(X_0(143)) = -24 enters a_1; for scalars a_1 = N chi/6 < 0.
3. **Numbers (NUMERICAL).** S1 = 2.107219 (P = 1/2) and 3.238510 (complex pair) — exactly Matter-05's one-loop sums. S1/S0 = 1.326 and 1.081; a scalar count would give
   A_X = 11.1 A_E or 9.06 A_E, the same order as V_cb's 18–36 but short by two to four.
4. **The sign (STRUCTURAL).** For the 8D vector multiplet (six scalar-like components chi/6 each, the curve 1-form -2chi/3, eight spinor components -chi/12 each) a_1 = -chi/3 = +8:
   wrong sign, no minimum, the curve collapses. The matter blocks' flux terms enter a_1 at the same order, so the verdict waits on the full twisted spectrum.
5. **The FI direction.** The blow-up sizes enter only at the next order (a_2): at this order the Casimir energy does not stop GS-1's relaxation either.

Open: a_1 of the full twisted spectrum with flux terms; the a_2 dependence on the blow-ups; anything else that stabilises the moduli.

## MS-02 (2026-10-05): the full twisted spectrum

6. **Dolbeault pairs (EXACT).** Every block's fields are forms on the curve in the block's bundle (degree d), paired (0,0) with (0,1); their traces differ by the index d + 1 - g:
   a_1 = chi/6 + d/2 and -chi/3 - d/2, sum -chi/6. Flux cancels pair by pair. Four complex pairs per block, bosons and fermions alike: a_1/N = -chi/12 = +2 for every block.
7. **No minimum (EXACT at this order).** The size-independent term follows c_1(S~).[omega] = chi(C) A_E - sum eps < 0 for genus >= 2; every blow-up makes it more negative. V has no
   minimum in A_X for any block or twist; the jet vacuum's curve collapses. A sphere would allow a minimum, a genus-13 curve does not. MS-4 superseded (confirmed).
8. **What this leaves (STRUCTURAL).** Nothing in the model stabilises the Kaehler moduli. The consistent reading is a fixed geometry: no GS-1 relaxation, a constant Casimir
   vacuum energy, and A_X/A_E a parameter measured by V_cb unless a canonical prescription fixes it.
