# PAR-01 — the parity audit (2026-10-08)

Filed as PA-01, which is Path A's wave label; relabelled PAR-01 by CC-62. Claim ids PA-1 to PA-4 are unchanged.

Module `mtft.matter.parity_audit`; tests `tests/test_pa01_parity_audit.py`.

| Odd source | First-order reach | Effect on 1–2 mixing |
|---|---|---|
| gamma mixed flux (Poincare type) | every gamma-charged block, the condensate | basis rotation — cancels in U_u^dagger U_d (TL-2); physical only via the condensate |
| condensate's odd components | u^c, e^c, nu^c, L, H_u | real, but up sector and leptons only (DS-1) |
| P2's torus displacement; 13 odd moduli H^1(O_X) (x) H^0(T_E) | the third families | tails only, ~1e-5 at the Fricke separation (FH-2) |
| torus Wilson lines | gamma-charged blocks | continuous ones sit at even values (VM-01); Z/2 ones are saddles (BL-3) |
| even | curve fluxes, exceptional charges, torus flux densities, lattice six-fold and two-unit breaking, even spin structures, Phi-type Higgs modes | — |

1. **No odd source reaches Q, d^c or H_d at first order with a physical effect (EXACT).** DS-2's floor applies.
2. **The one loophole (NUMERICAL).** An up-sector local arm large enough to lower c needs |v_u| >= 0.18 (z = 0.25) and makes m_u/m_c 24–56 times too large.
3. **The 13 odd deformations (EXACT).** 37 even, 13 odd; the odd ones move the torus origin along the curve and reach the jet point only through tails.
4. **Verdict (STRUCTURAL).** |V_ub| >= 0.0084 with tails capped at 13.8: the two-point design is ruled out as a flavour model. Redesigns: a B-stack field with an
   odd profile (absent from the spectrum) or a different home for the left-handed top.
