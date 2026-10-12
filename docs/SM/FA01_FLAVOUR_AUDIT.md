# FA-01 — the flavour audit (2026-10-06)

Module `mtft.matter.flavour_audit`; tests `tests/test_fa01_flavour_audit.py`.

| Observable | Status | How | Ledger |
|---|---|---|---|
| V_us = 0.225 | fit | eps_P from the gamma-tilt; needs r = 91–896 | CA-1, DD-01 |
| V_cb = 0.041 | fit | sets the P1–P2 separation and so r | DD-4, SP-2 |
| V_ub | fit | c = 0.40, consistent, not sharp | CA-4 |
| m_s/m_b, m_c/m_t with V_cb | failed (pre-registered) | symmetric tails miss by 1.9 | TO-01, FL-01, CD-01 |
| lepton pattern 3.09, 0.30 | failed | dressing, insertion, B–L Wilson line closed | QL-01, CI-01, BL-01 |
| y_tau/y_b = 1.40–1.53 | open | needs a 30–35% threshold | GC-01 |
| Pati–Salam scale R = 3.68 | open | flux designs quantise R | GC-01, TH-01 |
| two measurements of r | post-hoc pass | single-image artefact, see 2 | CA-2, FA-2 |

1. **Scoreboard (STRUCTURAL).** Eight observables, six continuous parameters introduced to meet them, two pre-registered tests (no hit, failed).
2. **The tail model was single-image (NUMERICAL).** On the flat torus the tail is an image sum, crossing over to 2D (V_cb -> 2 pi r1 r2^3 / s_X^2) beyond the torus size.
   V_cb's r rises by 2.5–3.4: 63–130 over the 2-torsion points (single image 22–34), 196–359 at the Fricke point — both overlap the Cabibbo window 91–896; the opposite
   cusp (49–90) touches it. Post hoc; flat-curve approximation.
3. **The curve is not flat on that scale (EXACT).** Gauss–Bonnet on genus 13: a canonical path of length l spans 12.28 l mean curvature radii at every A_X —
   1.3 (Fricke), 2.2–2.3 (2-torsion preimages), 2.6 (opposite cusp). The right tail object is the Green's function of the canonical (Arakelov) metric.
4. **Placement by arithmetic (EXACT, literature).** X0(143)(Q) = the four cusps (Faltings; Mazur–Kenku; cusps rational for squarefree N). P2's curve coordinate has been fitted.
5. **What survives (STRUCTURAL).** TO-01 (a ratio of tails), the lepton pattern (exact), b–tau normalisation and the Pati–Salam scale: four independent failures.

CC-55 supersedes CA-2. Next: the Arakelov Green's function of X0(143) between the blown-up cusp and the arithmetic homes, with a pre-registered test — r from V_cb must
land in the Cabibbo window 91–896. Closed forms for the cusp-to-cusp value (Abbes–Ullmo, N squarefree and prime to 6) give an independent check.
