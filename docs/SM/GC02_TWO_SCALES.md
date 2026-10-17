# GC-02 — two scales, not one: the couplings matched at the curve's compactification scale (2026-10-17)

Module `mtft.matter.two_scales` (`ps_b`, `ps_inverse_couplings`, `two_scale_ratio`, `compactification_scale`, `two_scale_report`);
`tests/test_gc02_two_scales.py`. Follows VE-5.

## The point

GC-01 matched the flux-split couplings — R = (alpha_R^-1 - alpha_4^-1)/(alpha_L^-1 - alpha_4^-1) equal to a quantised design ratio — at the
Pati–Salam breaking scale M_PS, running the Standard Model one loop from m_Z up to it. That identifies two scales: where Pati–Salam breaks,
and where the eight-dimensional theory is matched to four dimensions. VE-5 separates them. The curve's Kaluza–Klein tower starts at
M_C = 1/(2 R_C) = 0.063–0.070 M_KK, and the breaking is a four-dimensional condensate below the KK scales (SC-01's reading R1), so M_PS <= M_C.
The consistent one-loop chain is: Standard Model to M_PS; alpha_4 = alpha_3, alpha_L = alpha_2, alpha_Y^-1 = alpha_R^-1 + (2/3) alpha_4^-1 there;
Pati–Salam running to M_C; the tree-level flux splitting as the boundary condition at M_C. GC-01 is the special case M_C = M_PS.

## What Pati–Salam running does (GC-6, EXACT at one loop)

With PS-01's content, 3 (4,2,1) + 3 (4bar,1,2) + 6 (1,2,2) Weyl fermions and n_H Higgs bidoublets, b_4 = -32/3 and b_L = b_R = 2/3 + n_H/3:
left–right symmetric for any bidoublet Higgs content (only a (4bar,1,2) scalar would split them). Both differences alpha_{L,R}^-1 - alpha_4^-1
therefore shrink at the same rate above M_PS, and since alpha_R^-1 - alpha_4^-1 is the larger, R rises. Consequences:

- the designs with R below the Standard-Model value at M_PS — 1/9, 3/7, 1, 7/3 and 1/5 — are reached only with M_C = M_PS: GC-01's scales
  8.9e11–2.1e14 GeV stand for them;
- the designs GC-01 listed as unreachable, R = 5 (single point) and R = 9 (two point), are reached after Pati–Salam running, at a
  computable M_C(M_PS).

## The numbers (GC-7, NUMERICAL)

| M_PS (GeV) | R at M_PS | M_C for R = 5 | M_C for R = 9 |
|---|---|---|---|
| 1.12e6 (the b–tau crossing) | 3.678 | 1.7e7 | 2.8e8 |
| 1e7 | 3.57 | 1.5e8 | 1.7e9 |
| 1e9 | 3.25 | 1.4e10 | 7.7e10 |
| 8.9e11 (GC-01's lowest) | 2.33 | 1.1e13 | 2.1e13 |

With M_PS at the b–tau crossing: M_KK = M_C/(0.063–0.070) = 2.5e8–2.8e8 GeV (R = 5) or 3.9e9–4.4e9 GeV (R = 9), M_SS = 0.4943 M_KK =
1.2e8–1.4e8 or 2.0e9–2.2e9 GeV. y_b = y_tau at M_PS by construction, so GC-4's 30–35 per cent correction is not needed on these designs and
CC-48's caveat on 1.1e6 GeV becomes conditional on GC-01's single-scale matching. Content sensitivity: without the six bidoublet fermions
M_C = 7.2e7 and 4.9e9 GeV; a second Higgs bidoublet or a (4bar,1,2) scalar moves it by under 25 per cent.

Consequences elsewhere: SC-2's M_KK window (8.9e11–2.1e14 GeV) assumed the single-scale matching; at M_KK = 2.5e8–4.4e9 GeV, g_4^2 = 0.47
and lambda = 109–133, so every conclusion of SC-2 and SC-3 holds a fortiori (SC-2 carries the note). VE-6's GeV^4 ledger scales with M_KK^4 and
its conclusions do not change.

## Not included (as in GC-01)

The one-loop thresholds at M_C — the curve's own Kaluza–Klein towers, VI.2, now with a definite scale to be evaluated at; which background
realises R = 5 or 9 (the canonical matter04 background has R = 1); and the Pati–Salam leptoquark bound from K_L -> mu e, which puts M_PS near
or above 1e6 GeV (ESTIMATE), so the b–tau crossing sits at the edge of what is allowed.

## Package changes

New `mtft.matter.two_scales` (five functions on top of GC-01's running); ledger GC-6, GC-7 (wave GC-02; 296 claims); GC-4, CC-48 and SC-2 carry GC-02's
caveats; puzzle VI.4 reopened and VI.3 updated; wave registered under Piece VI; `tests/test_gc02_two_scales.py` (4 tests).
