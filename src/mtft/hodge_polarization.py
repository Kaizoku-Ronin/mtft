#!/usr/bin/env python3
"""
mtft.hodge_polarization — pinned polarization data for H1(X0(143)).
===================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

Source artifact: session "Modular curves and encryption" (2026-08-07):
star involution eta from (c:d) -> (-c:d) on Manin symbols; eta^2 = I to
6.7e-16, [eta, T_h] = 0 to 7.8e-16, splitting 26 = 13 + 13.
J = eta o (scale by lambda_f per Hecke line); g(x,y) = omega(x, Jy).

STATUS of LAMBDA_TABLE: computed-not-certified (that session's own caveat):
single path pair (Omega^- from {0,oo}, Omega^+ from Re{oo,1/3}); ratios
well-defined only up to rational lattice normalization; second-path E2 pending.
f1 anchor EXACT from elliptic curve 143a1: tau = 0.5 - 1.0232745927 i,
analytic rank 1 (L(f1,1) = 0) — independently re-certified by mtft.gl2_peel.

CORPUS CORRECTION LEDGER (append-only; every id unique, never reused):
CC-01: Paper 26 sec.8 complex eigenvalue retracted (totally real fields).
CC-02: w-series shift chain: sum w_n n^-s = F(s+1) with F(s) = -zeta(s) zeta'(s+1);
       Paper 1 Prop 1.5's -zeta'(s+1) is the summand series, not the weight
       series (AG Pr 4.1.4 mistranscription). E2 three-route (W1 study sec.1).
CC-03: finite_atom_curvature fixed precision floor silently wrong at beta=80;
       adaptive doubling now, rtol 1e-20 (curvature.py, test_area_geometry.py).
CC-04: cold amplitude closed form c = (9 L5^2)/(25 L2 L3) [1 - (9/5) L5/L3
       + (4/5) L5/L2] = 0.27012646530542495706433719670365; retracts the
       six-atom beta=200 extraction 0.270126465305424759517602 (atom-6
       contamination (5/6)^200 = 1.46e-16 vs observed 1.98e-16).
CC-05: MTFT_Chapter5 Table 1 misstates w_6, w_8, w_9, w_10; corrected in
       studies/v0150_2026aug/CC-05_chapter5_weight_table.md.
CC-06: M8 study mixed-basis bug: true dim commutant({Hecke, V}) = 2 = span{I, Z}
       (study claimed 1); canonical M8b amplitude 0.20610964892935077.
CC-07: du03 "5-dim commutant" was edge-metric-specific; joint (graph-local AND
Hodge-natural) admissible space is 1-dim (scalars), metric-independent
across the 43-dim invariant cone (8 sampled metrics, cond 201-425).
       (Renumbered from CC-02 on 2026-08-12; see CC-07_du03_metric_correction.md.)
"""
ETA_CERT = {'eta_sq_minus_I': 6.7e-16, 'commutator_Th': 7.8e-16,
            'split': (13, 13)}
COND_G_REPRESENTATIVE = 112.3

LAMBDA_TABLE = {          # |Omega^-/Omega^+| per Hecke line
 ('f1', 1): 1.0232745926964612,   # EXACT (curve 143a1)
 ('f2', 1): 0.9804748832, ('f2', 2): 5.7895654701,
 ('f2', 3): 0.4555644526, ('f2', 4): 2.2811927490,
 ('f3', 1): 6.9514450547, ('f3', 2): 1.9783354598,
 ('f3', 3): 1.7560581532, ('f3', 4): 2.4059869500,
 ('f3', 5): 0.0792244877, ('f3', 6): 0.1517936036,
}
LAMBDA_STATUS = "computed-not-certified: single path pair; second-path E2 pending"

SATAKE_EXAMPLES = {2: (0, 2), 3: (-1, 3), 7: (-2, 7)}   # a_p, |alpha_p|^2 = p
