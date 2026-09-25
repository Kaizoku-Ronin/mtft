# mtft v0.30.0 — the Standard Model campaign: flux models, anomaly ledger, arithmetic spin (2026-09-13)

Consolidates 0.29.1 (KK05 corrections, W11 parity theorem), 0.29.2 (`surface.arithspin`) and adds
`surface.smflux` (SM-01): exhaustive multi-stack flux search on X0(143) with hypercharge embedding,
representation tables, mixed-anomaly ledger with the anomaly-free U(1) span, Higgs-lifting analysis, and
Brill–Noether purity of chiral sectors with the even arithmetic spin structure S0.
Exact results: no four-stack model has three chiral families (class theorem); 30 five-stack models, all with
m_L = −3 and >= 6 vector-like doublet pairs; representative M1 = degrees (0,−3,3,3,0), y = (1/6,0,−1/2,1/2,−1/2)
with 3 families + 3 nu^c + 6 Higgs pairs, hypercharge anomaly-free and in the anomaly-free span; realised on
X0(143) by W13-fixed CM divisors (W11 broken as required), with pure chiral sectors.
Open (declared): Green–Schwarz sector, the SUSY-convention search (indices shifted by −12), Yukawa cup products
for the CM-twisted bundles (step 5 tool), scale and couplings.  Tests: `test_surface_smflux.py` (2, ~10 s).
Base: live 0.29.0 (SHA-256 a5a0d5e0…) + unpushed 0.29.1/0.29.2.  Pin four-way.
