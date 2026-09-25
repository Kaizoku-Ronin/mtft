# mtft v0.31.4 — response to the v0.31.1 review; exact anomaly extension; TRI-02 D8 lift (2026-09-16)

Review V0311 (Astra/Sol) adopted: relative rank tolerance (`rrspace._rank_rel`) and rank-<=2 wording for the W13 texture;
`smflux.abelian_anomaly_polynomial` / `mixed_nonabelian_anomalies` (cubic and mixed anomalies vanish exactly on
span{phase, Y, B-L}) and `smflux.majorana_obstruction` (nu^c has B-L = +1: bare Majorana mass forbidden while B-L is
unbroken); registers ship in `docs/SM/` (sdist); SM-14 amended.  TRI-02 (EXACT): `arithspin.al_lift_on_S0` — the
Atkin–Lehner V4 lifts to D8 on the half-forms H^0(S0) (A^2 = +1, B^2 = -1, commutator -1), answering the Pin-lift
question of the Klein-four thread.  Tests: `test_surface_review_v0311.py`.  Base: the 0.31.3 candidate.  Pin four-way.
Note for the auditor: the published 0.31.1 changelog header is dated 2026-09-12, before the 0.31.0 consolidation.
