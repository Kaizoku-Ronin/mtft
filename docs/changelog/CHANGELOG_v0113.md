# mtft v0.11.3 — deep structure wave (2026-08-05)

Three closures, zero new conjectures. Wave code: Claude Opus 5, run and
certified in-session; independently audited before push (Addendum BS,
Kimi K3 — every headline number re-derived on the auditor's own
machinery; two integration gaps completed by the auditor and disclosed
below).

## Closed — capture ceiling (Add. BQ leg 6 / AG-D5)

`studies/x0143_capture_parity.py`. The exact-1/2 ceiling is a theorem:
complex conjugation of X0(143) acts on the tessellation as
sigma:(c:d)->(-c:c+d), an involutive graph automorphism with 4 fixed
triangles (the real locus) fixing the width-1 nucleus. Odd states vanish
at fixed points => the odd sector cannot see the well; sigma-symmetric
dipoles obey <even|D|odd> = 0; the odd bottom (= the free spectral gap
0.2726) is a dark sink. Ceiling = ||P_even psi0||^2: exactly 1/2 from
the 52 moved triangles, exactly 1 from the 4 fixed ones; V0- and
gamma0-independent; 1% sigma-breaking of the dipole sends it to 1.
Certified to 2e-14 across all 56 starts.

Auditor confirmations (Addendum BS): sigma verified on an independent
P^1(Z/143) chain — involution on all 168 cosets, exactly 4 fixed cosets
[(0,1),(1,71),(11,1),(13,10)], one per fixed triangle; sigma =
conjugation . T^-1 exactly. Dark-sink identity verified to 2.4e-15
(odd-sector ground = free Fiedler mode, well amplitude 2.2e-17).
V0 sweep {2,4,8} and gamma0 sweep {0.01,0.05,0.20} leave ceilings fixed
to <= 9e-15; a 1% dipole asymmetry at one moved node sends every
ceiling to 1.000000 — the symmetry is load-bearing.

## Closed — f3 norm<->embedding pairing (the orbit-Zeno gate)

`studies/x0143_f3_pairing.py`. Measured, not assumed: Rankin slopes per
embedding, f1 point-count calibration, f2 known-pairing control. Shipped
sigma-order CONFIRMED as the unique optimal assignment; f3 pairing is
the minimizer over all 720 assignments with a 5.6x cost margin; near
pair 0.01008/0.01085 resolved (3.3 sigma in the session analysis;
auditor's independent estimates give >= 3 sigma under every noise model
tried, up to 25 sigma under window jitter). PET_F3 order:
DIAGNOSTIC -> Cert. Var(tau) = 4.1836 (point value; interval
[2.330, 4.252] retired — verified by the auditor as exactly the min/max
over all 720 pairings), Var(mu) = 2.2315. H4 verdict: ordering robust;
null-test caveat and lifetime non-correspondence unchanged.

Integration note: the wave shipped the a_p extraction only; the analysis
block (sieve + residue fits, calibration, 720-enumeration, Zeno point
values) was completed by the auditor on auditor-verified machinery and
reproduces every session number from the study's own extraction
(f2 control RMS 1.67%, f3 RMS 1.90%). The extraction code is unchanged.

## Added — the Hodge endomorphism tower

`tests/test_hodge_tower.py`. Fast tier: Hecke commutant dim 44,
dim_Q Q[T2] = 11 (= 1+4+6 = K1 x K2 x K3 = End(J0(143)_new) tensor Q),
Hecke fields totally real (Albert type I — no Weil-class pathology is
possible). Slow tier: J from the engine's own period functionals,
J^2 = -I, [T2,J] = 0, Hodge-cut commutant dim 22. Places the generation
decomposition on the proven (divisor-generated) side of the Hodge
Conjecture. Auditor's measured residuals (Addendum BS): ||J^2+I|| =
2.2e-16, ||[T2,J]||/||T2|| = 5.5e-15, nullities 44/22/11 confirmed.

Version: 0.11.3 three-way (pyproject / __init__.py / CITATION.cff).
Open items carried: Paper 32 tau erratum (paper-side); React drawn-loop
stage; weight >= 4 / Kuga-Sato frontier.
Audit lineage: BQ (Kimi) -> v0.11.1 (Claude) -> BR (Kimi) -> v0.11.2
(ChatGPT catch, Kimi fix) -> v0.11.3 wave (Claude) -> BS (Kimi).
