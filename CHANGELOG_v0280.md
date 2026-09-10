# mtft v0.28.0 — CC-23 sampler fix, marked Hecke geometry, oldsector hardening (2026-09-10)

## CC-23 — `lattice.metropolis_sweep` used the wrong staple contraction
The plaquettes containing U_mu(x) sum to Re Tr(U_mu(x)·staple); the sampler contracted with
staple^dagger, so dS was wrong (sign wrong in 18/24 random proposals, max error 4.2 at beta=2).
Fixed; `tests/test_cc23_staple.py` compares local dS with full recomputation (3e-13).  Found by
Astra (QCD-05), verified here.  With CC-21 (hot starts not in SU(N)), every sampler-derived
lattice number before 0.28.0 must be rerun; cold-start plaquette values from `wilson_action`
itself are unaffected.  Note also: the sweep accepts on the Wilson term only — arithmetic
effects require full-action acceptance or validated reweighting (open).

## `surface.marked` (EXACT) — what a spectrum cannot see
PLANCK-02's witness: an additive Hamiltonian isospectral to the connected arithmetic interaction
({−108,−12,12,108}), so no trace of a function of D distinguishes them; and the exact heat action
Z(x) = Tr e^{-x D²} is strictly decreasing, so it selects no scale (the quartic truncation's
minimum is an artefact).  Marked observables do distinguish them: eta = 1 vs 0.  Module carries
the decomposition (tau-norm), degeneracy maps with their gates (J♯J = (r+1)I, J0♯J1 = a_r I),
cross-map transfers (diagonal blind to alpha; K01 − K10 = 320 alpha R reads it), and the two
competing objectives side by side: least leakage 409600/81 + 7200 alpha² (min at 0) versus
Herm(Hecke) distance 57600/41 (alpha−1)² (zero at 1).  The choice of principle is explicit.

## `surface.oldsector` hardening (Astra's v0.27.4 audit)
Prime/distinctness/Hasse-bound validation; `skew_correction` (alias `interaction_term`);
`hermitian_hecke_selection` reports status unique / not_unique / not_identifiable /
no_member_in_family; `good_prime_replacement_control` replaces a local operator by the scalar
a_ell (T_ell on the fixed oldclass) and shows skew correction AND connected term vanish — the
correct representation-specific "bad-prime only" statement; the earlier test only checked skewness.

## Recorded study verdicts (tags)
Herm(ordinary Hecke) selects alpha = 1: EXACT within the family; 28/28 prime pairs (Astra);
conditional on the admissibility principle.  Spectral heat action selects no scale: EXACT.
Color transport (QCD-06): compact su(3) on the forgetful-complement triplet after trace
removal — conditional; a non-arithmetic control gives su(3) too (generic mechanism).
Dirac pairing: still fails for every balanced template.  No physical scale, dimension, or
particle identification is claimed anywhere in this release.
Base: live 0.27.4 tarball.  Tests: 5 new.  Release pin four-way.
