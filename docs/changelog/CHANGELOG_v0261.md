# mtft v0.26.1 — Ising route C and the three-route sum rule (2026-09-05)

`surface.ising.density_of_states(N, negative_edges=())`: exact D_k by variable
elimination with base-2^(n+1) packed polynomials (independent implementation
of the method Astra applied to v0.26.0; width 11 and 0.1 s at N = 143).
`even_subgraph_polynomial` (cut/cycle transform) and `sum_rule_gate(N)`:
the signed sum of all 4^g Kasteleyn Pfaffians must equal A(t).

Results: routes agree to <= 1e-14 at N = 6, 11, 15, 35, 55, 77 (genus 0-7;
brute force joins through genus 5).  At N = 143 the elimination reproduces
Astra's 85 integers D_k and A_j exactly; A(1/sqrt3) = 14.751840885840 is now
the exact target for `full_sum_job` (the 4^13 Pfaffian sum).  Cross-engine
E2 at genus 7: this package's Pfaffian sum vs Astra's polynomial, 8.7e-15.

Recorded negative control (Astra): N = 6 and N = 11 have identical dual
graphs, so no scalar zero-field thermal observable can detect genus; genus
enters only through the spin-structure terms, which the Pfaffian route
resolves and the elimination route sums.

## Release-pin fix (2026-09-06, from Kimi's audit)
The version gates in `test_origami_v0200.py`, `test_periods_v0210.py`,
`test_periods_v0220.py` and `test_v0240_kakeya.py` had shipped stale for four
releases.  They now read a single constant `tests/_release_pin.py::RELEASE_VERSION`,
and `tests/test_release_pin.py` checks it against pyproject.toml, `__init__`
and CITATION.cff.  A bump is now four files, one of them a test.
