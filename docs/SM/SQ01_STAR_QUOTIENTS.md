# SQ-01 — the star quotients X0(pq)* for supersingular p < q: is X0(143) special? (2026-10-08)

Module `mtft.matter.star_quotients`; tests `tests/test_sq01_star_quotients.py`; provenance `scripts/star_quotient_scan.py` (needs cypari2);
command line `python -m mtft matter levels [--level N] [--json]`.

X0(N)* = X0(N)/<W_p, W_q> is how X0(143) folds onto 143a1. The question asked: can another level built from two supersingular primes do the same?

1. **Two routes (EXACT).** Route 1 (PARI/GP 2.17.2, frozen): the level-N newforms whose Atkin–Lehner signs are all +1. Old forms add nothing, because
   every newform of supersingular prime level has w_p = -1 (Ogg). Route 2 (pure Python, computed by the tests): each involution's fixed points from
   class numbers of imaginary quadratic orders with Eichler symbols, plus the Z[i] term for W_2, then the holomorphic Lefschetz formula
   tr(W | H^0(K)) = 1 - nu/2 and g* = (1/|W|) sum tr. Route 2 covers all 105 levels and agrees with route 1 at all 99 that route 1 reached
   (2201, 2419, 2773, 2911, 3337 and 4189 have route 2 only), and at all 75 levels with three or four supersingular primes up to 700.
   At 143 the fixed-point counts are 0 (W_11), 4 (W_13) and 20 (W_143).
2. **The census (EXACT).** X0(N)* is a sphere at 21 levels (6, 10, 14, 15, 21, 22, 26, 33, 34, 35, 38, 39, 46, 51, 55, 62, 69, 87, 94, 95, 119), an
   elliptic curve at 13 (57, 58, 65, 77, 82, 91, 118, 123, 141, 142, 143, 145, 155) and of genus >= 2 at 71. Every elliptic quotient has analytic
   rank 1 (NUMERICAL, PARI lfunorderzero); odd rank is forced, since all signs +1 give Fricke sign +1 and root number -1 (SQ-1, SQ-2).
   **Honest negative N20:** X0(143) is not the only supersingular semiprime level that folds onto a rank-1 elliptic curve.
3. **What singles out 143 (EXACT, found after the fact).**

   | Property | Levels among the 105 |
   |---|---|
   | Twin primes, star quotient elliptic | 143 only (twin ladder 15, 35, 143, 323, 899 -> genus 0, 0, 1, 5, 12) |
   | Consecutive supersingular primes, star quotient elliptic | 77, 143 |
   | dim S2^new = p, the smaller prime | 26, 38, 39, 51, 57, 58, 65, 69, 85, 91, 143 |
   | dim S2^new = p with p, q consecutive | 143 only |
   | dim S2^new = q, the larger prime | 221, 247, 299, 377, 403, 533, 611, 767, 923 (all p = 13); with consecutive primes, 221 only |
   | Galois-orbit pattern (1, 4, 6) | 143 only |

   Signs are constant on Galois orbits: at 143 the orbits of dimension 1, 4 and 6 carry (w_11, w_13) = (+,+), (-,+), (+,-); with the two old forms
   from 11a the sectors (+,+), (+,-), (-,+), (-,-) have dimensions 1, 6, 5, 1 (SQ-3).
4. **The nearest sibling (EXACT).** 123 = 3 x 41 is the only other level of index 168 and, with 143, the only torsion-free level of genus 13; it has
   four cusps, area 56 pi and a rank-1 elliptic quotient. Its new space has dimension 7 (orbits 1, 1, 2, 3) against 143's 11 (1, 4, 6) (SQ-4).
5. **Moonshine cross-check (EXACT).** Every level in the scan whose star quotient is a sphere — the 21 above and, with three primes, 30, 42, 66, 70,
   78, 105, 110 — is the order of an element of the Monster (GAP 4.12.1, CTblLib 1.3.7: 73 element orders over 194 classes). This is Ogg's
   observation for primes carried to products. 143 is not an element order; among the semiprime levels the converse fails at 57 and 93 (SQ-5).
6. **Denser factorizations (EXACT).** With three or four supersingular primes and N <= 700, X0(N)* is an elliptic curve at 102, 114, 130, 138,
   174, 182, 190, 195, 210, 231, 238. Each extra prime doubles the cusps (8 or 16) and the Atkin–Lehner group, so these change the arena rather
   than extend it (SQ-6).
7. **Flavour (STRUCTURAL).** The exact obstructions behind pieces IX and X do not use the level: DS-1's stack pattern, the odd-mode count (an odd
   B-stack mode costs 3(2g - 2) chiral families, 24 to 96 over the torus levels and 72 at 143, as OM-3 found), QL-1's dressing bound and CI-2's
   universality. DS-4's tail cap and BL-01's Wilson-line numbers were computed on X0(143) x 143a1 only (SQ-7). Another level would not reopen IX.
8. **Pre-registration for IV.6 (STRUCTURAL).** A rule that selects the level must be stated before it is evaluated on the 13 torus levels. It passes
   only if it selects 143 and rejects the other twelve, 123 and 77 included. The distinctions of item 3 were found after the fact and are not such
   a rule (SQ-8).

Reading: the arithmetic characterises 143 (twin primes; new dimension p with consecutive primes; the orbit pattern) but does not derive it. Why the
programme lives at 143 is now an open sub-puzzle of piece IV, with its test fixed in advance.
