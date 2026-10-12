# PR-01 — sealed pre-registrations and the trials counter (2026-10-10; corrected the same day, CC-63; PR-S3, PR-S4 and X17 added by GF-01, 2026-10-11)

Module `mtft.matter.prereg`; tests `tests/test_pr01_prereg.py`; command line `python -m mtft matter prereg [--json] [--verify]`.

1. **Seals (STRUCTURAL).** A pre-registered rule is sealed by the SHA-256 of its id, its statement, the source code of the callable that evaluates
   it, its parameters and its pass window, dated before any candidate is evaluated. `verify()` recomputes the hash from the code as it stands;
   `evaluate()` refuses a rule whose seal no longer verifies. A rule can therefore only be changed by re-sealing it under a new id with the old
   entry kept — append-only, as the CC corrections (PR-1).
2. **Legacy (STRUCTURAL).** T0–T5 were stated in prose before evaluation (puzzle map) but never hash-sealed; they stay listed as legacy and are not
   re-sealed after the fact. PR-S1, the strong-CP rule of TR-01 (puzzle test T6), was the first sealed entry; on the day of sealing the
   independent review found that it rejected theta = pi mirrors outright although theta-bar = theta + arg det(M_u M_d) lets a theta = pi point
   with arg det = pi pass. The registry handled it as designed: PR-S1 stays in the record as superseded, with its hash, `evaluate()` refuses it,
   and PR-S2 carries the corrected rule under a new hash (CC-63). Its outcome is pending because no candidate exists (PR-2).
   GF-01 (2026-10-11) sealed two more and ran them the same day: PR-S3, DS-01's touch rule on the three supersymmetric parent block menus
   (outcome mixed: U(8) FAIL, SO(16) and E8 PASS — in columns that do not give the exact Standard Model, GF-11), and PR-S4, charge quantisation
   through the Z6 of the Standard Model (PASS on every line; passes for every odd T_3R coefficient, so necessary and not what fixes hypercharge).
   The independent review found PR-S4's statement read the global form backwards (Z6 trivial on the matter allows G/Z6 or any of its covers,
   not a further quotient); the criterion was unchanged, so PR-S4 stays as superseded with its hash and PR-S5 carries the corrected statement
   (CC-64), the registry's second use as designed. Hashes in `prereg.SEALS`.
3. **The trials counter (EXACT).** Every property tried as a reason for the level 143 is logged against the arena of IV.6, the thirteen torus
   levels, with the subset it selects:

   | id | source | property | selects | 143 |
   |---|---|---|---|---|
   | X01 | SQ-3 | twin primes | 143 | alone |
   | X02 | SQ-3 | consecutive supersingular primes | 77, 143 | included |
   | X03 | SQ-3 | new dimension = smaller prime | 57, 58, 65, 91, 143 | included |
   | X04 | SQ-3 | new dimension = p and consecutive primes | 143 | alone |
   | X05 | SQ-3 | new dimension = larger prime | none | no |
   | X06 | SQ-3 | orbit pattern (1, 4, 6) [self-referential: 8 of 13 levels have a unique pattern] | 143 | alone |
   | X07 | SQ-4 | index 168 [self-referential: 5 of 13 unique] | 123, 143 | included |
   | X08 | SQ-4 | genus 13 [self-referential: 2 of 13 unique] | 123, 143, 145 | included |
   | X09 | SQ-4 | torsion-free | 7 levels | included |
   | X10 | SQ-5 | Monster element order | 57 | no |
   | X11 | Holes and Charges note | N = 7 mod 8 (Legendre's exception; h(-4N) = h(-N)) | 143 | alone |
   | X12 | Holes and Charges note | quotient supersingular at 2 (a_2 even) | 57, 77, 91, 123, 141, 143, 155 | included |
   | X13 | TR-3 | quotient has one real circle | 9 levels | included |
   | X14 | TR-3 | quotient has no rational torsion | 10 levels | included |
   | X15 | TR-3 | optimal curve has modular degree 4 | 12 levels | included |
   | X16 | TR-3 | odd level | 9 levels | included |
   | X17 | GF-01 | the shortest closed loop of X0(N) is the twelve loop (trace 4, discriminant 12) | 118, 141, 142, 143 | included |

   Seventeen tried, none pre-registered; four single 143 out, fifteen include it (chance 5.77 for rules chosen blind). X17 was computed in the
   GF-01 reply before it was logged, like the rest. Every property was chosen by looking at 143, so "includes 143"
   is near-certain by construction, and the chance columns of the report (5.46 inclusions and 0.38 single-level hits expected for rules chosen
   blind to the level) apply only to sealed rules, of which none has been run. X12 was first logged with a_2 = 0; the independent review corrected
   it to a_2 even, the supersingular condition at 2, which adds 57a1 and 91a1 (a_2 = -2) and moves the baseline from 5.31 to 5.46 (CC-63). The record is a look-elsewhere account, not evidence (PR-3;
   AG-D5, SQ-8). The list is append-only, and the test suite refuses a distinction in `star_quotients.distinctions()` that is not logged.

Reading: the registry turns SQ-8 from a sentence into a mechanism. A rule that selects the level now has to be sealed here, with its code, before
it is run on the thirteen levels; until one is, IV.6 stays open whatever the counter says.
