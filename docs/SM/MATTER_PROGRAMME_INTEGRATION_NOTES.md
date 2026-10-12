# The matter programme — notes for the integration pass (2026-10-12)

Candidate: mtft 0.34.0 with the reconciled Matter-01 to Matter-03 additions, plus the matter programme (Matter-04 to PG-01, 2 to 12 October 2026)
and its consolidation into a programme map. Version strings stay at 0.34.0 until the integration pass; pyproject.toml, `mtft.__version__` and
CITATION.cff agree.

## What the candidate adds

- `mtft.matter`: 57 modules carrying ledger claims, plus `mtft.divisor_gas` on the arithmetic side (Legend entries), the ledger `mtft.matter.ledger` and the programme map `mtft.matter.puzzle`.
- Command line: `python -m mtft matter presets | background NAME | ledger | puzzle | levels | timerev | prereg | forces | atlas | epsilon | scan | ratio3 | leptons | phases`
  (also `python -m mtft.matter ...`). `puzzle` takes `--piece I...X`, `--open`, `--json`, `--markdown` and `--dot`; `levels` takes `--level N` and `--json`;
  `timerev` takes `--json` and `--candidate KIND VALUE` (the sealed I.4 test); `prereg` takes `--json` and `--verify` (exit 1 if any seal changed);
  `forces` takes `--json`, `--touch PARENT`, `--column PARENT` and `--box N`; `atlas` takes `--json`; `python -m mtft.divisor_gas` prints the arithmetic clock.
- Reports: one per wave in docs/SM; `MATTER_LEDGER.md` and `MATTER_PUZZLE_MAP.md` are generated and checked against their sources by tests.
- Provenance: `scripts/arakelov_kernel_x0143.py` reproduces AG-01's kernel, `scripts/star_quotient_scan.py` SQ-01's route-1 tables and
  `scripts/time_reversal_curves.py` TR-01's curve table (all need cypari2; the last also elldata); the modules freeze their numbers. One SQ-01
  test re-runs PARI at 15 levels and one TR-01 test at 3 curves; both skip without cypari2. SQ-01's list of the Monster's
  element orders comes from GAP 4.12.1 with CTblLib 1.3.7; no test needs GAP. GF-01's figures are produced by `scripts/gf01/` and PG-01's by `scripts/pg01/` (matplotlib only).

## The record

- Ledger: 278 claims, of which 220 active, 33 caveat, 19 superseded and 6 retracted. Corrections CC-36 to CC-65 (30) are append-only, and every
  retraction cites its correction (`ledger.check()`).
- Programme map: ten pieces, numbered in dependency order. I (parent action) open; II (arena) and III (twist and families) stand; IV (size and home),
  V (stack and breaking) and VII (Higgs, vertex and top) partial; VI (couplings and b–tau) and VIII (neutrinos) open; IX (quark flavour) and
  X (lepton flavour) closed with documented negatives. 73 sub-puzzles, 29 open; 22 recorded negatives. Every ledger claim belongs to exactly one piece (`puzzle.check()`).
- Pre-registered tests: T0 mixed, T1 no hit, T2 fail, T3 pass (modest evidence), T4 fail, T5 fail, T7 pass (charge quantisation, PR-S5; PR-S4 superseded by CC-64), T8 mixed
  (the touch rule on the parent menus, PR-S3), T6 pending (the strong-CP rule for I.4,
  the first hash-sealed test: `mtft.matter.prereg` PR-S2; PR-S1 was superseded the same day by CC-63 after independent review, and stays in the record). SQ-8 registers the test for IV.6 (why the level is 143); it has not been run, because
  no selection rule exists yet; PR-01's trials counter logs the seventeen properties tried after the fact (four single 143 out; X17 is GF-01's twelve-loop systole) and its seals are
  where any future rule must be registered before it runs.
- This consolidation's corrections: CC-61 retracts LEP-4 (the equal-overlap fits are not unique); CC-62 relabels the parity audit's wave PAR-01.
- SQ-01 (after the consolidation): X0(143) is one of 13 supersingular semiprime levels whose star quotient is a rank-1 elliptic curve (N20); what
  singles it out is arithmetic only. Its census is EXACT by two routes at 99 of 105 levels (route 2 alone at the other 6).
- TR-01 (2026-10-10): the thirteen quotients are identified as curves (the one of modular degree 4 in each class; 65a2 rather than the optimal
  65a1), nine with one real circle and four with two; time reversal on a gauge coupling has two mirror classes under duality (theta = 0 and
  theta = pi) that never mix, 143a1's lattice sits in the theta = pi class, and the root number is the t -> -t parity of Lambda.
- GF-01 (2026-10-11): the forces read as eight dots; two new seals with outcomes (PR-S3 mixed, PR-S4 pass); hypercharge by two routes; the
  couplings as the shape of the force torus; with PV-02's twist the U(8) column recovers PV-02 while SO(16) and E8 through U(8) are closed by
  the sum line's parity and the lepton rows (N21), which makes piece IX's down-sector wall parent-level within the supersymmetric menus; the clock
  of twelve and the systole of X0(143) (the twelve loop, not the golden one). Statements about SO(16) and E8 are conditional on the supersymmetric
  twist and the Standard Model embedded through U(8). CC-64 (same day) superseded four of the wave's claims and one seal after independent review.
- PG-01 (2026-10-12): the weights as the divisor-gas energy at the Hagedorn wall and the Glaisher identities (Legend, EXACT); the critical-point
  atlas of the vacuum (shape: a cosh well with a sink; volume: no critical point, given PV-01's positive coefficients; condensate: a flat direction
  lifted by a mass term with 0 < m^2 < xi; Wilson line: a sink and a source, exact in PA-01's product approximation); Hubble friction, the energy
  argument for capture and a DIAGNOSTIC toy run; the simulation gate recorded as I.7. Nothing maps the arithmetic clock beta to a physical
  temperature. CC-65 (same day) corrected the wave's RH sentence (the marked energy keeps the poles at the zeros; RH-blindness belongs to Z_D and
  the summatory function) and superseded PG-2 and PG-3 by PG-4 and PG-5 after independent review, which confirmed every number.

## For the auditor

- Status classes: the matter ledger uses EXACT, NUMERICAL, ESTIMATE, STRUCTURAL and MODEL; nothing in it is tagged CERTIFIED. NUMERICAL claims rest on
  one route unless the statement says otherwise. The map records the canonical kernel's second route (the Abbes–Ullmo comparison) as open (II.1).
- Environment: the full routine suite passes on Python 3.13, NumPy 2.5 and SciPy 1.18 (1333 passed after CC-65, 40 slow deselected). PARI/GP 2.15.4 arrived
  with the GAP install used for SQ-01, so the 4 PARI/GP tests that skipped before (1279 passed, 4 skipped) now run, and pass.
  Before this consolidation it passed on Python 3.12 with the earlier stack. The change of stack exposed one fit whose result depended on which exact
  solution SciPy reached first (CC-61); the fit now returns a stated canonical representative and the full solution set.
- Not run here: the 40 slow cases.
- SQ-01 had an independent check: a separate pass recomputed every census claim three ways (PARI modular forms, gp modular symbols on the full
  space, Riemann–Hurwitz with its own fixed-point counts) and found no discrepancy; it also found Mordell–Weil rank 1 (ellrank) at all 13 tori.
- Verification:
  - three-way version check (pyproject.toml, `src/mtft/__init__.py`, CITATION.cff);
  - `python -m pytest tests -m "not slow"`;
  - `python -c "from mtft.matter import ledger, puzzle; print(ledger.check(), puzzle.check())"` (every value True);
  - `python -m mtft matter puzzle --markdown | diff - docs/SM/MATTER_PUZZLE_MAP.md` (no output);
  - with cypari2: `python scripts/star_quotient_scan.py --max-level 400 --check` (prints 'differences from the frozen tables: none').
- Teaching material is kept out of the package: the plain-language board built from `puzzle --json` lives outside the release.
