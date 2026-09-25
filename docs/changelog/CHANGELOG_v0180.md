# v0.18.0 — `mtft.boundary` and `mtft.gprun`

Two new subpackages. No behaviour change elsewhere.

## `mtft.boundary` — the RT-1B cusp/elliptic boundary arc

The canonical quadratic differentials sit strictly inside the weight-4
Hecke module, and the gap is an arithmetic object:

    0 -> S_4^(2)/H^0(2K) -> S_4/H^0(2K) -> S_4/S_4^(2) -> 0
         elliptic quotient      B_N        cusp quotient

with dim B_N = nu_infinity + nu_2 + nu_3. Good primes saturate both
quotients; level primes degenerate along explicitly identifiable channels.

* `CENSUS` — the 22-level discovery table (N, p, q, genus, b_N, nu_2,
  nu_3, good signature, U_p, U_q).
* `RELATIONS_143` — the exact cusp-leakage description at level 143.
* `THEOREM_TARGETS` — the two separated open targets.
* `HYPOTHESES` — Q1/Q3 hold, **Q2/Q4 falsified**, H1 pre-registered and
  not yet resolved.
* `_data/` — S_4(Gamma_0(143)) basis to q^140, the four cusp functionals
  c_d(F) = a_1(F | W_d), and T_2 / U_11 / U_13 as exact matrices, with
  `PROVENANCE.txt`. Under `src/`, following the `canonical` pattern.

### Gates (`mtft.boundary.gates`)

| gate | re-derives |
|---|---|
| `gate_census_consistency` | b_N = 4 + nu_2 + nu_3 from Kronecker symbols, all 22 levels |
| `gate_hypotheses` | Q1/Q3 true and **Q2/Q4 false** — a pass here means the table was not edited |
| `gate_h0_2k` | dim H^0(2K) = 36, dim S_4 = 40 |
| `gate_cusp_functionals` | all four c_d vanish identically on H^0(2K) |
| `gate_leakage` | ranks 4 / 2 / 1, vanishing support, and c_1 = c_11 |
| `gate_t2_from_qexpansions` | rebuilds T_2 from the shipped q-expansions; discrepancy **0** |

The last one means a corrupted operator matrix cannot pass. The falsified
hypotheses are gated as *false* deliberately: silently reviving a dead
hypothesis is the failure this arc spent four rounds learning to catch.

## `mtft.gprun` — a job runner for multi-hour GP computations

    py -m mtft.gprun

A local browser UI on 127.0.0.1: drop a `.gp` file, Run/Stop toggle, live
log tail, download `output.txt`. Standard library only; nothing leaves the
machine.

Built for the RT-1B.6 holdout, which is a ~4 hour job, but the real
motivation is provenance. Every run:

* SHA-256s the script **before** launch and writes the hash into the
  output header *and* footer, so a result can never drift from the bytes
  that produced it;
* freezes a copy of the script inside the run directory;
* records gp's version, start/end times, wall clock and exit code;
* flushes continuously, so an interrupted run leaves a readable partial
  log, explicitly stamped `PARTIAL run, not a result`.

Runs land in `mtft_runs/<timestamp>_<name>/` as `script.gp`,
`output.txt`, `meta.json`. Set `MTFT_GP` if gp is not on PATH.

**Bug found by its own test suite.** `stop()` tested `if self.proc` and so
did nothing when Stop was pressed before `Popen` returned — leaving a job
running with no way to reach it. Now a cancel flag is set and honoured as
soon as the process exists. Regression test:
`test_stop_before_process_exists`.

## Tests

`tests/test_boundary_and_gprun.py` — 15 tests (7 fast, 6 `slow`, 2 gprun
integration), including a provenance test that SHA-1s the five boundary
data files.

## Studies

`studies/rt1b_2026aug/` — the full arc: RT-0 symmetry ceiling, the
census, the anatomy, the cusp leakage matrices, both pre-registrations,
`FREEZE_HASHES.txt`, and the partial holdout log.

## Not done

RT-1B.6 holdout (frozen, 8 of ~46 levels, no verdict); the two theorem
targets; `mtft.periods`; the CC-09 downstream sweep.
