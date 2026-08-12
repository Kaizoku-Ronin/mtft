# CC-06 (append-only) — M8 study mixed-basis bug

**Filed:** 2026-08-12, v0.15.0 integration, as a **proposal** (auditor correction of
a wave artifact; the wave's own files are byte-preserved). **Class:** [Cert] —
reproduced verbatim by the auditor's independent scripts (`m8_verify_kimi.py`,
`m8_deep_kimi.py`, JSONs alongside).

## The bug

`m8_hecke_commutant_study.py` mixes bases when intersecting the Hecke commutant
with the V-operator commutant, and reports dim = 1.

## Corrected results (six primes each, auditor reproduction)

    commutant(Hecke)         = 60   (study agrees)
    commutant(V)             = 82   (study agrees)
    full Hecke algebra dim   = 52
    commutant({Hecke, V})    =  2   (study printed 1)

The joint commutant is span{I, Z} with Z/2 = (I - iota*)/2 the canonical
involution projector.

## Downstream correction

The canonical M8b amplitude is **0.20610964892935077**, not the 0.3805 quoted in
the m8b ledger (which used the study's mixed basis). The m8b study files are
byte-preserved; the corrected value is the one wired into the certified
CURVATURE-facing tables.

## Evidence

- `m8_verify_kimi.py` / `m8_verify_kimi.json`: six-prime commutant dimensions.
- `m8_deep_kimi.py` / `m8_deep_kimi.json`: basis-mixing diagnosis; Z projector
  identification.
