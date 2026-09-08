# MTFT irregular primes and Legendre intervals

This standalone companion experiment adds exact regularity labels and all
Bernoulli witnesses to the finite prime-trace work for MTFT 0.26.0. The default
census covers every prime below 10,201 and the first 100 square intervals.
It preserves both the proved-infinite irregular family and the regular family,
whose infinitude is open. The prime 2 is explicitly unclassified.

Read **Irregular_Legendre_Report.md** for the findings and limitations.
The two PNG/SVG figures show the trace web and conditional null comparisons.

## Reproduce

Use Python 3.12 and install `requirements.txt`. From this directory:

```bash
python build_census.py
python analyze_traces.py
python validate_experiment.py
python render_results.py
```

The original run reconciled the overlap with the prior study using:

```bash
python build_census.py --prior-census ../MTFT_Prime_Traces_v0.1.0/prime_traces.npz
```

The previous census is optional; it is not needed for the computation or
validation. The modular arithmetic and analysis need NumPy. SymPy is used only
for independent rational validation, and Matplotlib for figures.

Default randomizations: 49,999 per model/range, with explicit PCG64 seeds.
`analyze_traces.py --repetitions 999` is a faster exploratory run; it changes the
Monte Carlo resolution and is not the archived experiment.
`build_census.py --n-max` accepts 32–200 for separate arithmetic experiments,
but this report and figure layout describe the fixed default n=1–100 design;
use the default for reproducing the released report.
The report renderer checks for that default design. Seeded statistics and exact
data values reproduce with the pinned versions; timing fields and image/archive
metadata can change their byte hashes on a fresh run.

## Inspect records

```bash
python inspect_trace.py --prime 1093
python inspect_trace.py --prime 3511
python inspect_trace.py --interval 9
```

All interval endpoints are strict. For example, interval 9 is (81,100).
The existing twin, Sophie Germain, safe and base-2 Wieferich marks are also
included and can overlap; they are not mutually exclusive prime types.

## Data schema

- `census.npz`: prime values, square indices, exact fractional positions,
  residue modulo 30, regularity labels, irregularity index, complete even
  Bernoulli residues, per-interval counts and the four earlier prime marks.
  `bernoulli_offsets[j]:bernoulli_offsets[j+1]` locates prime j's residues in
  the flattened `bernoulli_residues` array. Entry r corresponds to k=2+2r.
  The residue array is uint16; this supported range stays below 65,536.
  Prime 2 has irregularity index -1 and an empty residue slice. Prime 3 is
  regular with index 0 and also an empty slice, since its criterion is vacuous.
- `prime_witnesses.json`: every prime's full list of zero-residue indices.
- `family_catalogue.json`: domain, membership criteria and infinitude status for
  the two added families. Infinitude status is never assigned to one prime.
- `census_summary.json`: exact counts and overlap reconciliation.
- `analysis_results.json`: all statistics, exact hypergeometric probabilities,
  strata, pointwise interval envelopes, raw p-values and Holm adjustments.
- `permutation_statistics.npz`: the actual randomization statistics. Columns:
  irregular-free interval count, longest run, thin-minus-other irregular
  fraction. Keys identify the selected n range and null model.
- `validation_results.json`: finite arithmetic and null-model verification.
- `EXPERIMENT_PLAN.md`: scope and tests fixed before new labels were computed.
- `PROVENANCE.json`, `SHA256SUMS.txt`: source and artifact integrity records.

## What the certificate means

The full inverse identity in F_p[t]/(t^(p−2)) validates every even Bernoulli
residue for each odd prime in the census. All positive witnesses are also
checked by a separate modular-power-sum calculation. These are exact finite
checks, not a Lean kernel proof or an infinitude proof. The permutation models
condition on prime positions and therefore cannot test Legendre itself.
