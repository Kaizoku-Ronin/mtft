# MTFT prime traces v0.1.0

Start with `Prime_Trace_Guide.md`. The catalogue separates proved infinitude
from open infinitude, while retaining overlapping exact membership labels.

The default run covers the first 1,000 square intervals and all 78,650 primes
below 1,002,001. This is a standalone companion for the MTFT workflow; it
uses standard integer arithmetic and does not require importing MTFT.

## Reproduce

Python 3.12.13, numpy 2.3.5, and matplotlib 3.10.8 were used.

```sh
python -m pip install -r requirements.txt
python build_prime_traces.py --n-max 1000
python plot_prime_traces.py
python prime_trace_cli.py --prime 1093
python prime_trace_cli.py --interval 33
python prime_trace_cli.py --catalogue
```

The builder uses an in-memory sieve through 2(n_max+1)²+2 so that twin,
Sophie Germain, and safe-prime tests include the necessary partner values.
The supported n_max range is 60–3000. Memory and work grow with n_max²;
this first tool is not a distributed or segmented record-search engine.

## Data files

| File | Content |
|---|---|
| `family_catalogue.json` | Definitions, fixed bases/moduli, infinitude status, mathematical sources, necessary residue restrictions |
| `prime_traces.npz` | Exact integer and boolean arrays, described below |
| `trace_summary.json` | Counts, example primes, validation gates, computation bounds, runtime |
| `interface_checks.json` | Four completed CLI checks |
| `Prime_Trace_Web.png`, `.svg` | Two windows showing family memberships on the infinite residue classes |
| `Prime_Trace_Growth.png`, `.svg` | Observed cumulative counts and proportions through the full finite range |

`prime_traces.npz` stores arrays indexed either by prime or by square interval:

- Prime arrays: `p`, `square_n`, `position_numerator`, `position_denominator`,
  `residue_mod30`, `properties`, `fermat_quotient_base2`, `previous_prime`,
  `previous_same_family_gap`.
- Interval arrays: `n`, `upper_square`, `interval_total`, `interval_residue`,
  `interval_types`, `interval_exceptions`, `first_prime`, `cumulative_total`,
  `cumulative_residue`, `cumulative_types`.

The residue-column order is `[1,7,11,13,17,19,23,29]`. The property-column
order is `[twin_member,sophie_germain,safe,wieferich_base2]`. These constants
also appear in `prime_families.py` and `trace_summary.json`.

The exact position within a square interval is
`position_numerator / position_denominator`. The sentinel −1 means the
Fermat quotient is inapplicable at p=2, there is no previous prime, or a
same-family predecessor is absent/inapplicable. Boolean membership arrays
distinguish the latter cases. The CLI renders applicable missing values as
JSON null.

Residue rows partition the primes above 5. Property columns overlap; their
counts must not be summed as a partition. Infinitude of an intersection or
complement is not inherited from a parent family. No sampling, random seed,
floating-point primality decision, or automatic theorem-status update is used.
