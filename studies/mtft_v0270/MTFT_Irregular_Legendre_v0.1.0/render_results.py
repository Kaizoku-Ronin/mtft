"""Render standalone scientific figures and a report from saved experiment data."""
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
REGULAR = '#557b9d'
IRREGULAR = '#d98523'
GLOBAL = '#3c6588'
STRATIFIED = '#318879'
OBSERVED = '#a53343'
INK = '#263444'
RESIDUES = [1, 7, 11, 13, 17, 19, 23, 29]


def save_figure(fig, stem):
    fig.savefig(ROOT / f'{stem}.png', dpi=180, facecolor='white')
    fig.savefig(ROOT / f'{stem}.svg', facecolor='white')
    plt.close(fig)


def census_figure(data, summary):
    fig, axes = plt.subplots(2, 1, figsize=(12.5, 8.6), gridspec_kw={'height_ratios': [1, 1.18]})
    fig.subplots_adjust(left=.072, right=.978, top=.87, bottom=.17, hspace=.38)
    fig.suptitle('Regularity across the first 100 square intervals', x=.072, y=.975, ha='left', fontsize=19, weight='bold', color=INK)
    fig.text(.072, .93, f"{summary['irregular_count']:,} irregular primes · {summary['regular_count']:,} regular primes · prime 2 kept as an exception", fontsize=11.5, color=INK)
    n = data['n']; ax = axes[0]
    ax.axvspan(.5, 31.5, color='#eff1f3', zorder=0)
    ax.bar(n, data['interval_regular'], color=REGULAR, width=.84, label='Regular')
    ax.bar(n, data['interval_irregular'], bottom=data['interval_regular'], color=IRREGULAR, width=.84, label='Irregular')
    ax.bar(n, data['interval_exceptions'], bottom=data['interval_regular'] + data['interval_irregular'], color='#899099', width=.84, label='Prime 2')
    ax.axvline(31.5, color='#7a8791', linestyle='--', linewidth=1)
    ax.set(xlim=(.3, 100.8), ylabel='Primes in the interval', xlabel=r'Square interval index $n$: $n^2 < p < (n+1)^2$')
    ax.set_ylim(0, float(data['interval_total'].max()) * 1.21)
    ax.text(4, ax.get_ylim()[1] * .92, 'Earlier inspected: n ≤ 31', fontsize=9.5, color='#52606d')
    ax.text(35, ax.get_ylim()[1] * .92, 'New range: n ≥ 32', fontsize=9.5, color='#52606d')
    ax.legend(loc='upper right', ncol=3, frameon=False, fontsize=9)
    ax.grid(axis='y', alpha=.16)
    ax.set_axisbelow(True)
    ax = axes[1]
    keep = np.isin(data['residue_mod30'], RESIDUES)
    x = data['square_n'] + data['position_numerator'] / data['position_denominator']
    y = np.array([RESIDUES.index(int(r)) if r in RESIDUES else -1 for r in data['residue_mod30']])
    ax.axvspan(1, 32, color='#eff1f3', zorder=0)
    for row in range(len(RESIDUES)):
        ax.axhline(row, color='#dbe1e6', linewidth=.8, zorder=0)
    reg = keep & data['regular']; irr = keep & data['irregular']
    ax.scatter(x[reg], y[reg], s=13, color=REGULAR, alpha=.42, edgecolors='none', zorder=2)
    ax.scatter(x[irr], y[irr], s=17 + 10 * data['irregularity_index'][irr], color=IRREGULAR, alpha=.9, edgecolors='white', linewidths=.25, zorder=3)
    for p, offset in [(1093, (7, .43)), (3511, (5, .43))]:
        j = int(np.searchsorted(data['p'], p))
        ax.scatter([x[j]], [y[j]], marker='*', s=170, color=OBSERVED, edgecolors='white', linewidths=.8, zorder=5)
        family = 'irregular' if data['irregular'][j] else 'regular'
        ax.annotate(f'{p}: Wieferich, {family}', xy=(x[j], y[j]), xytext=(x[j] + offset[0], y[j] + offset[1]),
                    fontsize=9, color=OBSERVED, arrowprops={'arrowstyle': '-', 'color': OBSERVED, 'lw': .8},
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .93, 'pad': 1.5}, zorder=6)
    ax.axvline(32, color='#7a8791', linestyle='--', linewidth=1)
    ax.set(xlim=(1, 101), ylim=(-.55, 7.85), yticks=np.arange(8), yticklabels=[str(r) for r in RESIDUES],
           ylabel='Residue p mod 30', xlabel=r'$n+(p-n^2)/(2n+1)$: each unit is one square interval')
    ax.set_title('The prime web: regularity over the eight residue classes', loc='left', pad=13, fontsize=12, color=INK)
    fig.text(.072, .045, 'Orange marker size counts Bernoulli witnesses. The web omits the three residue exceptions 2, 3 and 5.\nThe regular/irregular split is exact in this finite census; it gives no guarantee of future interval coverage.', fontsize=9.4, color='#52606d', linespacing=1.5)
    save_figure(fig, 'Irregular_Legendre_Census')


def null_figure(analysis, draws):
    by_key = {s['key']: s for s in analysis['scenarios']}
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 8.3))
    fig.subplots_adjust(left=.073, right=.98, top=.825, bottom=.145, hspace=.66, wspace=.30)
    fig.suptitle('Do the irregular-prime traces depart from shuffled labels?', x=.073, y=.974, ha='left', fontsize=18, weight='bold', color=INK)
    fig.text(.073, .925, '49,999 shuffles per model and range · prime positions and interval prime counts stay fixed', fontsize=11, color=INK)
    fig.legend(handles=[Line2D([0], [0], color=GLOBAL, lw=2, label='Global shuffle'),
                        Line2D([0], [0], color=STRATIFIED, lw=2, label='Residue + size-band shuffle'),
                        Line2D([0], [0], color=OBSERVED, lw=2, label='Observed')], loc='upper left', bbox_to_anchor=(.064, .90), ncol=3, frameon=False, fontsize=10)
    names = ['Irregular-free intervals', 'Longest consecutive empty run', 'Thin-interval irregular fraction\nminus other-interval fraction']
    for row, n_min in enumerate([6, 32]):
        first = by_key[f'n{n_min}_global']; second = by_key[f'n{n_min}_stratified']
        for column in range(3):
            ax = axes[row, column]
            a = draws[first['key']][:, column].copy(); b = draws[second['key']][:, column].copy()
            obs = first['tests'][column]['observed']
            if column == 2:
                a *= 100; b *= 100; obs *= 100
                bins = np.linspace(min(a.min(), b.min()), max(a.max(), b.max()), 31)
                xlabel = 'Difference (percentage points)'
            else:
                bins = np.arange(-.5, max(a.max(), b.max(), obs) + 1.5, 1)
                xlabel = 'Number of square intervals'
                ax.set_xticks(np.arange(0, int(max(a.max(), b.max(), obs)) + 1))
            for values, color in [(a, GLOBAL), (b, STRATIFIED)]:
                weights = np.ones(len(values)) / len(values)
                ax.hist(values, bins=bins, weights=weights, histtype='stepfilled', color=color, alpha=.14)
                ax.hist(values, bins=bins, weights=weights, histtype='step', color=color, linewidth=1.5)
            ax.axvline(obs, color=OBSERVED, linewidth=1.8, zorder=5)
            ax.set_title(f'n = {n_min}–100 · {names[column]}', fontsize=10.4, pad=11, color=INK)
            ax.set_xlabel(xlabel, fontsize=9)
            if column == 0:
                ax.set_ylabel('Fraction of shuffles', fontsize=9)
            p1 = first['tests'][column]['p_value_mc']; p2 = second['tests'][column]['p_value_mc']
            ax.text(.97, .96, f'Observed {obs:.2f}' if column == 2 else f'Observed {int(obs)}', ha='right', va='top', transform=ax.transAxes, fontsize=9, color=OBSERVED)
            ax.text(.97, .83, f'p: {p1:.3f} global\n   {p2:.3f} stratified', ha='right', va='top', transform=ax.transAxes, fontsize=8.8, color=INK,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .78, 'pad': 2})
            ax.grid(axis='y', alpha=.15); ax.set_axisbelow(True)
    fig.text(.073, .06, 'Raw Monte Carlo p-values are shown; all 12 Holm-adjusted p-values are 1.00. Empty-count and run tests are upper-tail;\nthe fraction comparison is two-sided. These conditional checks assess label organization, not Legendre’s conjecture.', fontsize=9.4, color='#52606d', linespacing=1.5)
    save_figure(fig, 'Irregular_Legendre_Null_Comparisons')


def report(data, census, analysis, validation, draws):
    lookup = {s['key']: s for s in analysis['scenarios']}
    full = lookup['n6_global']; strat = lookup['n6_stratified']; fresh = lookup['n32_global']; fresh_s = lookup['n32_stratified']
    table = []
    labels = ['Irregular-free intervals', 'Longest run of irregular-free intervals', 'Thin-minus-other irregular fraction']
    for a, b in [(full, strat), (fresh, fresh_s)]:
        for j, label in enumerate(labels):
            t = a['tests'][j]; u = b['tests'][j]
            observed = f"{100*t['observed']:+.2f} percentage points" if j == 2 else str(int(t['observed']))
            band = t['null_central_95_percent']
            reference = f"[{100*band[0]:+.2f}, {100*band[1]:+.2f}] pp" if j == 2 else f"[{band[0]:.0f}, {band[1]:.0f}]"
            table.append(f"| {a['n_min']}–100 | {label} | {observed} | {reference} | {t['p_value_mc']:.5f} | {u['p_value_mc']:.5f} |")
    fresh_coverage = float(np.mean(draws[fresh['key']][:, 0] == 0))
    fresh_s_coverage = float(np.mean(draws[fresh_s['key']][:, 0] == 0))
    overlap_lines = []
    for j, mark in enumerate(['Twin members', 'Sophie Germain', 'Safe', 'Wieferich, base 2']):
        m = data['other_properties'][:, j] & (data['p'] > 2)
        overlap_lines.append(f"| {mark} | {int(np.sum(m & data['regular']))} | {int(np.sum(m & data['irregular']))} |")
    report_text = f'''# Irregular primes in Legendre intervals

Experiment v0.1.0 · 2026-09-06 · companion to MTFT 0.26.0

## Finding

This finite census detected no departure from either conditional label-shuffle
model in the specified comparisons. The observed irregular-free intervals, their
longest run, and the regularity mix in thin intervals are all compatible with
these reference models. This does not establish independence or randomness.

Every interval from n = 15 through n = 100 contains both regular and irregular
primes. The absence of irregular-free intervals in the previously unexamined
range n = 32–100 is common under the nulls: {100*fresh_coverage:.2f}% of global
shuffles and {100*fresh_s_coverage:.2f}% of stratified shuffles also have full coverage.
It therefore supplies little evidence for an eventual irregular-prime version
of Legendre. No new prime-gap theorem is claimed.

## Exact census

The range is 1 < p < 101² = 10,201, including {census['prime_count']:,} primes.
Of the {census['odd_prime_count']:,} odd primes, **{census['irregular_count']:,} are
irregular and {census['regular_count']:,} are regular**. The prime 2 is explicitly
unclassified. All {census['witness_pair_count']:,} irregular pairs (p,k) are retained.

| Number of Bernoulli witnesses | Primes |
|---|---:|
| 0, regular | {census['irregularity_index_histogram']['0']} |
| 1 | {census['irregularity_index_histogram']['1']} |
| 2 | {census['irregularity_index_histogram']['2']} |
| 3 | {census['irregularity_index_histogram']['3']} |

There are no prime-free square intervals in this finite range. The irregular-free
interval indices are **1, 2, 3, 4, 5, 9, 13, 14**. For 2 ≤ n ≤ 100 every interval has a
regular prime. Regular and irregular classifications here are finite exact facts;
their respective infinitude statuses are properties of the families.

![Counts and residue traces](Irregular_Legendre_Census.png)

For an odd prime p, the full witness set is

\\[
\\mathcal{{W}}(p)=\\{{k\\in\\{{2,4,\\ldots,p-3\\}}:p\\mid\\operatorname{{num}}(B_k)\\}}.
\\]

The empty set defines regularity. The irregular family is proved infinite;
infinitude of the regular family is open. These definitions and the known global
counting results are discussed by [Luca, Pizarro-Madariaga and Pomerance](https://math.dartmouth.edu/~carlp/irreg.pdf).
Infinitude of a family supplies no guarantee that each short interval contains
one of its members. [Chamberland and Straub](https://arxiv.org/html/2602.22502v1)
review the unresolved consecutive-square problem.

## Conditional experiments

The plan was written before calculating labels above 1,024. The earlier
conversation had already inspected the smaller census, so n = 6–100 is an
exploratory reuse of those data; n = 32–100 is reported separately. This split is
not an external preregistration and does not make the larger research program
confirmatory. The range and tests were not expanded after inspecting p-values.

Both nulls hold every prime value and every interval's prime count fixed.
The global null uniformly permutes the labels and conditions on their total.
The stratified null conditions additionally on irregular totals in each residue
modulo 30 and fixed band floor((n−6)/20). Each selected analysis range is shuffled
separately. The nulls are reference models, not established laws of prime regularity.

If a stratum contains M primes, K irregular labels, and m positions in an
interval, its exact probability of contributing no irregular prime is

\\[
P_0=\\frac{{\\binom{{M-K}}{{m}}}}{{\\binom{{M}}{{m}}}}.
\\]

Multiply across strata for an interval, then sum over intervals to obtain the
expected number of irregular-free intervals. This expectation does not assume
that different intervals are independent. It is {full['expected_irregular_free_count']:.4f}
under the global null and {strat['expected_irregular_free_count']:.4f} under the
stratified null for n = 6–100; the observed count is 3.

Thin intervals are the bottom quarter of m_n/μ_n, where
μ_n = (2n+1)/log(n²+n+1/2), with ties broken by n. This definition uses no
regularity labels. The logarithmic expression is a ranking reference, not a
claimed prime-count formula for each short interval. Fractions are calculated
among primes pooled across the thin and other intervals, respectively.

| Range n | Statistic | Observed | Global central 95% null range | Global p | Stratified p |
|---|---|---:|---|---:|---:|
{chr(10).join(table)}

There are 49,999 shuffles in each of four scenarios, using PCG64 seeds
20260906–20260909. The tests for excess empty intervals and long runs are
upper-tail; the fraction test is two-sided around its analytic conditional
expectation. Ties are included. Monte Carlo p-values use (1 + tail count)/50,000.
**All 12 Holm-adjusted p-values are 1.00.** The plotted 95% ranges describe null
distributions; they are not confidence intervals for a population parameter.

![Conditional comparison distributions](Irregular_Legendre_Null_Comparisons.png)

In the new range the expected number of irregular-free intervals is only
{fresh['expected_irregular_free_count']:.4f} globally and {fresh_s['expected_irregular_free_count']:.4f}
under stratification. Zero observed is therefore unsurprising and the upper-tail
tests have little opportunity to detect anything without an actual empty
interval. Nonsignificance is not an equivalence test.

## Overlaps with the earlier prime families

The previously computed four prime marks reproduce exactly on the overlap.
Among odd primes in this census their intersections are:

| Additional property | Regular | Irregular |
|---|---:|---:|
{chr(10).join(overlap_lines)}

These rows overlap; they must not be added as disjoint categories. In particular:

| Prime | Base-2 Wieferich | Regularity | Bernoulli witness indices |
|---:|---|---|---|
| 1093 | Yes | Regular | None |
| 3511 | Yes | Irregular | 1416, 1724 |

Thus Wieferich status alone does not determine regularity. Two examples do not
establish statistical independence. No infinitude claim is inherited by any
intersection in this table.

## Why the finite arithmetic is checkable

In the ring F_p[t]/(t^(p−2)), form
A(t) = Σ_(j=0)^(p−3) t^j/(j+1)!.
The constant term is 1, so the inverse is unique. From t/(exp(t)−1), its
coefficient of t^k is B_k/k!. All factorials required here are invertible modulo
p. Newton inversion doubles the retained degree using b ← b(2−Ab).

Every convolution uses int64 integers, with the conservative bound
(p−2)(p−1)² < 2^63−1 checked before arithmetic. No floating-point operation
enters a Bernoulli residue, witness decision, or prime count. Floating point is
used for probabilities, statistical summaries and plotting.

The saved even residues, together with B_0=1, B_1=−1/2 and the vanishing odd
Bernoulli numbers above B_1, reconstruct the full inverse. Multiplication by A
checks the entire series identity for every odd prime, including every negative
membership decision. This is a finite computational certificate, not a Lean proof.

For the separate power-sum check, Faulhaber's formula gives, for even
2 ≤ k ≤ p−3,

\\[
\\sum_{{a=1}}^{{p-1}}a^k \\equiv pB_k\\pmod{{p^2}}.
\\]

The terms discarded modulo p² have p-adically invertible denominators in this
range. Compute the left side using modular exponentiation, verify divisibility
by p, and divide to recover B_k modulo p.

| Validation | Coverage |
|---|---:|
| Primality checked by trial division | All {validation['integers_checked_by_trial_division']:,} integers in range |
| Complete modular inverse identities | All {validation['full_modular_series_product_certificates']:,} odd primes |
| Comparison to exact rational Bernoulli numbers | {validation['bernoulli_residues_checked_against_exact_rationals_below_1024']:,} residues below p=1,024 |
| Separate modular-power check | All {validation['all_positive_witnesses_checked_by_independent_power_sums']:,} positive witnesses |
| Separate modular-power negative checks | {validation['negative_entries_checked_by_independent_power_sums']} sampled nonzero residues |
| Exhaustive toy null assignments | {validation['exhaustively_enumerated_toy_null_assignments']} |
| Prime counts partition into the two families plus prime 2 | All 100 intervals |

## What this establishes for the Legendre program

For n ≥ 2, let I_n and R_n count the two families. The identity
L_n = I_n + R_n is exact, and Legendre asks whether L_n ≥ 1 for every n.
The n=1 case is checked directly. Our shuffles condition on L_n, so they cannot
test the event L_n=0: any prime-free interval would stay prime-free in every
shuffle. Their purpose is to test whether regularity contains extra organization
within the already observed primes.

This first experiment found no detectable additional organization in its three
statistics. It supports keeping both families and their witnesses in the trace
catalogue. A useful subsequent research question would require a specified
arithmetic mechanism connecting those witnesses to interval coverage; scaling
the same census alone cannot supply a proof.

## Files and reproduction

See README.md for commands. census.npz retains every even residue with offsets;
prime_witnesses.json records every prime's full witness set. The raw Monte Carlo
statistics, exact hypergeometric probabilities, analysis plan and all scripts
are included. Validation results are in validation_results.json. SHA256SUMS.txt
binds the packaged files, and PROVENANCE.json identifies the prior census.
'''
    (ROOT / 'Irregular_Legendre_Report.md').write_text(report_text)


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK})
    data = np.load(ROOT / 'census.npz')
    census = json.loads((ROOT / 'census_summary.json').read_text())
    analysis = json.loads((ROOT / 'analysis_results.json').read_text())
    validation = json.loads((ROOT / 'validation_results.json').read_text())
    draws = np.load(ROOT / 'permutation_statistics.npz')
    if int(data['n'][-1]) != 100 or any(s['repetitions'] != 49999 for s in analysis['scenarios']):
        raise ValueError('This report renderer describes the fixed 100-interval, 49,999-shuffle design. Reproduce with the default arguments.')
    census_figure(data, census)
    null_figure(analysis, draws)
    report(data, census, analysis, validation, draws)
    provenance = {'date_utc': '2026-09-06', 'experiment': 'MTFT Irregular Legendre v0.1.0',
                  'mtft_companion_version': '0.26.0',
                  'analysis_plan_sha256': hashlib.sha256((ROOT / 'EXPERIMENT_PLAN.md').read_bytes()).hexdigest(),
                  'prior_census': census['prior_census_reconciliation'],
                  'computed_data_sources': 'Local exact integer arithmetic; no downloaded irregular-prime table.',
                  'references': ['https://math.dartmouth.edu/~carlp/irreg.pdf', 'https://arxiv.org/html/2602.22502v1'],
                  'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'matplotlib': matplotlib.__version__}}
    (ROOT / 'PROVENANCE.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps({'figures': ['Irregular_Legendre_Census.png', 'Irregular_Legendre_Null_Comparisons.png'],
                      'report': 'Irregular_Legendre_Report.md', 'validation_status': validation['status']}))


if __name__ == '__main__':
    main()
