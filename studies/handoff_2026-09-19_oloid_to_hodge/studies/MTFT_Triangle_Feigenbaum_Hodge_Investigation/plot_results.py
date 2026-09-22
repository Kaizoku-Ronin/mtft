#!/usr/bin/env python3
"""Render the two numerical checks; these do not identify their operators."""
import csv
from pathlib import Path
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent

def read(relative):
    with (ROOT / relative).open(newline='') as stream:
        return list(csv.DictReader(stream))

scaling = [r for r in read('scaling/scaling_results.csv') if r['delta_ratio']]
periods = read('geometry/period_results.csv')
generation = [int(r['n']) for r in scaling]
delta_error = [abs(float(r['delta_ratio']) - 4.66920160910299) for r in scaling]
alpha_error = [abs(float(r['alpha_magnitude_ratio']) - 2.50290787509589) for r in scaling]
x = [float(r['log_inverse_epsilon']) for r in periods]
y = [float(r['period_integral']) for r in periods]
asymptote = [v / (2 * math.sqrt(2)) + math.log(8) / math.sqrt(2) for v in x]

plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 160})
fig, (left, right) = plt.subplots(1, 2, figsize=(12, 5.2))
fig.subplots_adjust(left=.08, right=.97, top=.77, bottom=.21, wspace=.29)
fig.suptitle('Two checks for the chosen quadratic family', x=.08, ha='left',
             y=.97, fontsize=19, weight='bold')
fig.text(.08, .865, r'$f_c(z)=z^2+c$  |  Iteration scaling and an elliptic period',
         color='#526171', fontsize=12)

left.semilogy(generation, delta_error, 'o-', color='#176b8b', label=r'$|\delta_n-\delta|$')
left.semilogy(generation, alpha_error, 's-', color='#bd5e35', label=r'$|\alpha_n-\alpha|$')
left.set(xlabel=r'Period-doubling generation $n$ (period $2^n$)',
         ylabel='Absolute difference from reference limit',
         title='Finite ratios approach Feigenbaum limits')
left.set_xticks([2,4,6,8,10])
left.grid(alpha=.2, axis='y')
left.legend(frameon=False)

right.plot(x, asymptote, color='#bd5e35', linewidth=1.8,
           label=r'$\log(1/\epsilon)/(2\sqrt{2})+\log(8)/\sqrt{2}$')
right.scatter(x, y, facecolors='white', edgecolors='#176b8b', s=45,
              zorder=3, label='AGM period calculation')
right.set(xlabel=r'$\log(1/\epsilon)$, with $c=-1-\epsilon$',
          ylabel=r'Branch-to-branch period integral $I(c)$',
          title='The first elliptic period grows logarithmically')
right.grid(alpha=.2, axis='y')
right.legend(frameon=False, loc='upper left', fontsize=9)
fig.text(.08, .08, '60/90-digit calculations; numerical diagnostics, not interval certificates.',
         fontsize=10, color='#526171')
fig.text(.08, .035, 'The two tests do not identify a common scaling operator or determine the MTFT parameter.',
         fontsize=10, color='#526171')
fig.savefig(ROOT / 'investigation.png', facecolor='white')
plt.close(fig)
print(ROOT / 'investigation.png')
