#!/usr/bin/env python3
"""Scientific summary of independently specified transport experiments."""
from pathlib import Path
import csv
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
def read(name):
    with (ROOT/name).open(newline='') as f:return list(csv.DictReader(f))

growth=read('monodromy/monodromy_growth.csv')
interval=read('root_checks/interval_results.csv')
plt.rcParams.update({'font.size':11,'axes.spines.top':False,
    'axes.spines.right':False,'figure.dpi':160})
fig,(left,right)=plt.subplots(1,2,figsize=(12.8,5.5))
fig.subplots_adjust(left=.075,right=.97,top=.77,bottom=.24,wspace=.3)
fig.suptitle('Transport growth depends on the geometry of the path',
    x=.075,y=.98,ha='left',fontsize=18,weight='bold')
fig.text(.075,.88,'A fixed loop can amplify periods; a chaotic interval map can have zero transport exponent.',
    fontsize=11,color='#4f6072')
labels={'elliptic_zero':(r'$S^n$: bounded','#757575'),
    'parabolic_minus_one':(r'$U^n$: linear','#2b7a78'),
    'parabolic_alternation':(r'$(US)^n$: linear','#6875ab'),
    'hyperbolic_word':(r'$(U^2S)^n$: exponential','#bd5e35')}
for word in dict.fromkeys(r['word'] for r in growth):
    rows=[r for r in growth if r['word']==word]
    label,color=labels.get(word,(word,'#6875ab'))
    left.plot([int(r['repetitions']) for r in rows],
        [math.log(float(r['operator_norm'])) for r in rows],'o-',
        label=label,color=color,markersize=4)
left.set(xlabel='Repetitions of the specified closed path',
    ylabel=r'$\log\|M^n\|_2$',title='Closed paths around singular parameters')
left.grid(alpha=.18);left.legend(frameon=False,fontsize=9)
seeds=list(dict.fromkeys(r['seed'] for r in interval))
for i,seed in enumerate(seeds):
    rows=[r for r in interval if r['seed']==seed]
    right.loglog([int(r['steps']) for r in rows],
        [float(r['transport_finite_exponent']) for r in rows],'o-',
        label=f'Transport, initial condition {i+1}',markersize=4)
right.axhline(math.log(2),color='#bd5e35',linestyle='--',
    label=r'Base reference: $\log 2$')
right.set(xlabel='Iterations of the imposed logistic base map',
    ylabel='Finite exponent',title=r'Flat transport on $c\in[-3,-2]$')
right.grid(alpha=.18);right.legend(frameon=False,fontsize=8,loc='upper right',bbox_to_anchor=(1,.8))
fig.text(.075,.105,'Left: exact monodromy characters establish the growth classes; plotted norms depend on the chosen frame.',
    fontsize=9,color='#4f6072')
fig.text(.075,.06,'Right: telescoping proves the transport exponent is zero. Numerical samples illustrate the theorem.',
    fontsize=9,color='#4f6072')
fig.savefig(ROOT/'transport_comparison.png',facecolor='white')
plt.close(fig)
print(ROOT/'transport_comparison.png')
