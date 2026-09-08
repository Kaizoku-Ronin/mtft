#!/usr/bin/env python3
import json
from pathlib import Path
from math import factorial
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dyadic_bounds import rational

ROOT=Path(__file__).resolve().parent
cert=json.loads((ROOT/'BOUND_CERTIFICATE.json').read_text())
data=np.load(ROOT/'short_horizon_pointwise.npz')
ts=data['times'];measured=data['relative_error']
row=next(r for r in cert['rows'] if rational(r['horizon'])==1 and r['order']==12)
den=float(rational(cert['denominator_lower']))
coeff=[float(rational(x)) for x in cert['difference_norm_uppers']]
bound=(sum(coeff[k]*ts**k/factorial(k) for k in range(13))+
       float(rational(row['column_tail_upper']))*ts**13)/den
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(11.5,4.8),gridspec_kw={'width_ratios':[1.25,1]})
for ax in axes:
    ax.axhline(1,color='#9b333b',ls='--',lw=1.2,label='1% target')
    ax.axvspan(0,.5,color='#2b927f',alpha=.12)
    ax.plot(ts,100*bound,color='#7a5cb5',lw=2,label='Proved envelope (order 12)')
    ax.plot(ts,100*measured,color='#176f8c',lw=2,label='Sampled matrix error')
    ax.set_xlabel('Dimensionless evolution parameter t')
    ax.grid(alpha=.18)
axes[0].set(xlim=(0,1),ylim=(0,12),ylabel='Relative Frobenius error (%)',title='The guarantee has a limited horizon')
axes[1].set(xlim=(0,.5),ylim=(0,1.1),title='Certified interval: 0 ≤ t ≤ 0.5')
axes[0].legend(loc='upper left',frameon=False,fontsize=9)
axes[1].annotate('Bound: 0.602%\nSampled max: 0.199%',xy=(.5,100*bound[500]),
                 xytext=(.16,.72),arrowprops={'arrowstyle':'-','color':'#555'},fontsize=10)
fig.suptitle('T2 memory reduction: five complement directions → four',fontsize=16,y=.98)
fig.text(.5,.018,'Exact rational certificate for frozen Hermitian matrices · No physical time identification',ha='center',fontsize=10,color='#555')
fig.tight_layout(rect=(0,.07,1,.94))
for ext in ['png','svg']:
    fig.savefig(ROOT/f'T2_Short_Horizon.{ext}',dpi=180)
