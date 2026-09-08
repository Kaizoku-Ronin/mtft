"""Scientific plots of finite observations with theorem status kept explicit."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from prime_families import RESIDUES,LABELS

ROOT=Path(__file__).resolve().parent
DATA=np.load(ROOT/'prime_traces.npz')
COLORS=['#147f9b','#d27828','#528e4b','#963f9f']
SYMBOLS=['o','^','s','*']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'savefig.dpi':180,'figure.facecolor':'white'})


def save(fig,name):
    fig.savefig(ROOT/(name+'.png'),bbox_inches='tight',pad_inches=.15)
    fig.savefig(ROOT/(name+'.svg'),bbox_inches='tight',pad_inches=.15)
    plt.close(fig)


def weave():
    fig,axes=plt.subplots(2,1,figsize=(9,7.4),layout='constrained')
    for ax,(start,stop) in zip(axes,((31,36),(57,62))):
        lo,hi=start*start,stop*stop
        ids=np.flatnonzero((DATA['p']>lo)&(DATA['p']<hi))
        p=DATA['p'][ids];rows=np.array([RESIDUES.index(int(r)) for r in DATA['residue_mod30'][ids]])
        for j in range(8):ax.axhline(j,color='#d2dcdf',lw=1,zorder=0)
        squares=np.arange(start,stop+1)**2
        for i in range(len(squares)-1):
            if i%2==0:ax.axvspan(squares[i],squares[i+1],color='#f2f5f6',zorder=-1)
        ax.scatter(p,rows,s=15,color='#82969b',zorder=2)
        for j,(offset,size) in enumerate(zip((-.24,-.08,.08,.24),(28,40,28,140))):
            seen=DATA['properties'][ids,j]
            ax.scatter(p[seen],rows[seen]+offset,s=size,marker=SYMBOLS[j],color=COLORS[j],zorder=4+j)
            if j==3:
                for x,y in zip(p[seen],rows[seen]+offset):
                    ax.annotate(str(x),xy=(x,y),xytext=(-8,0),textcoords='offset points',fontsize=10,
                                ha='right',va='center',fontweight='bold',color=COLORS[j])
        ax.set(xlim=(lo,hi),ylim=(7.55,-.6),yticks=range(8),yticklabels=[str(r) for r in RESIDUES],
               xticks=squares,xticklabels=[f'{n}²\n{n*n:,}' for n in range(start,stop+1)],
               ylabel='Residue row: p mod 30',xlabel='Prime p; shaded bands are consecutive-square intervals',
               title=f'Five square intervals: n = {start}–{stop-1}')
        ax.tick_params(axis='x',length=0);ax.spines['left'].set_visible(False)
    handles=[Line2D([],[],ls='',marker='o',color='#82969b',ms=4,label='Every prime in the window')]
    handles += [Line2D([],[],ls='',marker=m,color=c,ms=8 if m=='*' else 5,label=l) for m,c,l in zip(SYMBOLS,COLORS,LABELS)]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False,fontsize=10)
    fig.suptitle('Prime traces on residue classes modulo 30\nRows: infinitude proved · Colored labels: infinitude open',fontsize=14)
    save(fig,'Prime_Trace_Web')


def growth():
    fig,axes=plt.subplots(1,3,figsize=(14,4.6),layout='constrained')
    n=DATA['n'];x=DATA['upper_square'];total=DATA['cumulative_total']
    cmap=plt.get_cmap('tab10')
    for i,r in enumerate(RESIDUES):
        axes[0].plot(n,DATA['cumulative_residue'][:,i],lw=1,color=cmap(i),label=str(r))
    axes[0].set(title='Eight families proved infinite',ylabel='Observed prime count in each residue class',xlabel='Last square interval n')
    axes[0].legend(title='p mod 30',ncol=4,fontsize=8,frameon=False,loc='upper left')
    for j,(label,color) in enumerate(zip(LABELS,COLORS)):
        y=DATA['cumulative_types'][:,j].astype(float);y[y==0]=np.nan
        axes[1].plot(n,y,color=color,lw=1.8,label=label)
        axes[2].plot(n,100*y/total,color=color,lw=1.8)
    axes[1].set(yscale='log',title='Four families with infinitude open',ylabel='Observed cumulative members [log scale]',xlabel='Last square interval n')
    axes[2].set(yscale='log',title='Rarity relative to all primes',ylabel='Members / all primes [%; log scale]',xlabel='Last square interval n')
    axes[1].legend(frameon=False,loc='center right',fontsize=8)
    for ax in axes:
        ax.set_xlim(1,n[-1]);ax.grid(alpha=.14)
    fig.suptitle(f'Finite traces through n = {n[-1]:,}, prime bound {x[-1]:,} · no extrapolation',fontsize=15)
    save(fig,'Prime_Trace_Growth')


if __name__=='__main__':weave();growth();print('Created trace web and growth figures as PNG and SVG.')
