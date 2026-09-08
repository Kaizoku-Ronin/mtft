"""Scientific figure generated solely from the recorded numerical results."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter
from common import ROOT, read


def main():
    charge = read('charge_results.json')
    local = read('local_qed_results.json')
    controls = read('photon_atomic_results.json')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titleweight':'bold','axes.labelcolor':'#263445',
                         'text.color':'#172333','xtick.color':'#46556a','ytick.color':'#46556a'})
    colors = ['#007f82','#e1942c','#7964b6']
    fig,axes = plt.subplots(2,2,figsize=(13,10))
    fig.subplots_adjust(left=.09,right=.97,bottom=.12,top=.85,wspace=.31,hspace=.47)
    fig.suptitle('MTFT electron–photon investigation',fontsize=23,x=.09,y=.97,ha='left',weight='bold')
    fig.text(.09,.925,'Arithmetic internal operators; spacetime and QED dynamics supplied explicitly.',fontsize=12,color='#536477')
    ax=axes[0,0]
    for epsilon,color in zip([.05,.2,.4],colors):
        rows=[r for r in charge['scans'] if r['seed']==2026090700 and r['epsilon']==epsilon]
        ax.plot([r['charge_difference'] for r in rows],[r['commutator_norm'] for r in rows],
                marker='o',color=color,label=f'Mixing strength {epsilon:g}')
    ax.set(title='A  Neutral mixing constrains charge',xlabel='Old–quartic charge difference',
           ylabel=r'$\|[Q,M]\|_F$  (dimensionless)')
    ax.legend(frameon=False,fontsize=10,loc='upper left')
    ax.text(.02,.95,'16 seeds agree with the exact pattern',transform=ax.transAxes,va='top',fontsize=9,color='#536477')
    ax.set_ylim(-.025,.87); ax.set_xticks([0,.5,1,2]); ax.grid(alpha=.18)
    # Move legend away from the explanatory line.
    ax.legend(frameon=False,fontsize=9,loc='upper left',bbox_to_anchor=(0,.86))
    ax=axes[0,1]
    cases=['common_charge','equal_mixed_charges','different_mixed_charges']
    values=[[r['relative_covariance_defect'] for r in local['runs'] if r['case']==c] for c in cases]
    values.append([r['relative_spurion_covariance_defect'] for r in local['runs'] if r['case']=='different_mixed_charges'])
    for i,(v,color) in enumerate(zip(values,[colors[0],colors[0],'#bb4438',colors[2]])):
        ax.scatter(np.linspace(i-.1,i+.1,len(v)),v,color=color,s=35,zorder=3)
    ax.axhline(1e-10,color='#8491a1',linestyle='--',linewidth=1,label='Numerical gate')
    ax.set_yscale('log'); ax.set_ylim(6e-17,.1)
    ax.set_xticks(range(4),['Common\ncharge','Equal mixed\ncharges','Different mixed\ncharges','Transform\nmass source'])
    ax.tick_params(axis='x',labelsize=9)
    ax.set(title='B  A local spacetime test detects failure',ylabel='Relative fermion covariance defect')
    ax.legend(frameon=False,fontsize=9,loc='center right'); ax.grid(axis='y',alpha=.18)
    ax.text(.03,.96,r'Periodic $3^4$ lattice · 4,212 complex entries',transform=ax.transAxes,va='top',fontsize=9,color='#536477')
    ax=axes[1,0]
    for mode,label,color in [([1,0,0],'Modes (1,0,0) / (1,1,0)',colors[0]),([1,2,3],'Mode (1,2,3)',colors[1])]:
        rows=[r for r in controls['maxwell']['records'] if r['mode']==mode]
        ax.loglog([r['points_per_axis'] for r in rows],[r['relative_dispersion_error'] for r in rows],
                  marker='o',color=color,label=label)
    ax.set(title='C  Maxwell propagation converges',xlabel='Spatial grid points per direction',
           ylabel=r'Relative error in $\omega/|k|$')
    ax.set_xticks([16,32,64,128]); ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.set_ylim(7e-5,.09)
    ax.legend(frameon=False,fontsize=9,loc='lower left'); ax.grid(which='both',alpha=.18)
    ax.text(.03,.96,'Two transverse modes at each nonzero sample',transform=ax.transAxes,va='top',fontsize=9,color='#536477')
    ax=axes[1,1]
    for n,ell,label,color in [(1,0,'1s',colors[0]),(2,1,'2p',colors[1]),(3,2,'3d',colors[2])]:
        rows=[r for r in controls['coulomb']['records'] if r['n']==n and r['ell']==ell]
        ax.loglog([r['spacing_bohr'] for r in rows],[r['relative_continuum_error'] for r in rows],
                  marker='o',color=color,label=label)
    ax.set(title='D  Coulomb levels converge',xlabel='Radial step / Bohr radius',
           ylabel='Relative energy error vs continuum')
    ax.set_xticks([.0125,.025,.05,.1],['0.0125','0.025','0.05','0.1'])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(frameon=False,fontsize=10,loc='lower right'); ax.grid(which='both',alpha=.18)
    ax.text(.03,.96,'Point nucleus, infinite nuclear mass',transform=ax.transAxes,va='top',fontsize=9,color='#536477')
    fig.text(.09,.045,'A–B: compatibility constraints on trial Hecke mixing. C–D: numerical controls for imported electromagnetic laws.',
             fontsize=10,color='#536477')
    fig.text(.09,.022,'7 September 2026  •  Residuals are numerical diagnostics; the atomic solve does not resolve parts-per-billion coupling shifts.',
             fontsize=9,color='#536477')
    for suffix in ['png','svg']:
        fig.savefig(ROOT/f'Electron_Photon_Results.{suffix}',dpi=180,facecolor='white')
    plt.close(fig)
    print('Rendered Electron_Photon_Results.png and .svg')


if __name__ == '__main__':
    main()
