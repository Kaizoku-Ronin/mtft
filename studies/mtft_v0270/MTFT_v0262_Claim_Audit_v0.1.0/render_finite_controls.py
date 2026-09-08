"""Plot finite diffusion and exact-DOS thermal controls without exponent fitting."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent


def main():
    data=json.loads((ROOT/'finite_graph_results.json').read_text())['levels']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titleweight':'bold'})
    fig,axes=plt.subplots(1,2,figsize=(12,5.6))
    fig.subplots_adjust(left=.075,right=.98,bottom=.19,top=.76,wspace=.27)
    fig.suptitle('Finite geometry controls',x=.075,y=.97,ha='left',fontsize=23,weight='bold')
    fig.text(.075,.90,'MTFT 0.26.2 · exact graph inputs and numerical diagnostics · no critical exponent fitted',fontsize=11,color='#566475')
    for N,label,color in [('6','N=6 and N=11: identical graph','#8970ad'),('143','N=143: 56 vertices','#008b8b')]:
        row=data[N]
        trace=row['heat_trace']
        axes[0].semilogx([x['t'] for x in trace],[x['d_eff'] for x in trace],color=color,label=label,lw=2)
        betas=np.linspace(0,2.4,401)
        D=np.array(row['DOS_exact'],float); k=np.flatnonzero(D)
        energies=2*k-row['edges']; C=[]
        for beta in betas:
            logw=np.log(D[k])-2*beta*k
            w=np.exp(logw-max(logw));w/=w.sum()
            mean=w@energies
            C.append(beta*beta*(w@((energies-mean)**2))/row['vertices'])
        axes[1].plot(betas,C,color=color,lw=2,label=label)
    axes[0].set(title='Running diffusion dimension changes with scale',xlabel='Diffusion time t',ylabel=r'$d_{\rm eff}(t)$, including the zero mode',ylim=(-.05,2.5))
    axes[0].grid(alpha=.2);axes[0].legend(frameon=False,fontsize=9,loc='upper right')
    axes[1].set(title='Finite thermal peaks remain smooth',xlabel=r'Inverse temperature $\beta$  (J=1)',ylabel='Heat capacity per spin')
    axes[1].grid(alpha=.2);axes[1].legend(frameon=False,fontsize=9,loc='upper right')
    fig.text(.075,.07,'Levels 6 and 11 have surface genera 0 and 1, yet identical zero-field Ising thermodynamics.',fontsize=11,color='#566475')
    fig.text(.075,.025,'A limiting graph family, metric scale, and coupling prescription are needed before discussing universal critical exponents.',fontsize=10,color='#566475')
    for extension in ['png','svg']:
        fig.savefig(ROOT/f'Finite_Geometry_Controls.{extension}',dpi=180,facecolor='white')
    plt.close(fig)
    print('Rendered Finite_Geometry_Controls.png and .svg')


if __name__=='__main__':
    main()
