"""Standalone scientific figures from the saved finite-model experiments."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/'boundary_results.json').read_text())
COLORS={'ferro':'#007f86','one_reversed_bond':'#c2552f'}
LABELS={'ferro':'All bonds ferromagnetic','one_reversed_bond':'One boundary bond reversed'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
                     'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.dpi':180,'figure.facecolor':'white','axes.facecolor':'white'})


def save(fig,name):
    fig.savefig(ROOT/(name+'.png'),bbox_inches='tight',pad_inches=.18)
    fig.savefig(ROOT/(name+'.svg'),bbox_inches='tight',pad_inches=.18)
    plt.close(fig)


def temperature():
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.4),layout='constrained')
    for name,rows in DATA['temperature_rows'].items():
        b=np.array([r['beta'] for r in rows]);c=COLORS[name]
        axes[0].plot(b,[r['all_boundary_MI'] for r in rows],color=c,lw=2,label=LABELS[name])
        peak=DATA['MI_peaks'][name]
        axes[0].plot(peak['beta'],peak['mutual_information_bits'],'o',color=c)
        axes[1].plot(b,[r['H_A'] for r in rows],color=c,lw=1.5,ls='--')
        axes[1].plot(b,[r['all_boundary_conditional_entropy'] for r in rows],color=c,lw=2)
        axes[2].plot(b,[100*r['all_boundary_spin_accuracy'] for r in rows],color=c,lw=2)
    axes[0].axhline(1,color='#777777',lw=.8,ls=':')
    axes[0].annotate('Clean peak: 1.571 bits\nβ ≈ 0.612',xy=(.612,1.571),xytext=(1.27,1.65),
                     arrowprops={'arrowstyle':'-','color':'#555555'},fontsize=9)
    axes[0].set(title='Boundary information peaks',ylabel='Mutual information I(A : B) [bits]',ylim=(-.05,1.96))
    axes[1].set(title='Cold interiors become predictable',ylabel='Interior entropy [bits]',ylim=(-.3,14.6))
    axes[1].legend(handles=[Line2D([],[],color='#444444',ls='--',label='Before observing B'),
                           Line2D([],[],color='#444444',label='After observing B, averaged')],loc='upper right',fontsize=8,frameon=False)
    axes[2].set(title='Predict each interior spin',ylabel='Bayes accuracy, averaged [%]',ylim=(48,102))
    axes[2].legend(handles=[Line2D([],[],color=COLORS[k],lw=2,label=LABELS[k]) for k in COLORS],loc='lower right',fontsize=8,frameon=False)
    for ax in axes:
        ax.set(xlabel='Inverse temperature β [J = 1]',xlim=(0,3));ax.grid(alpha=.15)
    fig.suptitle('Eight boundary spins observing a fourteen-spin interior · MTFT N = 143',fontsize=15)
    save(fig,'Boundary_Temperature_Response')


def sensor_selection():
    fig,(ax,graph)=plt.subplots(1,2,figsize=(11.5,5.2),layout='constrained',gridspec_kw={'width_ratios':[1.05,1]})
    rows=DATA['sensor_selection_beta064']['ferro'];k=np.arange(9)
    ax.fill_between(k,[r['worst_MI'] for r in rows],[r['best_MI'] for r in rows],color='#007f86',alpha=.13,label='Range across all sensor subsets')
    ax.plot(k,[r['best_MI'] for r in rows],'-o',color=COLORS['ferro'],ms=4,label='Best subset at each budget')
    ax.plot(k,[r['greedy_MI'] for r in rows],'--s',color=COLORS['one_reversed_bond'],ms=4,label='Greedy additions')
    ax.set(xlabel='Number of observed boundary spins',ylabel='Mutual information I(A : observed B) [bits]',
           title='Evaluate every one of the 256 subsets',xticks=k,xlim=(-.15,8.2),ylim=(-.04,1.68))
    ax.grid(alpha=.15);ax.legend(frameon=False,loc='upper left',fontsize=9)
    part=json.loads((ROOT/'partition.json').read_text());pos=np.load(ROOT/'graph143_drawing_positions.npy')
    A=set(part['A']);B=set(part['B']);pair={6,22};singles={24,40}
    for u,v in part['edges']:
        if u==v:continue
        crossing=(u in A and v in B) or (v in A and u in B)
        graph.plot(pos[[u,v],0],pos[[u,v],1],color='#9ba7a8' if crossing else '#d5dddd',lw=1.2 if crossing else .7,zorder=1)
    for v in range(56):
        x,y=pos[v]
        if v in B:
            color='#007f86' if v in pair else '#c2552f' if v in singles else '#d9e1e2'
            graph.scatter(x,y,s=130 if v in pair|singles else 65,marker='s',c=color,zorder=4)
            dx=-.10 if v in {6,36,40,53} else .09
            graph.text(x+dx,y+.045,str(v),ha='right' if dx<0 else 'left',va='bottom',fontsize=10,zorder=5)
        else:graph.scatter(x,y,s=37 if v in A else 12,c='#546b6f' if v in A else '#c1cccc',zorder=3)
    graph.set_aspect('equal');graph.axis('off');graph.set_title('The best pair excludes the best singles')
    graph.legend(handles=[Line2D([],[],marker='s',ls='',color='#007f86',label='Best pair: 6 + 22'),
                          Line2D([],[],marker='s',ls='',color='#c2552f',label='Tied best singles: 24 or 40')],
                 loc='lower center',bbox_to_anchor=(.5,-.07),frameon=False,fontsize=9)
    fig.suptitle('Sensor placement · all ferromagnetic bonds · β = 0.64',fontsize=15)
    save(fig,'Boundary_Sensor_Placement')


def outcomes():
    fig,(ax,defect)=plt.subplots(1,2,figsize=(11.5,4.8),layout='constrained',gridspec_kw={'width_ratios':[1.2,1]})
    examples=[r for r in DATA['observation_examples_beta064'] if r['model']=='ferro']
    prior=next(r for r in examples if r['mask']==0)
    aligned=next(r for r in examples if r['mask']==255 and r['values']==255)
    mixed=next(r for r in examples if r['mask']==255 and r['values']==85)
    row=next(r for r in DATA['temperature_rows']['ferro'] if r['beta']==.64)
    vals=[prior['H_A_given_this_observation'],row['all_boundary_conditional_entropy'],aligned['H_A_given_this_observation'],mixed['H_A_given_this_observation']]
    labels=['No sensors','All sensors,\naveraged','All +1\nP = 21.582%','Alternating bits\nP = 0.0168%']
    bars=ax.bar(range(4),vals,color=['#8b9da0','#007f86','#007f86','#c2552f'],width=.65)
    ax.bar_label(bars,fmt='%.3f',padding=4)
    ax.set(xticks=range(4),xticklabels=labels,ylabel='Interior entropy [bits]',ylim=(0,10),title='A rare outcome can increase uncertainty')
    ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    for name in COLORS:
        r=next(r for r in DATA['temperature_rows'][name] if r['beta']==.64)
        y=0 if name=='ferro' else 1
        defect.barh(y,r['all_boundary_MI'],color=COLORS[name],height=.55)
        defect.barh(y,r['all_boundary_conditional_entropy'],left=r['all_boundary_MI'],color=COLORS[name],height=.55,alpha=.23)
        defect.text(r['all_boundary_MI']/2,y,f"{r['all_boundary_MI']:.3f}",ha='center',va='center',color='white',fontsize=9)
        defect.text(r['all_boundary_MI']+r['all_boundary_conditional_entropy']/2,y,f"{r['all_boundary_conditional_entropy']:.3f}",ha='center',va='center',fontsize=10)
    defect.set(yticks=[0,1],yticklabels=['Clean','Reversed\nboundary bond'],xlabel='Interior entropy, before observing B [bits]',
               xlim=(0,8.8),ylim=(1.7,-.7),title='More information, more uncertainty left')
    defect.legend(handles=[Line2D([],[],color='#657d80',lw=8,label='Information revealed, averaged'),
                           Line2D([],[],color='#d3dfe0',lw=8,label='Remaining uncertainty, averaged')],
                  loc='lower left',frameon=False,fontsize=9)
    defect.grid(axis='x',alpha=.15);defect.set_axisbelow(True)
    fig.suptitle('Boundary observations and a model defect · β = 0.64',fontsize=15)
    save(fig,'Boundary_Outcomes_and_Defect')


if __name__=='__main__':
    temperature();sensor_selection();outcomes()
    print('Created three scientific figures as PNG and SVG.')
