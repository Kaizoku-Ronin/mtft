"""Scientific figures and a labeled equilibrium animation from stored experiment data."""
from __future__ import annotations
from collections import Counter
import csv
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle, Rectangle
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ising_exact import thermodynamics
from playground import bilayer_observables

ROOT=Path(__file__).resolve().parent
TEAL='#087e8b';ORANGE='#cd6a25';PURPLE='#7956a6';DARK='#263543';GRAY='#bac5cd'


def read(name):return json.loads((ROOT/name).read_text())


def save(fig,name):
    fig.savefig(ROOT/(name+'.png'),dpi=185,facecolor='white')
    fig.savefig(ROOT/(name+'.svg'),facecolor='white')
    plt.close(fig)


def header(fig,title,subtitle,left=.075):
    fig.suptitle(title,x=left,y=.965,ha='left',fontsize=21)
    fig.text(left,.905,subtitle,fontsize=11,color='#53616e')


def graph_layout(n,edges,steps=900):
    """Deterministic spring drawing; no metric or physical meaning assigned."""
    rng=np.random.default_rng(143);pos=rng.uniform(-1,1,(n,2))
    nonloops=np.array([(u,v) for u,v in edges if u!=v],dtype=int)
    k=math.sqrt(4/n)
    for step in range(steps):
        delta=pos[:,None,:]-pos[None,:,:]
        d2=np.maximum(np.sum(delta*delta,axis=2),1e-8)
        disp=np.sum(delta*(k*k/d2)[:,:,None],axis=1)
        u,v=nonloops.T;difference=pos[u]-pos[v]
        distance=np.maximum(np.linalg.norm(difference,axis=1),1e-8)
        force=difference*(distance/k)[:,None]
        np.add.at(disp,u,-force);np.add.at(disp,v,force)
        disp-=.025*pos
        temp=.075*(1-step/steps)+.0003
        norm=np.maximum(np.linalg.norm(disp,axis=1),1e-12)
        pos+=disp*(np.minimum(norm,temp)/norm)[:,None]
        pos-=pos.mean(axis=0)
    return pos/np.max(np.abs(pos))


def draw_graph(ax,pos,edges,colors,sizes=45,special_edge=None):
    for i,(u,v) in enumerate(edges):
        special=i==special_edge
        color=PURPLE if special else GRAY
        if u==v:
            direction=pos[u]/max(np.linalg.norm(pos[u]),1e-8)
            ax.add_patch(Circle(pos[u]+.05*direction,.05,fill=False,color=color,lw=1))
        else:
            ax.plot(*pos[[u,v]].T,color=color,lw=3.8 if special else .9,zorder=2 if special else 1,alpha=1 if special else .8)
    ax.scatter(*pos.T,c=colors,s=sizes,edgecolors='white',linewidths=.5,zorder=3)
    ax.set_aspect('equal');ax.set_xlim(pos[:,0].min()-.14,pos[:,0].max()+.14)
    ax.set_ylim(pos[:,1].min()-.14,pos[:,1].max()+.14);ax.axis('off')


def geometry_figure():
    with np.load(ROOT/'field_geometry_grid.npz') as z:
        beta,eta=z['beta'],z['eta'];mag=z['magnetization_per_spin'];K=z['gaussian_curvature']
    fig,ax=plt.subplots(1,2,figsize=(13.2,6.3))
    fig.subplots_adjust(left=.075,right=.92,bottom=.23,top=.78,wspace=.40)
    header(fig,'A field changes the statistical geometry','N = 143  ·  56 spins  ·  joint energy–magnetization counts are exact integers')
    a=ax[0].pcolormesh(beta,eta,mag.T,shading='auto',cmap='RdBu_r',vmin=-1,vmax=1,rasterized=True)
    cb=fig.colorbar(a,ax=ax[0],pad=.025,fraction=.05);cb.set_label('Mean magnetization per spin')
    ax[0].set_title('Field response',loc='left',pad=12)
    a=ax[1].pcolormesh(beta,eta,K.T,shading='auto',cmap='PuOr_r',norm=TwoSlopeNorm(vmin=-2,vcenter=0,vmax=.1),rasterized=True)
    ax[1].contour(beta,eta,K.T,levels=[0],colors='#344454',linewidths=1)
    cb=fig.colorbar(a,ax=ax[1],pad=.025,fraction=.05,ticks=[-2,-1,0,.05,.1]);cb.set_label('Gaussian curvature K')
    ax[1].set_title('Fisher curvature; dark line marks K = 0',loc='left',pad=12)
    for label,et in [('1',0),('2',.03)]:
        ax[1].plot(.64,et,'o',ms=7,mfc='white',mec=DARK,zorder=5)
        ax[1].annotate(label,(.64,et),xytext=(8,1),textcoords='offset points',fontsize=10,fontweight='bold',color=DARK)
    for a in ax:
        a.set(xlabel='Inverse temperature β',ylabel='Natural field η = βh',xlim=(.05,1),ylim=(-.08,.08))
        a.set_yticks([-.08,-.04,0,.04,.08])
    fig.text(.075,.12,'At β = 0.64:  ① η = 0 gives K = +0.08718;  ② η = 0.03 gives K = −0.05531.',fontsize=11)
    fig.text(.075,.055,'Metric: full Cov(S, M), with S = Σedges σuσv and M = Σspins σv. Curves and colors are numerical finite sums.\nSelected points agree across determinant, Brioschi, and Christoffel–Riemann formulas at 60-digit working precision.',fontsize=9.5,color='#53616e')
    save(fig,'MTFT_Field_Geometry')


def topology_layers_figure():
    twists=next(r for r in read('surface_twists.json')['levels'] if r['N']==143)
    layer=next(r for r in read('bilayer.json')['levels'] if r['N']==35)
    beta=np.linspace(.02,1.4,200)
    z0=np.array([thermodynamics(twists['baseline_counts'],56,84,b)['log_Z'] for b in beta])
    fig,ax=plt.subplots(1,3,figsize=(15.3,6.0))
    fig.subplots_adjust(left=.065,right=.975,top=.77,bottom=.25,wspace=.34)
    header(fig,'Two ways to add structure','Flat surface twists on N = 143; an explicit second spin field on N = 35',left=.065)
    for t in twists['basis_twists']:
        vals=(z0-np.array([thermodynamics(t['counts'],56,84,b)['log_Z'] for b in beta]))/beta
        ax[0].plot(beta,vals,color=GRAY,lw=.8,alpha=.55)
    for idx,col in [(0,TEAL),(6,ORANGE),(25,PURPLE)]:
        t=twists['basis_twists'][idx]
        vals=(z0-np.array([thermodynamics(t['counts'],56,84,b)['log_Z'] for b in beta]))/beta
        ax[0].plot(beta,vals,lw=2.1,color=col,label=f'Basis {idx}; min. cost {t["min_frustrated"]}')
    ax[0].set(title='Selected twist costs',xlabel='Inverse temperature β',ylabel='Free-energy cost ΔF / J',xlim=(0,1.4))
    ax[0].legend(frameon=False,fontsize=9,loc='upper left');ax[0].grid(axis='y',alpha=.15)
    C=np.array(twists['two_cycle_projection']['class_even_counts'],dtype=float)
    weights=np.array([C@(math.tanh(b)**np.arange(85)) for b in beta])
    probs=weights/weights.sum(axis=1)[:,None]
    for h,col,ls in [(0,DARK,'-'),(1,TEAL,'-'),(2,ORANGE,'--'),(3,PURPLE,'-')]:
        ax[1].plot(beta,probs[:,h],lw=2,color=col,ls=ls,label=f'{h:02b}')
    ax[1].set(title='Two homology coordinates',xlabel='Inverse temperature β',ylabel='Even-subgraph class probability',xlim=(0,1.4),ylim=(0,1.02))
    ax[1].legend(frameon=False,ncol=2,fontsize=9,title='Retained parity bits',title_fontsize=9)
    ax[1].grid(axis='y',alpha=.15)
    kappas=np.linspace(0,.4,161)
    bilayer_table=[]
    for b,col in [(.2,PURPLE),(.64,ORANGE),(1.,TEAL)]:
        rows=[bilayer_observables(layer,b,k) for k in kappas]
        bilayer_table.extend(rows)
        ax[2].plot(kappas,[r['overlap_per_site'] for r in rows],lw=2.2,color=col,label=f'β = {b:.2f}')
    ax[2].set(title='Coupling two 16-spin layers',xlabel='Inter-layer coupling κ',ylabel='Mean overlap ⟨Σi σiτi⟩ / 16',xlim=(0,.4),ylim=(0,1.02))
    ax[2].legend(frameon=False,fontsize=9,loc='lower right');ax[2].grid(axis='y',alpha=.15)
    fig.text(.065,.13,'26 flat twists checked under gauge changes.      All 2²⁹ even subgraphs checked independently.      At κ = 0: Zdouble = Zsingle² exactly.',fontsize=10)
    fig.text(.065,.05,'The tree/cotree basis is a recorded coordinate choice; only two homology bits are retained in the middle panel.\nThe bilayer interaction κ is introduced explicitly. All curves evaluate exact finite counting tables numerically.',fontsize=9.5,color='#53616e')
    save(fig,'MTFT_Twists_and_Layers')
    with (ROOT/'bilayer_response.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(bilayer_table[0]));w.writeheader();w.writerows(bilayer_table)


def defect_figure(pos,edges):
    r=read('defect_census.json');w=max(r['rows'],key=lambda row:row['minimum_rearrangement'])
    ref=np.array(r['baseline_ground']['spin_bits']);new=np.array(w['witness_bits']);changed=ref!=new
    counts=Counter(row['minimum_rearrangement'] for row in r['rows']);dist=sorted(counts)
    fig,ax=plt.subplots(1,2,figsize=(12.9,6.6),gridspec_kw={'width_ratios':[1.2,1]})
    fig.subplots_adjust(left=.075,right=.965,top=.79,bottom=.22,wspace=.21)
    header(fig,'One reversed bond can force nine spins to move','N = 143 antiferromagnet  ·  exhaustive census of all 84 single-bond reversals')
    draw_graph(ax[0],pos,edges,[ORANGE if c else '#8da6ae' for c in changed],np.where(changed,105,45),w['edge'])
    for v in w['endpoints']:
        ax[0].annotate(str(v),pos[v],xytext=(5,7),textcoords='offset points',fontweight='bold',fontsize=10,color=PURPLE)
    ax[0].set_title(f'Edge {w["edge"]}: orange nodes must change in this closest witness',loc='left',fontsize=11,pad=11)
    bars=ax[1].bar(range(len(dist)),[counts[d] for d in dist],color=[TEAL if d<9 else ORANGE for d in dist],width=.68)
    for bar,d in zip(bars,dist):
        ax[1].text(bar.get_x()+bar.get_width()/2,bar.get_height()+.8,str(counts[d]),ha='center',fontsize=11)
    ax[1].set_xticks(range(len(dist)),dist)
    ax[1].set(xlabel='Minimum Hamming distance to the original ground state',ylabel='Number of bond reversals',ylim=(0,43))
    ax[1].set_title('Closest new ground state, allowing global spin flip',loc='left',fontsize=11,pad=11)
    ax[1].grid(axis='y',alpha=.15);ax[1].set_axisbelow(True)
    fig.text(.075,.12,'For edge 4, the two new ground states lie at distances 9 and 47. Both have energy −66; the original ground energy is also −66.',fontsize=10.5)
    fig.text(.075,.055,'Exact ground counts and min-plus optimization agree. A separate joint-count calculation certifies the distances 9 and 47.\nPurple marks the reversed bond; orange marks changes in the displayed witness. Positions are a drawing convention.',fontsize=9.5,color='#53616e')
    save(fig,'MTFT_One_Bond_Defect')
    with (ROOT/'defect_census.csv').open('w',newline='') as f:
        fields=['edge','self_loop','was_frustrated','min_frustrated','ground_energy','ground_degeneracy','minimum_rearrangement']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        writer.writerows({k:row[k] for k in fields} for row in r['rows'])


def zeros_figure():
    f=read('fisher_zeros.json');lee=read('lee_yang.json')['rows']
    roots=np.array([complex(float(z['real']),float(z['imag'])) for z in f['Q45_roots']])
    fig,ax=plt.subplots(1,3,figsize=(15.3,6.2),gridspec_kw={'width_ratios':[1,1,1.05]})
    fig.subplots_adjust(left=.065,right=.975,top=.76,bottom=.27,wspace=.37)
    header(fig,'Where the finite partition functions vanish','N = 143  ·  complex temperature x = exp(−2β) and complex field y = exp(2η)',left=.065)
    for a in ax[:2]:
        a.axhline(0,color=GRAY,lw=.8);a.axvline(0,color=GRAY,lw=.8)
        a.scatter(roots.real,roots.imag,s=26,color=TEAL,zorder=3)
        a.set(xlabel='Re x',ylabel='Im x');a.set_aspect('equal',adjustable='datalim')
    ax[0].scatter([-1,0,0],[0,1,-1],marker='D',s=[52,32,32],color=ORANGE,zorder=4)
    ax[0].annotate('−1 (×28)',(-1,0),xytext=(-6.4,-1.6),arrowprops={'arrowstyle':'-','color':ORANGE},fontsize=9)
    ax[0].set_title('Temperature zeros: full root set',loc='left',fontsize=11,pad=12)
    ax[0].add_patch(Rectangle((-.05,-.46),.5,.92,fill=False,ec=PURPLE,lw=1))
    ax[1].set_aspect('auto');ax[1].set_xlim(-.05,.45);ax[1].set_ylim(-.46,.46)
    ax[1].plot([0,.45],[0,0],color=DARK,lw=2,zorder=2)
    ax[1].text(.2,.025,'Real temperatures',ha='center',fontsize=9,color=DARK)
    ax[1].set_title('Temperature zeros: positive-axis detail',loc='left',fontsize=11,pad=12)
    ax[1].annotate('0.27188 + 0.10008i',(.271884,.100081),xytext=(.06,.32),arrowprops={'arrowstyle':'-','color':TEAL},fontsize=9)
    theta=np.linspace(0,2*np.pi,500)
    ax[2].plot(np.cos(theta),np.sin(theta),color=GRAY,lw=1.2)
    for row,col,marker in zip(lee,[PURPLE,ORANGE,TEAL],['o','x','+']):
        rr=row['positive_imaginary_roots'];re=np.array([z['re_y'] for z in rr]);im=np.array([z['im_y_positive'] for z in rr])
        ax[2].scatter(np.r_[re,re],np.r_[im,-im],color=col,marker=marker,s=23,label=f'a/d = {row["satisfied_weight"]}')
    ax[2].set(xlabel='Re y',ylabel='Im y',xlim=(-1.14,1.14),ylim=(-1.14,1.14))
    ax[2].set_aspect('equal');ax[2].set_title('Field zeros: certified unit circle',loc='left',fontsize=11,pad=12)
    ax[2].legend(frameon=False,fontsize=8.5,loc='center')
    fig.text(.065,.15,'Exact factorization:  D(x) = 2(1 + x)²⁸(1 + x²)Q₄₅(x).       Closest field angles: 0.41934, 0.06761, 0.05684 radians for a/d = 2, 4, 8.',fontsize=10)
    fig.text(.065,.06,'Temperature: 45 Q roots computed at 40 and 65 digits; maximum normalized residual ≈ 4.6 × 10⁻⁶⁷. Orange diamonds are exact roots.\nField: exact reciprocal reduction plus rational root isolation certifies all 56 zeros on |y| = 1 at each displayed coupling. Coordinates are numerical.',fontsize=9.5,color='#53616e')
    save(fig,'MTFT_Complex_Zeros')


def sampling_figure(pos,edges):
    joint=read('joint_143.json');D=np.array(joint['counts'],dtype=float);n=56;E=84
    samples=np.load(ROOT/'equilibrium_samples.npz')
    meta=read('sampling.json')['rows']
    fig,ax=plt.subplots(2,4,figsize=(14.4,7.4),gridspec_kw={'height_ratios':[1.25,1]})
    fig.subplots_adjust(left=.065,right=.975,top=.81,bottom=.23,wspace=.27,hspace=.16)
    header(fig,'Independent equilibrium draws','N = 143  ·  exact integer conditional probabilities  ·  20,000 configurations at each coupling',left=.065)
    m=np.arange(n+1);M=(2*m-n)/n
    for j,r in enumerate(meta):
        q=r['a'];bits=samples[f'q{q}_spin_bits'][0]
        draw_graph(ax[0,j],pos,edges,[ORANGE if b else TEAL for b in bits],26)
        ax[0,j].set_title(f'a/d = {q}   ·   β = {r["beta"]:.4f}',fontsize=11,pad=6)
        logedge=r['beta']*(E-2*np.arange(E+1));w=D*np.exp(logedge-logedge.max())[:,None]
        exact=w.sum(axis=0)/w.sum();observed=samples[f'q{q}_magnetization_histogram']/r['sample_count']
        ax[1,j].plot(M,exact,lw=1.8,color=TEAL,label='Exact distribution')
        ax[1,j].scatter(M,observed,color=ORANGE,s=13,marker='x',label='20,000 draws',zorder=3)
        ax[1,j].set(xlabel='Magnetization per spin',xlim=(-1.06,1.06),ylim=(0,max(exact.max(),observed.max())*1.16))
        ax[1,j].set_xticks([-1,0,1]);ax[1,j].grid(axis='y',alpha=.15)
        if j==0:ax[1,j].set_ylabel('Probability mass')
    handles,labels=ax[1,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.52,.135),ncol=2,frameon=False,fontsize=10)
    fig.text(.065,.09,'The top row shows one recorded draw per coupling; colors denote opposite spins. Sampling uses exact integer weights and recorded pseudorandom seeds.',fontsize=10)
    fig.text(.065,.045,'No equilibration burn-in is used. These configurations are samples from finite Gibbs distributions, not a physical time trajectory.',fontsize=9.5,color='#53616e')
    save(fig,'MTFT_Equilibrium_Samples')
    animate_samples(pos,edges,samples,meta)


def animate_samples(pos,edges,samples,meta):
    from matplotlib import font_manager
    fontfile=font_manager.findfont('DejaVu Sans')
    titlefont=ImageFont.truetype(fontfile,27);font=ImageFont.truetype(fontfile,17);small=ImageFont.truetype(fontfile,14)
    frames=[]
    for frame in range(32):
        im=Image.new('RGB',(1200,395),'white');draw=ImageDraw.Draw(im)
        draw.text((35,18),'Independent equilibrium draws  ·  N = 143',font=titlefont,fill=DARK)
        draw.text((1065,25),f'Draw {frame+1:02d}',font=font,fill=DARK)
        for j,r in enumerate(meta):
            q=r['a'];bits=samples[f'q{q}_spin_bits'][frame]
            draw.text((j*292+45,68),f'a/d = {q}   β = {r["beta"]:.4f}',font=font,fill=DARK)
            xy=pos.copy();xy[:,0]=j*292+151+112*pos[:,0];xy[:,1]=224-107*pos[:,1]
            for u,v in edges:
                if u==v:
                    x,y=xy[u];draw.ellipse((x-8,y-8,x+5,y+5),outline=GRAY,width=1)
                else:draw.line([tuple(xy[u]),tuple(xy[v])],fill=GRAY,width=1)
            for (x,y),b in zip(xy,bits):
                draw.ellipse((x-4,y-4,x+4,y+4),fill=ORANGE if b else TEAL,outline='white',width=1)
        draw.text((35,365),'Each frame is a fresh equilibrium draw; playback speed has no physical meaning.',font=small,fill='#53616e')
        frames.append(im)
    frames[0].save(ROOT/'MTFT_Equilibrium_Draws.gif',save_all=True,append_images=frames[1:],duration=600,loop=0,optimize=True)


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.labelcolor':DARK,'text.color':DARK,'xtick.color':'#53616e','ytick.color':'#53616e','svg.fonttype':'none'})
    joint=read('joint_143.json');edges=joint['edges'];pos=graph_layout(56,edges)
    np.save(ROOT/'graph143_drawing_positions.npy',pos)
    geometry_figure();topology_layers_figure();defect_figure(pos,edges);zeros_figure();sampling_figure(pos,edges)
    with (ROOT/'joint_counts_143.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['disagreeing_edges_k','up_spins_m','S','M','exact_count'])
        w.writerows((k,m,84-2*k,2*m-56,d) for k,row in enumerate(joint['counts']) for m,d in enumerate(row) if d)
    print('Wrote five PNG/SVG figure pairs, equilibrium GIF, drawing positions, and three CSV tables.')


if __name__=='__main__':main()
