#!/usr/bin/env python3
"""Plot exact Matter-02 bookkeeping, jets and a holomorphic residue minor."""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
                         'axes.titlesize':14,'axes.labelsize':11})
    fig, axs=plt.subplots(2,2,figsize=(14,11.5),facecolor='#f7f9fc')
    fig.subplots_adjust(left=.075,right=.965,top=.865,bottom=.18,wspace=.29,hspace=.80)
    fig.suptitle('Matter-02  |  Exact modes and a holomorphic interaction',
                 x=.075,y=.965,ha='left',fontsize=21,fontweight='bold',color='#17283e')
    fig.text(.075,.92,'Specified cusp blow-up of X₀(143) × 143a1 • split background Φ = 0',
             fontsize=12,color='#42536a')
    qs=np.arange(-5,6)
    h2=np.array([13 if q==0 else max(13-q,0) if q>0 else max(14+2*q,0) for q in qs])
    h1=h2+(9*qs**2+3*qs)//2
    h1[qs==0]=14
    phi=np.array([13 if q==0 else max(13+q,0) if q<0 else max(14-2*q,0) for q in qs])
    ax=axs[0,0]
    ax.plot(qs,h1,'o-',color='#2575b8',label='h¹(ℒᑫ)')
    ax.plot(qs,h2,'s--',color='#8961b2',label='h²(ℒᑫ)')
    ax.plot(qs,phi,'^:',color='#d18b18',label='h⁰(Kℒᑫ)')
    ax.set(title='A. Actual dimensions exceed net indices',xlabel='Hom exponent q',ylabel='Complex dimension')
    ax.set_xticks(qs);ax.set_ylim(-3,143);ax.grid(alpha=.18);ax.legend(frameon=False,fontsize=10)
    ax=axs[0,1]
    selected={(0,1),(0,2),(1,1)}
    for z in range(4):
        for w in range(4):
            if z+w>=3:
                ax.scatter(z,w,s=330,color='#e5eaf1',marker='s')
                continue
            if (z,w)==(1,0): color='#a7afb9';text='z\nremoved'
            elif (z,w) in selected:
                color='#e7b545';text={(0,1):'w',(0,2):'w²',(1,1):'zw'}[(z,w)]
            else: color='#8ec6df';text='1' if (z,w)==(0,0) else 'z²'
            ax.scatter(z,w,s=1030,color=color,marker='s',edgecolor='white',linewidth=2)
            ax.text(z,w,text,ha='center',va='center',fontsize=10,color='#15293e')
    ax.set(title='B. Fat jets at the blow-up point',xlabel='Power of z',ylabel='Power of w',
           xlim=(-.45,3.45),ylim=(-.45,3.45))
    ax.set_xticks(range(4));ax.set_yticks(range(4));ax.set_aspect('equal')
    ax.text(.5,-.29,'Total degree < 3: six jets, one global image.\nGold: the three classes used in the rank witness.',
            transform=ax.transAxes,ha='center',fontsize=10,color='#42536a')
    ax=axs[1,0]
    matrix=np.array([[1,0,0],[0,0,1],[0,1,0]])
    ax.imshow(matrix,cmap='Blues',vmin=0,vmax=1)
    for i in range(3):
        for j in range(3):
            ax.text(j,i,str(matrix[i,j]),ha='center',va='center',fontsize=23,
                    color='white' if matrix[i,j] else '#56687c',fontweight='bold')
    ax.set_xticks(range(3),['z⁻⁶w⁻²','z⁻⁷w⁻²','z⁻⁶w⁻³'],fontsize=10)
    ax.set_yticks(range(3),['w','w²','zw'])
    ax.set_title('C. An exact residue minor: det = −1',pad=16)
    ax.text(.5,-.28,'A slice of the 18 × 44 × 10 holomorphic tensor.\nRank ≥ 3; neither a Higgs vev nor physical masses.',
            transform=ax.transAxes,ha='center',fontsize=10,color='#42536a')
    ax=axs[1,1];ax.axis('off')
    ax.set_title('D. Physical gates still to resolve',loc='left',pad=16)
    rows=[['Holomorphic rank lower bound','≥ 3'],['Linear Y–flux projection','0'],
          ['Quadratic Y–flux projection','54'],['Three light families / lifting','OPEN'],
          ['Normalized Yukawa and Higgs vacuum','OPEN'],['Parent axion/photon mass kernel','OPEN']]
    table=ax.table(cellText=rows,colWidths=[.76,.24],cellLoc='left',bbox=[0,.17,1,.73])
    table.auto_set_font_size(False);table.set_fontsize(10)
    for (i,j),cell in table.get_celld().items():
        cell.set_edgecolor('#dce3ec');cell.set_facecolor('#edf3fa' if i<3 else '#fff4e1')
        if j==1: cell.set_text_props(weight='bold')
    ax.text(0,.025,'Flux projections are exact traces. They are not\na derived Stückelberg action or a photon mass.',
            transform=ax.transAxes,fontsize=10,color='#42536a')
    fig.text(.075,.032,'Exact conditional on the archived basis, specified holomorphic representative and chosen bulk parent.\n'
             'Local jet coordinates are not a spacetime image. No observed particle masses are derived.',
             fontsize=10,color='#42536a')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=160,facecolor=fig.get_facecolor())
    print('Scientific figure:',args.output)


if __name__=='__main__': main()
