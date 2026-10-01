#!/usr/bin/env python3
"""Standalone scientific figure of the conditional Matter-03 results."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from mtft.matter.charged_modes import certified_paired_charged_report
    from mtft.matter.higgs_vacuum import paired_higgs_witness
    from mtft.matter.collinear import yukawa_witness
    from mtft.matter.kinetics import canonical_yukawa, mass_singular_values
    report=certified_paired_charged_report()
    witness=paired_higgs_witness()
    # The figure is tied to certified endpoint invariants, not a handpicked scan.
    assert witness.product_exceptional_valuation == 7
    assert report['a_map_rank'] == 41
    blue='#2468a2'; orange='#dc8747'; green='#2c876c'; ink='#19304a'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(13.5,9.5),layout='constrained')
    fig.suptitle('Matter-03: a Higgs branch exists; family selection remains unresolved',
                 fontsize=19,fontweight='bold',color=ink)
    ax=axes[0,0]; ax.set(xlim=(0,1),ylim=(0,1)); ax.axis('off')
    ax.set_title('A. Classical paired branch',loc='left',fontweight='bold')
    for x,label in ((.22,r'$\mathcal{L}$'),(.78,r'$\mathcal{L}^{-1}$')):
        ax.text(x,.72,label,ha='center',va='center',fontsize=22,
                bbox={'boxstyle':'round,pad=.6','facecolor':'#edf3f9','edgecolor':blue})
    for p1,p2,rad,label,y in (((.68,.76),(.32,.76),.25,r'$H_d$',.96),
                            ((.32,.64),(.68,.64),.25,r'$H_u$',.45)):
        ax.add_patch(FancyArrowPatch(p1,p2,connectionstyle=f'arc3,rad={rad}',
                     arrowstyle='-|>',mutation_scale=17,linewidth=2,color=blue))
        ax.text(.5,y,label,ha='center',fontsize=16,color=blue)
    ax.text(.5,.27,'Exceptional valuation 7 ⇒ irreducible spectral polynomial',
            ha='center',fontsize=11,color=ink)
    ax.text(.5,.16,'Stable active pair + wall spectators ⇒ classical metric exists',
            ha='center',fontsize=10,color=ink)
    ax.text(.5,.045,'Common nonzero scaling is free; no selected Higgs scale',
            ha='center',fontsize=11,fontweight='bold',color=orange)
    ax=axes[0,1]
    ax.set_title('B. Charged degree-one classes per Q sign',loc='left',fontweight='bold')
    labels=['Split background','Certified pair','Any aligned pair\n(lower bound)']
    vals=[279,237,160]; colors=[blue,green,'#aebcc8']
    bars=ax.barh(labels,vals,color=colors,height=.55)
    for bar,value in zip(bars,vals):
        ax.text(value+5,bar.get_y()+bar.get_height()/2,str(value),va='center',fontweight='bold')
    ax.set_xlim(0,325); ax.invert_yaxis(); ax.set_xlabel('Holomorphic classes, not observed-particle counts')
    ax.text(.02,.02,'Specified pair: 279 − 1 gauge quotient − 41 aa rank = 237',
            transform=ax.transAxes,fontsize=10,color=ink)
    ax=axes[1,0]
    ax.set_title('C. Six charged aa blocks',loc='left',fontweight='bold')
    ranks=[2,2,5,3,17,12]; caps=[14,14,18,15,32,25]
    positions=np.arange(6)
    ax.bar(positions,caps,color='#d7e0e9',label='Maximum from block size')
    ax.bar(positions,ranks,color=green,label='Rank of certified pair')
    for pos,val in zip(positions,ranks): ax.text(pos,val+.6,str(val),ha='center',fontsize=10)
    ax.set_xticks(positions,['L−\ndown','L−\nup','A1\ndown','A1\nup','A2\ndown','A2\nup'])
    ax.set_ylabel('Holomorphic rank'); ax.set_ylim(0,40)
    ax.legend(frameon=False,fontsize=10,loc='upper left')
    ax.text(.98,.98,'Actual total: 41\nCapacity total: 118',transform=ax.transAxes,
            ha='right',va='top',fontsize=10,color=ink)
    ax=axes[1,1]
    ax.set_title('D. INPUT metrics change coupling ratios',loc='left',fontweight='bold')
    y=np.array(yukawa_witness().residue_matrix,dtype=complex)[:,:,None]
    simple=canonical_yukawa(y,np.eye(3),np.eye(3),[[1]],kappa=1)
    changed=canonical_yukawa(y,np.diag([1,4,25]),np.diag([1,9,16]),[[1]],kappa=1)
    s1=mass_singular_values(simple[:,:,0]); s2=mass_singular_values(changed[:,:,0])
    assert np.allclose(s2,[1,1/8,1/15])
    xs=np.arange(3)
    ax.bar(xs-.18,s1,width=.36,color=blue,label='Identity Grams')
    ax.bar(xs+.18,s2,width=.36,color=orange,label='Supplied diagonal Grams')
    for x,value,label in zip(xs+.18,s2,['1','1/8','1/15']): ax.text(x,value+.03,label,ha='center',fontsize=10)
    ax.set_xticks(xs,['σ1','σ2','σ3']); ax.set_ylim(0,1.5)
    ax.set_ylabel('Normalized witness coupling; arbitrary units')
    ax.legend(frameon=False,fontsize=10)
    ax.text(.5,.77,'Rank 3 in both examples\nNo selected families or mass predictions',
            transform=ax.transAxes,ha='center',va='top',fontsize=10,color=ink)
    fig.supxlabel('A–C: exact conditional geometry of the fixed B3 bulk parent.  D: supplied metric example.\n'
                  'Actual Hermitian metric, Higgs scale, light-family lifting and physical photon remain open.',
                  fontsize=10,color=ink)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=180,facecolor='white')
    plt.close(fig)
    print('Matter-03 scientific figure:',args.output)


if __name__=='__main__': main()
