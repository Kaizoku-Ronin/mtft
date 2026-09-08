"""Generate the report and two standalone scientific figures from saved results."""
from collections import Counter
from fractions import Fraction
import hashlib
import platform
import numpy as np
import scipy
import sympy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from common import ROOT, BLOCKS, DIMS, CLASSES, CLASS_LABELS, read, dump

INK='#243447'
MUTED='#526577'
COLORS=['#466a9f','#16847d','#b36c30','#87559b','#bf4e61']

def pattern(p, mathtext=False):
    if mathtext:
        parts=[str(d) if n==1 else f'{d}^{{{n}}}' for d,n in sorted(Counter(p).items())]
        return '$'+r'\,'.join(parts)+'$'
    return '+'.join(map(str,p))

def save(fig,stem):
    fig.savefig(ROOT/f'{stem}.png',dpi=180,facecolor='white')
    fig.savefig(ROOT/f'{stem}.svg',facecolor='white')
    plt.close(fig)

def atlas_figure(rows,types,summary):
    q4=[c['pattern'] for c in types['marginals']['q4']]
    q6=[c['pattern'] for c in types['marginals']['q6']]
    counts=np.array([c['count'] for c in types['joint']]).reshape(5,11)
    densities=np.array([100*float(Fraction(c['density'])) for c in types['joint']]).reshape(5,11)
    fig=plt.figure(figsize=(14,11))
    gs=fig.add_gridspec(2,1,height_ratios=[1,1.18],left=.08,right=.97,top=.86,bottom=.14,hspace=.4)
    fig.suptitle('Two Hecke fields, 55 infinite prime traces',x=.08,y=.975,ha='left',fontsize=21,weight='bold',color=INK)
    fig.text(.08,.935,f"1,248 jointly unramified primes below 10,201 · {summary['joint_types_seen']} types seen · four types absent in this finite sample",fontsize=11.5,color=INK)
    ax=fig.add_subplot(gs[0]); ax.imshow(counts+1,cmap='Blues',norm=LogNorm(vmin=1,vmax=counts.max()+1),aspect='auto')
    ax.set_xticks(range(11),[pattern(p,True) for p in q6],fontsize=10)
    ax.set_yticks(range(5),[pattern(p,True) for p in q4],fontsize=11)
    ax.set_xlabel('Sextic factor degrees',labelpad=8); ax.set_ylabel('Quartic factor degrees',labelpad=10)
    ax.set_title('Each cell: observed prime count / exact limiting percentage',loc='left',fontsize=12,pad=13,color=INK)
    for i in range(5):
        for j in range(11):
            color='white' if counts[i,j]>=25 else INK
            ax.text(j,i,f'{counts[i,j]}\n{densities[i,j]:.3g}%',ha='center',va='center',fontsize=10,color=color,linespacing=1.6)
    ax.set_xticks(np.arange(-.5,11,1),minor=True); ax.set_yticks(np.arange(-.5,5,1),minor=True)
    ax.grid(which='minor',color='white',linewidth=1.5); ax.tick_params(which='minor',bottom=False,left=False)
    ax=fig.add_subplot(gs[1])
    for y in range(11): ax.axhline(y,color='#dfe5eb',lw=.8,zorder=0)
    for ci,part in enumerate(q4):
        chosen=[r for r in rows if r['joint_unramified'] and r['q4']['pattern']==part]
        xs=[r['square_n']+r['square_position_numerator']/r['square_position_denominator'] for r in chosen]
        ys=[q6.index(r['q6']['pattern']) for r in chosen]
        sizes=[29 if r['irregularity_index']>0 else 10 for r in chosen]
        ax.scatter(xs,ys,s=sizes,c=COLORS[ci],alpha=.78,edgecolors='white',linewidths=.18,zorder=2)
    for p,dx,dy in [(1093,3.5,.45),(3511,3.5,.45)]:
        r=next(r for r in rows if r['p']==p)
        x=r['square_n']+r['square_position_numerator']/r['square_position_denominator']; y=q6.index(r['q6']['pattern'])
        ax.scatter([x],[y],marker='*',s=125,color=INK,edgecolors='white',linewidths=.5,zorder=4)
        ax.annotate(str(p),(x,y),(x+dx,y+dy),fontsize=9,color=INK,arrowprops={'arrowstyle':'-','color':INK,'lw':.8},
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':1},zorder=5)
    ax.set(xlim=(1,101),ylim=(-.7,10.9),yticks=range(11),yticklabels=[pattern(p,True) for p in q6],
           ylabel='Sextic factor degrees',xlabel='Square-interval position: each unit is one interval n² < p < (n+1)²')
    ax.set_title('Prime traces: color = quartic type; larger dots = irregular primes',loc='left',fontsize=12,pad=12,color=INK)
    handles=[Line2D([0],[0],marker='o',ls='',color=c,label=pattern(p,True),markersize=6) for c,p in zip(COLORS,q4)]
    fig.legend(handles=handles,title='Quartic factor degrees',loc='lower left',bbox_to_anchor=(.073,.045),ncol=5,frameon=False,fontsize=10,title_fontsize=10)
    fig.text(.49,.073,'Stars: base-2 Wieferich primes 1093 and 3511.\nThe irregular overlay carries no new infinitude claim.',fontsize=9.5,color=MUTED,linespacing=1.5)
    fig.text(.08,.019,'Excluded from the joint atlas: 5, 7, 19, 103. Level primes 11 and 13 are unramified in these two coefficient fields and remain included.',fontsize=9,color=MUTED)
    save(fig,'Hecke_Prime_Atlas')

def stability_figure(summary,certificate):
    fig=plt.figure(figsize=(14,10.6))
    gs=fig.add_gridspec(2,6,height_ratios=[1,1.1],left=.074,right=.969,top=.845,bottom=.155,hspace=.52,wspace=.9)
    fig.suptitle('Which parts of the Hecke structure survive?',x=.074,y=.975,ha='left',fontsize=21,weight='bold',color=INK)
    fig.text(.074,.934,'400 conservative flows · 16 paired seeds · equal Frobenius perturbation size · dimensionless time 0–80',fontsize=11.5,color=INK)
    ax=fig.add_subplot(gs[0,:3]); ranks=np.array(certificate['AL_intersection_real_dimensions'])
    ax.imshow(ranks,cmap='Blues',vmin=0,vmax=12,aspect='auto')
    ax.set_xticks(range(4),['(+,+)','(+,−)','(−,+)','(−,−)']); ax.set_yticks(range(4),['Elliptic (2)','Oldspace (4)','Quartic (8)','Sextic (12)'])
    ax.set_xlabel('(W11, W13) symmetry character',labelpad=9)
    ax.set_title('Exact dimensions: oldspace and quartic overlap',loc='left',fontsize=11.8,pad=12,color=INK)
    for i in range(4):
        for j in range(4): ax.text(j,i,str(ranks[i,j]),ha='center',va='center',fontsize=14,color='white' if ranks[i,j]>=6 else INK)
    ax.set_xticks(np.arange(-.5,4,1),minor=True); ax.set_yticks(np.arange(-.5,4,1),minor=True)
    ax.grid(which='minor',color='white',lw=2); ax.tick_params(which='minor',bottom=False,left=False)
    ax=fig.add_subplot(gs[0,3:]); eps=summary['epsilon_values']
    for ci in [0,1,2,3,4]:
        c=summary['class_results'][ci]
        values=np.maximum(0,c['mean_peak_aggregate_leakage_by_epsilon'])*100
        ax.plot(eps,values,color=COLORS[ci],marker=['s','o','o','D','o'][ci],ms=4,
                ls='--' if ci<2 else '-',lw=1.7,label=CLASS_LABELS[ci],alpha=.9)
    ax.set(xlabel='Perturbation strength ε',ylabel='Mean peak total leakage (%)',ylim=(-1.5,65),xlim=(-.006,.406))
    ax.set_title('Less symmetry permits more exchange in this model',loc='left',fontsize=11.8,pad=12,color=INK)
    ax.legend(loc='upper left',fontsize=8.7,frameon=False); ax.grid(alpha=.18)
    for panel,ci in enumerate([2,3,4]):
        ax=fig.add_subplot(gs[1,2*panel:2*panel+2])
        values=100*np.array(summary['class_results'][ci]['time_and_seed_mean_transfer'])
        ax.imshow(values,cmap='Blues',vmin=0,vmax=100)
        ax.set_xticks(range(4),['Ell','Old','Q4','Q6']); ax.set_yticks(range(4),['Ell','Old','Q4','Q6'])
        ax.set_xlabel('Source sector'); ax.set_ylabel('Destination sector' if panel==0 else '')
        ax.set_title(CLASS_LABELS[ci]+'\nMean transfer at ε = 0.4 (%)',fontsize=11,pad=14,color=INK)
        for i in range(4):
            for j in range(4):
                text='0' if values[i,j]<1e-7 else f'{values[i,j]:.1f}'
                ax.text(j,i,text,ha='center',va='center',fontsize=10.4,color='white' if values[i,j]>=55 else INK)
        ax.set_xticks(np.arange(-.5,4,1),minor=True); ax.set_yticks(np.arange(-.5,4,1),minor=True)
        ax.grid(which='minor',color='white',lw=1.3); ax.tick_params(which='minor',bottom=False,left=False)
    fig.text(.074,.07,'The middle perturbation has zero diagonal block matrices, yet still mixes oldspace and quartic modes.\nExact AL selection rule: a uniform oldspace source can lose at most 50%; elliptic and sextic sectors remain protected.',fontsize=10,color=MUTED,linespacing=1.6)
    fig.text(.074,.022,'Top-right: average over seeds of each flow’s maximum on the fixed time grid. Bottom: averages over both time and seeds; every column sums to 100%.',fontsize=9,color=MUTED)
    save(fig,'Hecke_Sector_Stability')

def report(rows,types,atlas,stability,fields,sectors,validation):
    q4table='\n'.join(f"| {pattern(c['pattern'])} | {c['density']} | {c['count']} |" for c in types['marginals']['q4'])
    q6table='\n'.join(f"| {pattern(c['pattern'])} | {c['density']} | {c['count']} |" for c in types['marginals']['q6'])
    unseen='\n'.join(f"| {pattern(c['q4'])} | {pattern(c['q6'])} | {c['density']} |" for c in types['joint'] if c['count']==0)
    results=[]
    for label,c in zip(CLASS_LABELS,stability['class_results']):
        peak=c['mean_peak_aggregate_leakage']*100; old=c['mean_peak_leakage_by_source'][1]*100
        results.append(f"| {label} | {'0 (roundoff)' if peak<1e-9 else f'{peak:.2f}%'} | {'0 (roundoff)' if old<1e-9 else f'{old:.2f}%'} | {c['largest_frequency_shift']:.4f} |")
    resulttable='\n'.join(results)
    complete6=[r['p'] for r in rows if not r['q6']['ramified'] and r['q6']['pattern']==[1]*6]
    text=r'''# Hecke fields and sector stability

MTFT 0.26.0 companion experiment · 6 September 2026 · v0.1.0

The experiments give the “tissue” idea two precise forms: an algebraic web of **55 provably infinite prime-splitting classes**, and an exact map of which Hecke sectors can exchange amplitude when some symmetries are preserved. The quartic and sextic splitting fields are linearly disjoint. On X0(143), the oldspace and quartic sector nevertheless share an Atkin–Lehner character, so those involutions alone do not separate all Hecke sectors. These are different relationships: field disjointness concerns arithmetic extensions; sector mixing concerns operators on homology.

## 1. The prime-splitting atlas

**EXACT arithmetic with an accompanying mathematical argument.** The T2 factors are

\[
g_4(x)=x^4-3x^3-x^2+5x+1,
\qquad
h_6(x)=x^6-10x^4+2x^3+24x^2-7x-12.
\]

| Coefficient field | Galois group of its splitting field | Discriminant |
|---|---|---|
| Quartic | S4, order 24 | 1,957 = 19 × 103 |
| Sextic | S6, order 720 | 194,616,205 = 5 × 7 × 5,560,463 |

Both discriminants are squarefree, so they are also the number-field discriminants of the corresponding monogenic fields. The coefficient fields themselves have degrees 4 and 6; the orders 24 and 720 refer to their Galois closures.

The groups were computed with SymPy and independently established from modular factorizations:

| Field | Irreducible reduction | A cycle of length n−1 | A transposition |
|---|---|---|---|
| Quartic | p = 2: (4) | p = 3: (1,3) | p = 43: (1,1,2) |
| Sextic | p = 19: (6) | p = 3: (1,5) | p = 307: (1,1,1,1,2) |

An irreducible reduction makes the group transitive. The long cycle makes a point stabilizer transitive on all the other roots; conjugating a transposition then supplies every transposition. This proves Sn in each case. The reduction-to-cycle theorem and this group criterion are [Milne, Theorem 4.28 and Lemma 4.31](https://www.jmilne.org/math/CourseNotes/FT.pdf). All factor coefficients and independent finite-field irreducibility checks are retained in the certificate.

The ramification supports are disjoint. An intersection of the two Galois closures would therefore be unramified at every finite prime over Q. The Minkowski bound forces that intersection to be Q; see [Sutherland, Corollary 14.27](https://math.mit.edu/classes/18.785/2021fa/LectureNotes14.pdf). Thus the compositum has group S4 × S6 and degree 17,280. An alternative check is that the only possible nontrivial common quotient of S4 and S6 is C2, but their discriminant quadratic fields differ. In particular, the coefficient fields have no shared nontrivial subfield.

For a cycle partition with m_j cycles of length j, the proportion of permutations is

\[
\delta(\lambda)=\frac{1}{\prod_j j^{m_j}m_j!}.
\]

The product group gives joint densities δ(λ4)δ(λ6). Applying the [Chebotarev density theorem](https://math.mit.edu/classes/18.785/2015fa/LectureNotes25.pdf) makes every one of the 5 × 11 = **55** joint types infinite. This is a deduction from exact field structure and a theorem, not an extrapolation from the finite counts.

### Finite census

We factored both polynomials at all **1,252 primes below 10,201**. The joint unramified atlas contains **1,248** primes after excluding **5, 7, 19, 103**. The sextic ramification prime 5,560,463 lies outside the range. The level primes **11 and 13 are retained**: bad reduction of the modular curve and ramification of these coefficient fields are distinct properties.

A partition such as 1+1+2 means two linear factors and one irreducible quadratic. Each marginal below has 1,250 unramified primes. Densities are asymptotic targets; multiplying one by the finite sample size gives a reference count, not a proven finite-sample expectation under independent sampling.

Quartic:

| Factor degrees | Exact density | Count |
|---|---|---:|
__Q4_TABLE__

Sextic:

| Factor degrees | Exact density | Count |
|---|---|---:|
__Q6_TABLE__

The sextic completely splits at just **__COMPLETE6__** in this sample. No sampled prime completely splits in both fields; that joint class still has positive density **1/17,280**. Exactly **51 of the 55 joint classes** appear. The four missing classes are:

| Quartic type | Sextic type | Proved positive density |
|---|---|---|
__UNSEEN_TABLE__

This supplies the requested contrast between sparse observations and known infinitude. The atlas preserves the square-interval position and earlier irregularity labels of each prime. For example, **1093** is regular with types (1,1,2) and (1,2,3); **3511** has irregularity index 2 with types (1,1,2) and (1,1,2,2). Both are base-2 Wieferich primes. These two marked examples carry no statistical inference.

The “proved infinite” label belongs to each splitting class. This experiment does **not** establish the infinitude of its intersection with irregular primes, statistical independence of irregularity, or a prime in every square interval. A positive global density gives no such interval guarantee. No post-hoc dependence test was fitted to these labels.

![Joint splitting counts and prime traces](Hecke_Prime_Atlas.png)

## 2. Symmetry and sector stability

**EXACT sector decomposition.** Rational CRT projectors were constructed from the squarefree minimal polynomial

\[
m(x)=x(x+2)g_4(x)h_6(x).
\]

They were verified to be idempotent, pairwise orthogonal, complete, and equal to the projectors obtained independently from the packaged integral block bases. Their real ranks are 2, 4, 8, 12. The block called “ghost” in the package is the level-11 oldspace, recorded here as “old”; the quartic block is a distinct newform sector.

The simultaneous Atkin–Lehner intersections are:

| Hecke sector | (+,+) | (+,−) | (−,+) | (−,−) | Total |
|---|---:|---:|---:|---:|---:|
| Elliptic | 2 | 0 | 0 | 0 | 2 |
| Oldspace | 0 | 0 | 2 | 2 | 4 |
| Quartic | 0 | 0 | 8 | 0 | 8 |
| Sextic | 0 | 12 | 0 | 0 | 12 |

The signs refer to (W11,W13); dimensions are real homology dimensions, twice the corresponding complex dimensions. An operator commuting with both involutions preserves each sign space. The **(−,+)** space contains two oldspace dimensions and all eight quartic dimensions. It can therefore mix those two Hecke sectors while preserving both involutions. Elliptic and sextic sectors have no such overlap with another Hecke sector.

### Controlled conservative model

**DIAGNOSTIC numerical experiment.** We whitened the frozen Hodge metric G, then used

\[
A_0=I+\frac{T_2}{4\|T_2\|_{\mathrm{op}}},\qquad
A=A_0+\epsilon V,\qquad F=JA,\qquad U(t)=e^{tF}.
\]

Here V is real symmetric, commutes with J, and has Frobenius norm 1 in the whitened coordinates. Then F is skew-symmetric and Hamiltonian, so U is orthogonal and symplectic. This deliberately restricts the study to conservative linear mixing. All time and frequency units are dimensionless; A0, perturbation strength, random ensemble, and horizon are chosen experiment parameters, not physics predictions.

The five classes use the same random input per seed: scalar offsets on Hecke sectors, arbitrary within-sector perturbations, an average over the W11/W13 symmetry group, its oldspace–quartic cross-block part only, and an unrestricted perturbation within the stated J-commuting symmetric class. “Hecke commuting” is this specific family of scalar block offsets, not an exhaustive sample of the full commutant. A within-sector perturbation need not commute with every Hecke operator; it still preserves the four coarse sectors.

The planned grid contains 16 paired seeds, five ε values from 0 to 0.4, and 161 times from 0 to 80: **400 flows**. Every energy matrix stayed positive; the smallest eigenvalue was **__MIN_EIGEN__**. For a uniform source in sector i, we measured

\[
M_{ji}(t)=\frac{\|P_jU(t)P_i\|_F^2}{d_i}.
\]

Every column sums to one. Total leakage averages 1−M_ii with weights d_i/26. “Mean peak” means each seed's maximum on the fixed time grid, averaged across seeds. It is not a continuous-time supremum or a confidence bound.

| Perturbation, ε = 0.4 | Mean peak total leakage | Mean peak oldspace leakage | Largest frequency shift |
|---|---:|---:|---:|
__RESULT_TABLE__

The last column is the largest change among the 13 sorted frequencies over all 16 seeds; it is in the model's dimensionless units. Frequency movement can occur even when sector transfer is exactly forbidden. The numerical magnitude of leakage is model-dependent; the forbidden connections follow from the exact projectors and symmetry characters.

For both Atkin–Lehner-preserving classes, all forbidden transfers remained below **__FORBIDDEN__**. Under any orthogonal flow preserving these characters, a uniform oldspace source can lose at most **50%** of its population: its other two dimensions have character (−,−) and cannot enter the quartic sector. A uniform quartic source can lose at most **25%**, and dimension-weighted total leakage is at most **4/26 ≈ 15.38%**. These are rank bounds, not fitted percentages.

The general within-sector perturbations provide a useful additional distinction: preserving the four projectors is sufficient for block protection even when the individual T2 eigenmodes inside a block mix. Thus full Hecke commutation is sufficient but stronger than necessary for this particular coarse protection question.

![Exact symmetry intersections and controlled transfer](Hecke_Sector_Stability.png)

## 3. What the two images add

The recurrence image suggests a matrix lift of z_next = z² + c:

\[
Z_{k+1}=Z_k^2+C.
\]

For any fixed Hecke projector P,

\[
[P,Z^2+C]=[P,Z]Z+Z[P,Z]+[P,C].
\]

If Z0 and C preserve every sector, induction shows that every iterate does too. If Z0 and C belong to Q[T2], all iterates remain in that commutative algebra. The nonlinear square therefore preserves an existing block decomposition under these assumptions. Cross-sector terms require incompatible initial data or forcing. This is an exact algebraic observation about the proposed matrix lift; it is not a claim that scalar Mandelbrot dynamics encodes Hecke geometry. No escape-time or fractal experiment has been run.

The ham-sandwich image suggests a different question: can one spatial cut balance several sector densities simultaneously? The standard theorem bisects n finite absolutely continuous measures in R^n; its familiar three-object case is in R³. The image's unrestricted “three in R^n” wording needs that dimensional qualification. [Tao's exposition](https://terrytao.wordpress.com/2008/11/27/the-kakeya-conjecture-and-the-ham-sandwich-theorem/) also explains the polynomial version, which could offer more degrees of freedom for balancing multiple embedded densities. A discrete atomic embedding needs a specified rule for boundary mass or smoothing before importing such a statement.

This run tested an algebraic balance analogy, not a spatial cut: the “balanced old–quartic” perturbation has **P_i V P_i = 0 for every sector**, a stronger condition than zero block traces, and still transfers amplitude. An elementary example is

\[
P=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
V=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]

Both diagonal blocks vanish, while [P,V] ≠ 0 and the transition probability under exp(−itV) is sin²t. Equal balance of scalar quantities does not enforce invariant subspaces. This distinction makes a later spatial-balance experiment more informative: it should measure both density balance and operator coupling.

## 4. Verification, scope, and files

The atlas checked **2,504 factor products**, **2,500 unramified factor-degree patterns by independent polynomial gcd/Frobenius arithmetic**, **2,498 discriminant-parity identities**, and 14 individual witness-factor irreducibility checks. All 24 and 720 permutations were enumerated to verify the density formula independently. Rational identities certify the sector ranks and sign intersections. Input hashes and the original experiment plan match.

The largest residual among the simulation's conservation, commutation, pairing, and independent matrix-exponential checks was **__RESIDUAL__**, against tolerance 10⁻⁹. Recomputing transfers directly in the original integral cycle frame with the Hodge norm gave maximum error **__FRAME_ERROR__**. Seven inspector cases, including invalid requests, passed. These are exact arithmetic certificates plus numerical validation, not a Lean formalization or interval-arithmetic error enclosure.

The bundle includes scripts, frozen input data, exact projectors, field certificates, all prime records, all transfer arrays, pinned dependencies, validation results, figures in PNG/SVG, and provenance. The prior irregularity labels were reused from the earlier verified Bernoulli census, not reclassified here. See README.md for the reproduction order and inspect_hecke.py for individual primes, intervals, or sectors.

The main continuation suggested by these results is to ask which **actual geometric or arithmetic perturbations** populate the allowed oldspace–quartic channel. The present Gaussian experiment establishes the selection rule and its measurable consequence; a subsequent experiment can replace the chosen random perturbations with graph defects or other specified operators and compare their projections into that channel.
'''
    replacements={'__Q4_TABLE__':q4table,'__Q6_TABLE__':q6table,'__UNSEEN_TABLE__':unseen,
                  '__COMPLETE6__':', '.join(map(str,complete6)),'__RESULT_TABLE__':resulttable,
                  '__MIN_EIGEN__':f"{stability['minimum_perturbed_eigenvalue']:.6f}",
                  '__FORBIDDEN__':f"{stability['forbidden_AL_transfer_max']:.3g}",
                  '__RESIDUAL__':f"{validation['maximum_simulation_residual']:.3g}",
                  '__FRAME_ERROR__':f"{validation['original_cycle_frame_transfer_error']:.3g}"}
    for key,value in replacements.items(): text=text.replace(key,value)
    assert '__RESULT_' not in text
    (ROOT/'Hecke_Atlas_and_Stability_Report.md').write_text(text)

def main():
    rows=read('prime_atlas.json'); types=read('splitting_types.json'); atlas=read('atlas_summary.json')
    stability=read('stability_summary.json'); fields=read('field_certificate.json'); sectors=read('sector_certificate.json')
    validation=read('validation_results.json'); assert validation['status']=='passed'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelcolor':INK,
                         'xtick.color':INK,'ytick.color':INK,'axes.spines.top':False,'axes.spines.right':False})
    atlas_figure(rows,types,atlas); stability_figure(stability,sectors)
    report(rows,types,atlas,stability,fields,sectors,validation)
    dump('PROVENANCE.json',{'experiment':'MTFT Hecke atlas and stability','version':'0.1.0','date':'2026-09-06',
         'python':platform.python_version(),'dependencies':{'numpy':np.__version__,'scipy':scipy.__version__,
         'sympy':sympy.__version__,'matplotlib':matplotlib.__version__},'source':read('inputs/manifest.json'),
         'plan_erratum_sha256':hashlib.sha256((ROOT/'PLAN_ERRATUM.md').read_bytes()).hexdigest(),
         'prime_labels_sha256':hashlib.sha256((ROOT/'inputs/prime_labels.json').read_bytes()).hexdigest(),
         'epistemic_status':{'field_and_sector_certificates':'exact arithmetic plus stated mathematical arguments',
                            'flows':'double-precision diagnostic; 1e-9 validation tolerance',
                            'irregularity_overlay':'reused exact prior census; descriptive only',
                            'image_connections':'proved matrix lemma and balance counterexample; no spatial cut or fractal computation'},
         'sources':['https://www.jmilne.org/math/CourseNotes/FT.pdf',
                    'https://math.mit.edu/classes/18.785/2021fa/LectureNotes14.pdf',
                    'https://math.mit.edu/classes/18.785/2015fa/LectureNotes25.pdf',
                    'https://docs.sympy.org/latest/modules/polys/numberfields.html',
                    'https://terrytao.wordpress.com/2008/11/27/the-kakeya-conjecture-and-the-ham-sandwich-theorem/']})
    print('Rendered report and two figures (PNG and SVG).')

if __name__=='__main__':
    main()
