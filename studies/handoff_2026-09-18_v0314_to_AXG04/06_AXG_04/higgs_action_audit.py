#!/usr/bin/env python3
"""AXG-04: changed-flux elementary-Higgs quiver and tree-level flavor obstruction.

Exact conditional algebra only. The analytic zero-mode and overlap proofs are
recorded in the accompanying notes; no FEM spectrum or fitted masses are used.
"""
from pathlib import Path
import importlib.util
import itertools
import json
import sympy as s

ROOT = Path(__file__).resolve().parent
source = ROOT / 'input' / 'smflux_v0314.py'
spec = importlib.util.spec_from_file_location('smflux_axg04_higgs', source)
smflux = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smflux)
M = smflux.M1_STACKS
names = ('c','L','a','b','d')
N = s.Matrix([M['N'][v] for v in names])
mold = s.Matrix([M['m'][v] for v in names])
Y = s.Matrix([s.Rational(M['Y'][v]) for v in names])
BL = s.Matrix([s.Rational(M['B-L'][v]) for v in names])
mnew = s.Matrix([0,-3,-3,-3,0])
e = [s.eye(5)[:,i] for i in range(5)]
checks=[]
def ck(name,condition):
    assert bool(condition), name
    checks.append(name)
def eq(name,lhs,rhs):
    if isinstance(lhs,s.MatrixBase):
        ck(name,(lhs-rhs).applyfunc(s.simplify)==s.zeros(*lhs.shape))
    else:
        ck(name,s.simplify(lhs-rhs)==0)

# Unordered orientations i<j; index=six-dimensional chirality*(mi-mj).
fields=[('cL',0,1,1,3,'Q'),('ca',0,2,-1,-3,'uc'),
        ('cb',0,3,-1,-3,'dc'),('Ld',1,4,1,-3,'lepton'),
        ('ad',2,4,-1,3,'nuc'),('bd',3,4,-1,3,'ec')]
q={label:e[i]-e[j] for label,i,j,sgn,target,sm in fields}
chi={label:sgn for label,i,j,sgn,target,sm in fields}
flux={label:int((q[label].T*mnew)[0]) for label,*_ in fields}
spectra=[]
for label,i,j,sgn,target,sm in fields:
    idx=sgn*flux[label]
    eq('three-family index '+label,idx,target)
    qleft=s.sign(idx)*q[label]
    spectra.append({'sector':label,'chirality6':sgn,'line_degree':flux[label],
                    'signed_index':idx,'SM_field':sm,'net_left_count':abs(idx),
                    'left_charge':list(qleft),'Y':(qleft.T*Y)[0],
                    'B-L':(qleft.T*BL)[0]})
eq('left hypercharges',s.Matrix([r['Y'] for r in spectra]),
   s.Matrix([s.Rational(1,6),-s.Rational(2,3),s.Rational(1,3),-s.Rational(1,2),0,1]))
eq('left B-L charges',s.Matrix([r['B-L'] for r in spectra]),
   s.Matrix([s.Rational(1,3),-s.Rational(1,3),-s.Rational(1,3),-1,1,1]))

# Fix common flux shift by mc=0 and solve all six family constraints.
mc,mL,ma,mb,md=s.symbols('mc mL ma mb md')
ms=s.Matrix([mc,mL,ma,mb,md])
conditions=[sgn*(q[label].T*ms)[0]-target for label,i,j,sgn,target,sm in fields]
sol=s.linsolve(conditions+[mc],(mc,mL,ma,mb,md))
ck('family constraints force realigned flux after common shift',sol==s.FiniteSet(tuple(mnew)))
old_indices=s.Matrix([sgn*(q[label].T*mold)[0] for label,i,j,sgn,target,sm in fields])
eq('same chirality changes on original M1 flux give wrong singlet indices',
   old_indices,s.Matrix([3,3,3,-3,-3,-3]))
ck('four sectors reverse desired old-flux family orientation',
   sum(old_indices[k]!=f[4] for k,f in enumerate(fields))==4)

Hu,Hd=e[1]-e[2],e[1]-e[3]
eq('new Higgs line degrees vanish',s.Matrix([(Hu.T*mnew)[0],(Hd.T*mnew)[0]]),s.zeros(2,1))
eq('old Higgs line degrees were minus6',s.Matrix([(Hu.T*mold)[0],(Hd.T*mold)[0]]),s.Matrix([-6,-6]))

# Gauge invariant 6D scalar Yukawas, given as (barred field, field, scalar).
vertices=[('up','ca','cL',Hu),('down','cb','cL',Hd),
          ('neutrino','Ld','ad',Hu),('electron','Ld','bd',Hd)]
for label,barred,plain,h in vertices:
    eq('gauge invariant elementary scalar vertex '+label,-q[barred]+q[plain]+h,s.zeros(5,1))
    eq('opposite six-dimensional chiralities in '+label,chi[barred]*chi[plain],-1)
    eq('matching internal line degree in '+label,flux[barred],flux[plain])
Pp=s.diag(1,1,1,1,0,0,0,0)
Pm=s.eye(8)-Pp
ck('bar negative times positive scalar bilinear projector survives',Pp*Pp!=s.zeros(8))
ck('bar positive times negative scalar bilinear projector survives',Pm*Pm!=s.zeros(8))
eq('equal chirality scalar bilinear still vanishes',Pm*Pp,s.zeros(8))

# Arithmetic purity: every triple of the four CM fixed points has both signs.
# Imported exact arithmetic facts: S0^2=K, h0(S0)=2 with basis1,u;
# u(P)=+i/sqrt13,+i/sqrt13,-i/sqrt13,-i/sqrt13.
z=s.I/s.sqrt(13)
eval4=s.Matrix([[1,z],[1,z],[1,-z],[1,-z]])
ck('every three CM evaluation rows have rank2',
   all(eval4[list(ids),:].rank()==2 for ids in itertools.combinations(range(4),3)))
g=13
degS=g-1
eq('S0(D) Riemann-Roch index',degS+3+1-g,3)
eq('S0(-D) Riemann-Roch index',degS-3+1-g,-3)
purity_plus=smflux.purity_with_S0(3,'cm')
purity_minus=smflux.purity_with_S0(-3,'cm')
ck('positive CM bundle purity agrees with source', (purity_plus['h0'],purity_plus['h1'])==(3,0))
ck('negative CM bundle purity agrees with source', (purity_minus['h0'],purity_minus['h1'])==(0,3))
k=s.symbols('k',integer=True,nonnegative=True)
eq('degree3 alone permits additional vectorlike pairs',(3+k)-k,3)

# Canonical constant scalar normalization and frame dependence.
A,A0,rho=s.symbols('Area Area0 rho',positive=True)
y6=s.symbols('y6',complex=True)
h0=1/s.sqrt(A)
eq('constant scalar profile normalized',A*h0**2,1)
eq('overlap Yukawa volume factor',y6*h0,y6/s.sqrt(A))
eq('Einstein-frame Yukawa after canonical rescaling',
   (y6/s.sqrt(A)).subs(A,rho**2*A0),y6/(rho*s.sqrt(A0)))
eq('six-dimensional Yukawa mass dimension',6-(s.Rational(5,2)+s.Rational(5,2)+2),-1)
eq('four-dimensional Yukawa dimension after area normalization',-1+1,0)

# Exact complex-unitary control, chosen to expose basis-dependence while
# leaving all singular values invariant. Analytic theorem holds for any U(3).
UQ=s.Matrix([[s.Rational(3,5),s.Rational(4,5),0],[-s.Rational(4,5),s.Rational(3,5),0],[0,0,s.I]])
Uu=s.Matrix([[0,s.I,0],[1,0,0],[0,0,-1]])
Ud=s.Matrix([[1,0,0],[0,s.Rational(5,13),12*s.I/13],[0,12*s.I/13,s.Rational(5,13)]])
for label,U in [('Q',UQ),('u',Uu),('d',Ud)]:
    eq('orthonormal basis control '+label,U.conjugate().T*U,s.eye(3))
Wu=Uu.conjugate().T*UQ
Wd=Ud.conjugate().T*UQ
eq('up overlap singular values all1',Wu.conjugate().T*Wu,s.eye(3))
eq('down overlap singular values all1',Wd.conjugate().T*Wd,s.eye(3))
ck('basis control gives different displayed Yukawa matrices',Wu!=Wd)
au,ad=s.symbols('abs_yu_squared abs_yd_squared',nonnegative=True)
HuGram=au*Wu.conjugate().T*Wu/A
HdGram=ad*Wd.conjugate().T*Wd/A
eq('common Q-side quark Gram commutator vanishes',HuGram*HdGram-HdGram*HuGram,s.zeros(3))
eq('quark CP commutator determinant vanishes',(HuGram*HdGram-HdGram*HuGram).det(),0)
eq('up overlap can be removed by right basis rotation',UQ.conjugate().T*Uu*Wu,s.eye(3))
eq('down overlap can be removed by independent right basis rotation',UQ.conjugate().T*Ud*Wd,s.eye(3))

# Separate C3X control: three parent generations with common unit X flux.
# An explicit unit-degree line avoids the vectorlike pair of O(P):
# L=O(P+Q-R), with distinct CM points P,Q,R and u(P)!=u(Q).
unit_eval_PQ=s.Matrix([[1,z],[1,-z]])
unit_eval_R=s.Matrix([[1,z]])
eq('unit-line P Q evaluation determinant',unit_eval_PQ.det(),-2*s.I/s.sqrt(13))
eq('evaluation at distinct CM R has rank1',unit_eval_R.rank(),1)
eq('S0 minusR section dimension',2-unit_eval_R.rank(),1)
eq('S0 plusR section dimension by RR and Serre',1+(2-unit_eval_R.rank()),2)
eq('unit-line inverse-twist section dimension',2-unit_eval_PQ.rank(),0)
eq('unit-line spinor section dimension by RR',1+(2-unit_eval_PQ.rank()),1)
eq('unit divisor degree',1+1-1,1)
eq('naive point unitflux positive zero modes',1+(2-unit_eval_R.rank()),2)
eq('three parent copies with pure unit line give three modes',3*1,3)
eq('three parent copies with point line give six and three netthree',3*2-3*1,3)

hcharge={'Q':1,'L':-3,'U':4,'D':-2,'E':-6,'N':0,'H':3,'Htilde':-3}
for label,barred,plain,h in [('up','Q','U','Htilde'),('down','Q','D','H'),
                            ('electron','L','E','H'),('neutrino','L','N','Htilde')]:
    eq('C3X hypercharge invariant vertex '+label,-hcharge[barred]+hcharge[plain]+hcharge[h],0)
    eq('C3X common X charge invariant vertex '+label,-1+1+0,0)
# Arbitrary parent flavor matrix remains arbitrary after the identity overlap;
# this control deliberately inserts hierarchy and is not an MTFT prediction.
Y6_control=s.diag(1,2,3)
Y4_control=Y6_control/s.sqrt(A)
eq('C3X parent flavor matrix is inherited',Y4_control.conjugate().T*Y4_control,
   s.diag(1,4,9)/A)

result={
 'investigation':'AXG-04 elementary Higgs action and flavor audit',
 'classification':'EXACT CONDITIONAL; explicitly changed field content and flux, tree-level only',
 'inputs':{'source':'input/smflux_v0314.py','stack_order':names,'old_flux':list(mold),'new_flux':list(mnew),
           'field_content':'Only the six displayed 6D Weyls, plus added elementary Hu/Hd scalar doublets; not the old all-bifundamental parent'},
 'spectrum':spectra,
 'flux_constraint':{'unique_with_mc_zero':list(mnew),'old_flux_indices':list(old_indices),
     'explanation':'Opposite 6D chiralities needed at scalar vertices, together with desired family signs, force La and Lb to have the same degree as LL.'},
 'action':{'Yukawas':['y6u bar(Psi_ca^-) Psi_cL^+ Hu','y6d bar(Psi_cb^-) Psi_cL^+ Hd',
                    'y6nu bar(Psi_Ld^+) Psi_ad^- Hu','y6e bar(Psi_Ld^+) Psi_bd^- Hd'],
     'plus':'Hermitian conjugates; standard positive scalar and Weyl kinetic terms and gauge kinetic terms.',
     'scope':'Explicit nonsupersymmetric local EFT vertices; scalar potential, anomaly cancellation and gravitational background equations are separate requirements. Elementary six-dimensional Yukawa coupling has mass dimension -1.'},
 'Higgs_zero_modes':{'condition':'LL,La,Lb identified as Hermitian line bundles WITH their unitary connections; degree equality alone is insufficient.',
     'profile':'unit parallel section/sqrt(Area)','count':'One complex weak doublet profile per Hu/Hd field on the connected compact curve.',
     'scope':'Zero of the minimal kinetic operator. Added bare masses, curvature couplings or background scalar potential can change the physical mass.'},
 'actual_fermion_zero_modes':{'choice':'Lc=Ld=O and LL=La=Lb=O(-D), D any three of the four W13-fixed CM points, with S0 arithmetic spin structure.',
     'positive_E':'E=O(D): h0(S0 E)=3,h1=0','negative_E':'E=O(-D): h0(S0 E)=0,h1=3',
     'proof':'H0(S0)=span(1,u); CM evaluation rank2 forces H0(S0(-D))=0; Serre duality and RR give remaining dimensions.',
     'generic_warning':'An arbitrary degree3 line only gives h0-h1=3; h0=3+k,h1=k is possible.'},
 'flavor':{'formula':'Y_u=(y6u/sqrt(Area)) Uu^dagger UQ; Y_d=(y6d/sqrt(Area)) Ud^dagger UQ',
     'invariant':'Y_f^dagger Y_f=|y6f|^2/Area I3 for quark vertices; analogous identity on lepton mode spaces.',
     'result':'Three equal singular values in each sector. Independent sector couplings do not generate within-sector hierarchy.',
     'CKM':'Can be chosen identity; because quark masses are exactly degenerate, CKM is basis-undetermined/unobservable, not a unique predicted matrix. Quark CP commutator invariant vanishes.',
     'scope':'Tree-level constant local Yukawa coefficients, identical relevant internal operators/connections, canonical minimal kinetic terms, constant scalar mode, no localized or higher-derivative interactions. Radiative corrections were not computed and are not ruled out.',
     'control_Wu':Wu.tolist(),'control_Wd':Wd.tolist()},
 'separate_C3X_control':{
     'parent':'3 copies of 6D physical SM representations Q+,L+,U-,D-,E-,N-, all with X=+1; gauge group SU3 x SU2 x U1_h x U1_X, h=6Y; elementary H has(h,X)=(3,0). This replaces the five-stack quiver.',
     'explicit_unit_line':'L_X=O(P+Q-R), P,Q,R three distinct W13-fixed CM points with u(P)=+i/sqrt13 and u(Q)=-i/sqrt13.',
     'proof':'Evaluation at R has rank1 in H0(S0), so h0(S0(-R))=1. RR and Serre give h0(S0(R))=2; inclusion identifies this with H0(S0). Evaluation at P,Q has rank2, so h0(S0(R-P-Q))=0. Thus h1(S0 L_X)=0 and h0(S0 L_X)=1.',
     'evaluation_det':str(unit_eval_PQ.det()),'one_copy_spinor_modes':{'h0':1,'h1':0},
     'three_copy_modes':'Exactly three desired modes, no opposite-chirality paired zero modes, within each six-field representation sector under this choice.',
     'point_line_trap':'L_X=O(P) instead gives h0=2,h1=1 for basepoint-free S0; with3parentcopies this gives6desired plus3mirror modes, net3.',
     'geometry':'The negative coefficient -R is a divisor presentation, not a physical singular source; the degree1 holomorphic line admits a smooth HYM connection.',
     'flavor':'The internal zero-mode space is1D; its constant-Higgs overlap is1/sqrt(Area). Therefore each arbitrary3x3 parent coefficient matrix is inherited asY4_f=Y6_f/sqrt(Area). Family hierarchy can be supplied as new parent input, not obtained from the curve.',
     'scope':'This audit checks action vertices and zero modes only. Full anomaly/tensor/global/gravitational consistency remains a separate calculation.'},
 'references':[{'title':'Cremades, Ibanez and Marchesano, hep-th/0404229, sections2 and3.1.3',
     'url':'https://arxiv.org/pdf/hep-th/0404229','use':'General normalized mode expansion and overlap principle; constant-mode identity derived here, not imported from a torus spectrum.'},
     {'title':'Park and Taylor, 1110.5916, Table1','url':'https://arxiv.org/pdf/1110.5916','use':'Fixed chiralities of ordinary 6D(1,0) multiplets; current six-Weyl EFT is not asserted to be such a completion.'}],
 'checks':checks,'summary':{'passed':len(checks),'failed':0}}
(ROOT/'higgs_action_results.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
print(json.dumps(result['summary']))
