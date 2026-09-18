#!/usr/bin/env python3
"""AXG-04 exact six-sector chiral-parent anomaly controls.

All spectra, gauge subgroups and fluxes below are explicit new choices.
No vacuum, UV completion or uniqueness claim is inferred from factorization.
"""
from pathlib import Path
import itertools
import json
import sympy as s

ROOT=Path(__file__).resolve().parent
names=('c','L','a','b','d')
N=(3,2,1,1,1)
f=s.Matrix(s.symbols('f_c f_L f_a f_b f_d'))
C,c3,L,p,q=s.symbols('c2_color c3_color c2_weak p1 p2')
x,y,z,w,r,h,u=s.symbols('x y z w r h u')
m=s.Matrix([0,-3,-3,-3,0])
mu=m/3
mold=s.Matrix([0,-3,3,3,0])
k1=s.Matrix([[0,-2,1,1,0]])
k2=s.Matrix([[3,-4,0,0,1]])
hyper=s.Matrix([s.Rational(1,6),0,-s.Rational(1,2),s.Rational(1,2),-s.Rational(1,2)])
bl=s.Matrix([s.Rational(1,3),0,0,0,-1])
signs={(0,1):1,(0,2):-1,(0,3):-1,(1,4):1,(2,4):-1,(3,4):-1}
checks=[]

def check(label,ok):
    assert bool(ok),label
    checks.append(label)

def eq(label,a,b):
    if isinstance(a,s.MatrixBase):
        check(label,(a-b).applyfunc(s.simplify)==s.zeros(*a.shape))
    else:
        check(label,s.expand(a-b)==0)

def matrix(M):
    return [[str(v) for v in row] for row in M.tolist()]

def push(P,flux):
    return s.expand(sum(flux[i]*s.diff(P,f[i]) for i in range(5)))

ch=[[3,0,-C,c3/2,C*C/12],[2,0,-L,0,L*L/12],
    [1,0,0,0,0],[1,0,0,0,0],[1,0,0,0,0]]
W=(7*p*p-4*q)/s.Integer(5760)
sector={}
for i,j in itertools.combinations(range(5),2):
    delta=f[i]-f[j]
    coeff=[]
    for degree in range(5):
        coeff.append(s.expand(sum(ch[i][a]*(-1)**b*ch[j][b]*delta**e/s.factorial(e)
                    for a in range(degree+1) for b in range(degree-a+1)
                    for e in [degree-a-b])))
    sector[i,j]={'I8':s.expand(coeff[4]-p*coeff[2]/24+N[i]*N[j]*W),
                 'ch3_minus_p_ch1':s.expand(coeff[3]-p*coeff[1]/24)}
I8=s.expand(sum(chi*sector[i,j]['I8'] for (i,j),chi in signs.items()))
I6=push(I8,m)
ledger=[]
for (i,j),chi in signs.items():
    old_index=mold[i]-mold[j]
    new_index=chi*(m[i]-m[j])
    eq('unchanged net chiral sector '+names[i]+names[j],new_index,old_index)
    eq('three-copy unit-flux index '+names[i]+names[j],3*chi*(mu[i]-mu[j]),new_index)
    eq('sector index pushforward '+names[i]+names[j],
       push(chi*sector[i,j]['I8'],m),new_index*sector[i,j]['ch3_minus_p_ch1'])
    ledger.append({'sector':names[i]+names[j],'chirality_6D':chi,
                   'dimension':N[i]*N[j],'net_index_4D':new_index,
                   'natural_hypercharge':hyper[i]-hyper[j]})
eq('signed six-dimensional dimension zero',sum(chi*N[i]*N[j] for (i,j),chi in signs.items()),0)
eq('no pure gravitational polynomial',I8.subs(dict(zip(f,[0]*5))).subs({C:0,L:0,c3:0}),0)
eq('no irreducible p2',s.diff(I8,q),0)
eq('color cubic native scalar direction',s.diff(I8,c3),(k1*f)[0]/2)
eq('native scalar flux admissibility',(k1*m)[0],0)
eq('native scalar preserves hypercharge',(k1*hyper)[0],0)
eq('native scalar preserves B-L',(k1*bl)[0],0)
eq('native scalar is determinant character',k1[1],-N[1])

# Independent Cartan-weight reconstruction, without Chern-character products.
alpha,beta,gamma=s.symbols('alpha beta gamma')
roots=[[alpha,beta,-alpha-beta],[gamma,-gamma],[0],[0],[0]]
cartan={C:-(alpha*alpha+alpha*beta+beta*beta),c3:-alpha*beta*(alpha+beta),L:-gamma*gamma}
independent=0
for (i,j),chi in signs.items():
    weights=[f[i]-f[j]+a-b for a in roots[i] for b in roots[j]]
    part=chi*sum(t**4/24-p*t*t/48+W for t in weights)
    eq('Cartan independent sector '+names[i]+names[j],chi*sector[i,j]['I8'].subs(cartan),part)
    independent+=part
eq('Cartan independent total',I8.subs(cartan),independent)

# Exact quotient by the only color-cubic scalar direction. Common phase drops.
relative={f[0]:x,f[1]:0,f[2]:y,f[3]:z,f[4]:w}
quotient={f[0]:x,f[1]:0,f[2]:y,f[3]:-y,f[4]:w}
F=L+y*y
G=C+L/3+p/6-3*x*x/2-w*w/2-y*y/3
Q=s.expand(I8.subs(quotient))
eq('full quotient factorization',Q,F*G)
Q6=s.expand(I6.subs(quotient))
eq('four-dimensional quotient anomaly',Q6,-3*F*(3*x+w))
eq('factor flux derivative first factor',3*s.diff(F,x)+3*s.diff(F,w),0)
eq('factor flux derivative second factor',3*s.diff(G,x)+3*s.diff(G,w),-3*(3*x+w))
relative_r=s.expand(I8.subs(relative).subs(z,r-y))
X6=s.cancel((relative_r-F*G)/r)
check('scalar remainder is polynomial',s.denom(X6)==1 or not s.denom(X6).has(r,x,y,w))
eq('exact scalar plus tensor decomposition',relative_r,r*X6+F*G)
pure=s.factor(Q.subs({x:0,y:0,w:0}))
eq('pure weak color gravity product',pure,L*(C+L/3+p/6))
pureM=s.hessian(pure,[C,L,p])/2
eq('core pure block rank',pureM.rank(),2)

# Spectator control: two index-zero positive weak bifundamentals and four
# negative neutral complex Weyls. They are new matter, not silently omitted.
S8=s.expand(I8+sector[1,2]['I8']+sector[1,3]['I8']-4*W)
eq('spectator La index zero',m[1]-m[2],0)
eq('spectator Lb index zero',m[1]-m[3],0)
eq('spectator leaves complete 4D chiral polynomial',push(S8,m),I6)
eq('spectator pure gravity zero',S8.subs(dict(zip(f,[0]*5))).subs({C:0,L:0,c3:0}),0)
eq('spectator color scalar direction unchanged',s.diff(S8,c3),(k1*f)[0]/2)
SG=C+L/2-(3*x*x+w*w)/2-3*y*y/2+p/4
SQ=s.expand(S8.subs(quotient))
eq('spectator two-product decomposition',SQ,F*SG+y*y*(4*y*y/3-p/6))
Y=s.symbols('Y')
restricted=SQ.subs({x:0,w:0,y**4:Y*Y,y**2:Y})
MS=s.hessian(restricted,[C,L,p,Y])/2
eq('spectator rank-four witness determinant',MS.det(),s.Rational(1,576))
gaugeMS=MS.extract([1,0,3],[1,0,3])
lg,dg=gaugeMS.LDLdecomposition(hermitian=False)
eq('spectator fixed gauge inertia certificate',dg,s.diag(s.Rational(1,2),-s.Rational(1,2),s.Rational(4,3)))
eq('spectator LDL reconstruction',gaugeMS,lg*dg*lg.T)
check('spectator has two positive directions independent of pure gravity',sum(bool(v>0) for v in dg.diagonal())==2)

# Traditional SU2 conditions interpreted as quantum GS/string-charge tests,
# rather than an independent torsion bordism anomaly after GS cancellation.
n2_core=sum(chi*(N[j] if i==1 else N[i]) for (i,j),chi in signs.items() if 1 in (i,j))
eq('one parent has four net complex SU2 doublets',n2_core,4)
eq('core conventional SU2 mod6 residue',n2_core%6,4)
eq('spectator conventional SU2 mod6 residue',(n2_core+2)%6,0)
eq('three-copy conventional SU2 mod6 residue',(3*n2_core)%6,0)
eq('net SU3 triplets vanish',2-1-1,0)
eq('three copies with unit flux preserve target I6',push(3*I8,mu),I6)
eq('three copies with old magnitude flux would give nine families',3*abs(m[0]-m[1]),9)

# Primitive original-product cocharacter a-b, in the native scalar kernel.
cochar=s.Matrix([0,0,1,-1,0])
eq('primitive cocharacter kills native scalar',(k1*cochar)[0],0)
abelian_sub={C:0,L:0,c3:0,p:0,q:0,**dict(zip(f,cochar*h))}
c_core=s.expand(I8.subs(abelian_sub)).coeff(h,4)
c_spec=s.expand(S8.subs(abelian_sub)).coeff(h,4)
eq('core primitive Abelian quartic',c_core,-s.Rational(1,3))
eq('spectator primitive Abelian quartic',c_spec,-s.Rational(1,6))
eq('P3 primitive Abelian quartic',3*c_core,-1)
check('core required string norm nonintegral',not (8*c_core).is_Integer)
check('spectator required string norm nonintegral',not (8*c_spec).is_Integer)
check('P3 passes only this norm witness',(24*c_core).is_Integer)

# Constructive P3 subgroup candidate. This is a NEW chosen gauge group, with
# only the allowed subgroup backgrounds; no inherited U-center quotient.
# All physical 6D fermions have X=1; Ld/ad/bd are conjugated to standard L/N/E.
embedding=6*hyper*h+mu*u
eq('subgroup kills native scalar identically',(k1*embedding)[0],0)
A4=L+9*h*h
B4=3*C+L+p/2-27*h*h-6*u*u
P3=s.expand(3*I8.subs(dict(zip(f,embedding))))
eq('P3 subgroup exact factorization',P3,A4*B4)
eq('P3 subgroup gravitational polynomial zero',P3.subs({C:0,L:0,h:0,u:0}),0)
P3I6=s.expand(s.diff(P3,u))
eq('P3 unit X flux anomaly',P3I6,-12*u*A4)
eq('P3 reduced polynomial agrees full six-sector flux push',P3I6,I6.subs(dict(zip(f,embedding))))
eq('P3 tensor flux source first component',s.diff(A4,u),0)
eq('P3 tensor flux source second component',s.diff(B4,u),-12*u)
Omega=s.Matrix([[0,1],[1,0]])
eq('tensor pairing unimodular',Omega.det(),-1)
eq('tensor pairing hyperbolic signature trace',s.trace(Omega),0)
av=s.Matrix([0,2])
bc=s.Matrix([0,-3])
bL=s.Matrix([-1,-1])
bh=s.Matrix([18,-54])
bX=s.Matrix([0,-12])
source=av*p/4-bc*C-bL*L+bh*h*h/2+bX*u*u/2
eq('integral anomaly coefficient source vector',source,s.Matrix([A4,B4]))
eq('lattice factorization normalization',(source.T*Omega*source)[0]/2,P3)
v1,v2=s.symbols('v1 v2',integer=True)
vec=s.Matrix([v1,v2])
eq('characteristic gravitational vector parity expression',
   (av.T*Omega*vec)[0]-(vec.T*Omega*vec)[0],2*v1*(1-v2))
check('all diagonal Abelian b vectors even',all(int(v)%2==0 for vv in [bh,bX] for v in vv))
physical_charges={
 'Q':{'h':1,'X':1,'chi':1,'dimension':6},
 'U':{'h':4,'X':1,'chi':-1,'dimension':3},
 'D':{'h':-2,'X':1,'chi':-1,'dimension':3},
 'L':{'h':-3,'X':1,'chi':1,'dimension':2},
 'E':{'h':-6,'X':1,'chi':-1,'dimension':1},
 'N':{'h':0,'X':1,'chi':-1,'dimension':1}}
for name,pair,orient in [('Q',(0,1),1),('U',(0,2),1),('D',(0,3),1),
                         ('L',(1,4),-1),('E',(3,4),-1),('N',(2,4),-1)]:
    qphys=orient*(embedding[pair[0]]-embedding[pair[1]])
    eq('physical subgroup charge '+name,qphys,physical_charges[name]['h']*h+u)

result={
 'investigation':'AXG-04 explicit six-sector chiral parents',
 'class':'EXACT conditional algebra; spectra and fluxes chosen explicitly',
 'conventions':{'N':N,'stack_order':names,'positive_Weyl_gravity':str(W),
  'complex_Weyl_chirality':'The displayed +/- weights multiply the full Chern character anomaly.',
  'tensor_normalization':'I8=1/2 Omega(X4,X4); X4=a*p1/4-b_color*c2_color-b_weak*c2_weak+(b_hh/2)h^2+(b_XX/2)u^2'},
 'core':{'flux':list(m),'ledger':ledger,'I8':str(I8),'I6':str(I6),
  'native_scalar_charge':list(k1),'quotient_coordinates':{str(k):str(v) for k,v in quotient.items()},
  'F4':str(F),'G4':str(G),'native_X6':str(X6),
  'native_X6_coordinates':'f_L=0, f_c=x, f_a=y, f_b=r-y, f_d=w; r=k1.f',
  'pure_gravitational_anomaly':0,'net_SU2_complex_doublets':n2_core,
  'pure_block_matrix':matrix(pureM),'pure_block_rank':pureM.rank(),
  'quartic_norm_witness':{'cocharacter':list(cochar),'coefficient_f4':str(c_core),'required_b_dot_b':str(8*c_core),'verdict':'fails integral string norm necessary condition'}},
 'spectator_control':{'added_fields':['La positive Weyl','Lb positive Weyl','four neutral negative Weyls'],
  'I8':str(S8),'quotient':str(SQ),'F4':str(F),'G4':str(SG),
  'second_product':[str(y*y),str(4*y*y/3-p/6)],
  'pure_gravitational_anomaly':0,'net_SU2_complex_doublets':6,
  'restricted_Gram_matrix':matrix(MS),'restricted_Gram_inertia':[2,2,0],
  'fixed_gauge_basis':['c2_weak','c2_color','y^2'],'fixed_gauge_LDL':matrix(dg),
  'quartic_norm_witness':{'coefficient_f4':str(c_spec),'required_b_dot_b':str(8*c_spec),'verdict':'fails despite conventional SU2 count passing'}},
 'P3_original_product_control':{'copies':3,'flux':list(mu),'families_4D':3,
  'I8':str(3*I8),'quotient':str(3*F*G),'net_SU2_complex_doublets':12,
  'primitive_a_minus_b_norm':str(24*c_core),
  'scope':'Passes the displayed SU2 and one U1 necessary tests; full original-product charge lattice is not claimed consistent.'},
 'P3_subgroup_candidate':{
  'gauge_group':'SU(3) x SU(2) x U(1)_h x U(1)_X, direct product',
  'h_normalization':'h is primitive U1 class; physical electric h charges are 6 times SM hypercharge',
  'embedding':list(embedding),'physical_6D_charges_one_copy':physical_charges,
  'copies':3,'X_flux_per_copy':1,'families_4D':3,'I8':str(P3),
  'tensor_F4':str(A4),'tensor_G4':str(B4),'Omega':matrix(Omega),
  'a':list(av),'b_color':list(bc),'b_weak':list(bL),'b_hh':list(bh),'b_XX':list(bX),'b_hX':[0,0],
  'gravitational_characteristic':'a=(0,2) is characteristic in the even hyperbolic U lattice',
  'I6':str(P3I6),'tensor_flux_source':[0,-12],
  'continuous_survivor':'hypercharge; X can become massive via tensor reduction',
  'conditional_discrete_remnant':'Z12 for one primitive period-2pi axion with charge12, before additional fields/global identifications',
  'status':'Explicit integral local-GS candidate with a nonchiral tensor; not a proved global quantum completion, UV parent, SUSY model or stable vacuum.',
  'family_provenance':'Three 6D parent copies and unit flux are choices. They replace one parent and index-three flux.'},
 'global_condition_notes':[
  'The traditional net-complex-SU2-doublet congruence is mod6; it counts six-dimensional fields, not flux-generated four-dimensional families.',
  'Lee-Tachikawa explain that these traditional conditions are recovered from globally quantized GS cancellation; standalone torsion bordism is trivial for BSU2.',
  'The primitive U1 norm check follows by placing c1=u+v on spin S2xS2, integral c1^2=2; X4=b*c1^2/2 has string charge b and I8 coefficient b.b/8.',
  'Positive kinetic terms, global product-group GS definition/inflow, the scalar action and vacuum remain separate requirements.'
 ],
 'references':[
  {'url':'https://arxiv.org/abs/hep-ph/0102010','use':'known nonsupersymmetric Q+,L+,U-,D-,E-,N- assignment and traditional generation/doublet congruence'},
  {'url':'https://arxiv.org/abs/2012.11622','use':'modern GS quantization interpretation of six-dimensional global gauge conditions'},
  {'url':'https://arxiv.org/abs/1711.04777','use':'cocharacter and string-charge quantization; full supergravity sufficiency is not imported'}],
 'checks':checks,'summary':{'passed':len(checks),'failed':0}
}
(ROOT/'chiral_parent_results.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
print(json.dumps(result['summary']))
