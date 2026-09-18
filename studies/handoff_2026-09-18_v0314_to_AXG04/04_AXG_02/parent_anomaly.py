#!/usr/bin/env python3
"""AXG-02: exact anomaly audit of a specified 6D matter lift, not a UV completion."""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import sympy as s

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'input/smflux_v0314.py'
spec=importlib.util.spec_from_file_location('axg02_smflux',SOURCE)
sm=importlib.util.module_from_spec(spec)
spec.loader.exec_module(sm)
names=('c','L','a','b','d')
N=tuple(sm.M1_STACKS['N'][i] for i in names)
m=s.Matrix([sm.M1_STACKS['m'][i] for i in names])
Y=s.Matrix([s.Rational(sm.M1_STACKS['Y'][i]) for i in names])
BL=s.Matrix([s.Rational(sm.M1_STACKS['B-L'][i]) for i in names])
f=s.Matrix(s.symbols('f_c f_L f_a f_b f_d'))
c2c,c3c,c2L,p1,p2=s.symbols('c2_color c3_color c2_weak p1 p2')
t=s.symbols('t')
K=s.Matrix([[0,-2,1,1,0],[3,-4,0,0,1]])
checks=[]
def check(label,condition):
    assert bool(condition), label
    checks.append(label)
def eq(label,lhs,rhs):
    if isinstance(lhs,s.MatrixBase):
        check(label,(lhs-rhs).applyfunc(s.simplify)==s.zeros(*lhs.shape))
    else: check(label,s.expand(lhs-rhs)==0)
def push(P):
    return s.expand(sum(m[i]*s.diff(P,f[i]) for i in range(5)))

check('M1 stack dimensions',N==(3,2,1,1,1))
eq('M1 flux',m,s.Matrix([0,-3,3,3,0]))
# Coefficients of the SU(N) Chern characters, degree k means a 2k-form.
ch=[
    [3,0,-c2c,c3c/2,c2c**2/12],
    [2,0,-c2L,0,c2L**2/12],
    [1,0,0,0,0],[1,0,0,0,0],[1,0,0,0,0],
]
pairs=list(itertools.combinations(range(5),2))
sector={}
for i,j in pairs:
    delta=f[i]-f[j]
    coeff=[]
    for k in range(5):
        coeff.append(s.expand(sum(ch[i][a]*(-1)**b*ch[j][b]*delta**e/s.factorial(e)
                                 for a in range(k+1) for b in range(k-a+1)
                                 for e in [k-a-b])))
    I8=s.expand(coeff[4]-p1*coeff[2]/24+N[i]*N[j]*(7*p1**2-4*p2)/5760)
    I6=s.expand((m[i]-m[j])*(coeff[3]-p1*coeff[1]/24))
    eq('sector pushforward '+names[i]+names[j],push(I8),I6)
    sector[i,j]={'I8':I8,'I6':I6,'dim':N[i]*N[j]}

I8=s.expand(sum(v['I8'] for v in sector.values()))
I6=s.expand(sum(v['I6'] for v in sector.values()))
eq('total pushforward identity',push(I8),I6)
eq('complex Weyl representation dimension',sum(v['dim'] for v in sector.values()),24)
C0=5*f[0]-2*f[1]-f[2]-f[3]-f[4]
U2=sum(N[i]*N[j]*(f[i]-f[j])**2 for i,j in pairs)
U4=sum(N[i]*N[j]*(f[i]-f[j])**4 for i,j in pairs)
Sc=sum(N[j]*(f[0]-f[j])**2 for j in range(1,5))
SL=sum(N[j]*(f[1]-f[j])**2 for j in range(5) if j!=1)
closed=(U4/24-c2c*Sc/2-c2L*SL/2+C0*c3c/2
        +s.Rational(5,12)*c2c**2+c2L**2/2+c2c*c2L
        -p1*U2/48+p1*(5*c2c+6*c2L)/24+(7*p1**2-4*p2)/240)
eq('compact hand-replication formula for I8',I8,closed)
eq('six-dimensional color cubic obstruction',s.diff(I8,c3c),C0/2)
eq('six-dimensional irreducible gravitational coefficient',s.diff(I8,p2),-s.Rational(1,60))
eq('color obstruction flux contraction zero',push(C0),0)
eq('pure color cubic vanishes after reduction',s.diff(I6,c3c),0)
eq('p2 vanishes after reduction',s.diff(I6,p2),0)
eq('C0 hypercharge overlap',C0.subs(dict(zip(f,Y))),s.Rational(4,3))
eq('C0 B-L overlap',C0.subs(dict(zip(f,BL))),s.Rational(8,3))

# Independent Cartan-weight construction.
u,v,w=s.symbols('u v w')
roots=[[u,v,-u-v],[w,-w],[0],[0],[0]]
I8_weights=0
for i,j in pairs:
    weights=[f[i]-f[j]+a-b for a in roots[i] for b in roots[j]]
    term=sum(z**4/s.Integer(24)-p1*z**2/s.Integer(48)
             +(7*p1**2-4*p2)/s.Integer(5760) for z in weights)
    sub={c2c:-(u*u+u*v+v*v),c3c:-u*v*(u+v),c2L:-w*w}
    eq('independent weights '+names[i]+names[j],sector[i,j]['I8'].subs(sub),term)
    I8_weights+=term
eq('independent full Cartan route',I8.subs(sub),I8_weights)

# AXG-01 four-dimensional normalization, rebuilt from the index formula.
P=sum((m[i]-m[j])*N[i]*N[j]*(f[i]-f[j])**3 for i,j in pairs)
Agr=sum((m[i]-m[j])*N[i]*N[j]*(f[i]-f[j]) for i,j in pairs)
A3=sum(s.Rational(1,2)*N[j]*(m[0]-m[j])*(f[0]-f[j]) for j in range(1,5))
A2=sum(s.Rational(1,2)*N[j]*(m[1]-m[j])*(f[1]-f[j]) for j in range(5) if j!=1)
eq('reproduces AXG01 full I6',I6,P/6-2*A3*c2c-2*A2*c2L-Agr*p1/24)

# A zero-index chirality change leaves the chiral pushforward invariant.
I8_flip=s.expand(I8-2*sector[0,4]['I8'])
Cdelta=(K*f)[1]-(K*f)[0]
eq('cd sector has zero index',m[0]-m[4],0)
eq('ab sector has zero index',m[2]-m[3],0)
eq('flipped cd still reproduces full I6',push(I8_flip),I6)
eq('flipped color anomaly aligns with K',s.diff(I8_flip,c3c),Cdelta/2)
eq('flipped color direction preserves Y',Cdelta.subs(dict(zip(f,Y))),0)
eq('flipped color direction preserves B-L',Cdelta.subs(dict(zip(f,BL))),0)
eq('flipped color direction preserves flux',push(Cdelta),0)
eq('flipped gravitational p2 coefficient',s.diff(I8_flip,p2),-s.Rational(1,80))
remainder=s.expand(I8_flip-Cdelta*c3c/2)
check('color six-form removed from remainder',not remainder.has(c3c))

# Ordinary two-form GS products contain only p1,c2c,c2L and Abelian quadratic
# four-forms, never the independent c3c six-form or irreducible p2 eight-form.
fourform_basis=[p1,c2c,c2L]+[f[i]*f[j] for i in range(5) for j in range(i,5)]
check('ordinary GS basis has no c3 or p2',all(not z.has(c3c,p2) for z in fourform_basis))
check('literal parent outside ordinary GS span',s.diff(I8,c3c)!=0 and s.diff(I8,p2)!=0)

# General bifundamental-only ordinary-GS obstruction, allowing signed 6D multiplicities.
n=s.symbols('n_cL n_ca n_cb n_cd')
Cgeneral=s.expand(sum(N[j]*n[j-1]*(f[0]-f[j]) for j in range(1,5)))
solutions=s.solve([s.diff(Cgeneral,fi) for fi in f],n,dict=True)
check('ordinary GS forces all net color bifundamentals zero',solutions==[dict.fromkeys(n,0)])
eq('quark index forced zero in that class',((m[0]-m[1])*n[0]).subs(solutions[0]),0)

result={
 'investigation':'AXG-02 parent anomaly',
 'class':'EXACT conditional algebra',
 'assumptions':['one complex 6D Weyl fermion of common chirality per unordered bifundamental pair',
                'central line fluxes only; no internal R-symmetry twist or boundary',
                'gauge bosons and ordinary metric are nonchiral; no gravitino, gaugino or tensor spectrum assumed'],
 'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
 'stack_order':names,'N':N,'m':list(m),
 'I8':str(I8),'I6':str(I6),
 'sector_polynomials':{names[i]+names[j]:{k:str(v) for k,v in data.items()} for (i,j),data in sector.items()},
 'obstruction':{'color_cubic_coefficient':str(C0/2),'p2_coefficient':'-1/60',
                'flux_contraction_color':0,'C0_Y':'4/3','C0_BL':'8/3',
                'ordinary_GS_scope':'products of gauge/Lorentz invariant four-forms from two-form tensor fields'},
 'zero_index_flip':{'changed_sector':'(c,dbar): 6D chirality + to -',
                    'I8':str(I8_flip),'color_cubic_coefficient':str(Cdelta/2),
                    'p2_coefficient':'-1/80','four_dimensional_chiral_pushforward':'unchanged',
                    'limitation':'not unchanged full zero-mode or massive spectrum, not anomaly-free parent'},
 'general_bifundamental_obstruction':{'coefficient':str(Cgeneral),
   'solution':{str(k):str(v) for k,v in solutions[0].items()},
   'scope':'independent five U1 factors, only ordinary two-form GS, color-charged net matter only standard bifundamentals'},
 'checks':checks,'summary':{'passed':len(checks),'failed':0},
}
(ROOT/'parent_anomaly_results.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
print(json.dumps(result['summary']))
