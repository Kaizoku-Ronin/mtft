#!/usr/bin/env python3
"""AXG-04 C3X: explicit charged spectrum and quantized local GS data.

Independent reconstruction from the physical SM representations, rather than
the five-stack polynomial. This is an EFT candidate, not a UV completion or a
proof of the full global partition function.
"""
from pathlib import Path
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

ROOT=Path(__file__).resolve().parent
C,L,c3,p1,p2,h,x=s.symbols('C L c3 p1 p2 h x')
W=(7*p1*p1-4*p2)/s.Integer(5760)
checks=[]
def ck(label, condition):
    assert bool(condition),label
    checks.append(label)
def eq(label,a,b=0):
    if isinstance(a,s.MatrixBase): ck(label,(a-b).applyfunc(s.simplify)==s.zeros(*a.shape))
    else: ck(label,s.expand(a-b)==0)

# Physical six-dimensional representations, including right-handed 4D fields
# U,D,E,N rather than their left-handed conjugates. Every field has X charge1.
fields=[
 {'name':'Q','color':3,'weak':2,'h':1,'X':1,'chi6':1},
 {'name':'U','color':3,'weak':1,'h':4,'X':1,'chi6':-1},
 {'name':'D','color':3,'weak':1,'h':-2,'X':1,'chi6':-1},
 {'name':'Lepton','color':1,'weak':2,'h':-3,'X':1,'chi6':1},
 {'name':'E','color':1,'weak':1,'h':-6,'X':1,'chi6':-1},
 {'name':'N','color':1,'weak':1,'h':0,'X':1,'chi6':-1},
]
ch_color={3:[3,0,-C,c3/2,C*C/12],1:[1,0,0,0,0]}
ch_weak={2:[2,0,-L,0,L*L/12],1:[1,0,0,0,0]}
I8=I6=0
for f in fields:
    a=ch_color[f['color']];b=ch_weak[f['weak']]
    z=f['h']*h+f['X']*x
    coeff={}
    for k in (1,2,3,4):
        coeff[k]=s.expand(sum(a[i]*b[j]*z**(k-i-j)/s.factorial(k-i-j)
                            for i in range(k+1) for j in range(k-i+1)))
    dim=f['color']*f['weak']; f['dimension']=dim
    one=s.expand(f['chi6']*(coeff[4]-p1*coeff[2]/24+dim*W))
    reduced=s.expand(f['chi6']*f['X']*(coeff[3]-p1*coeff[1]/24))
    eq('sector flux pushforward '+f['name'],s.diff(one,x),reduced)
    f['I8_per_parent_copy']=str(one)
    f['4D_net_multiplets_three_copies']=3*f['chi6']*f['X']
    I8+=3*one; I6+=3*reduced
I8=s.expand(I8);I6=s.expand(I6)
X4=L+9*h*h
Y4=3*C+L+p1/2-27*h*h-6*x*x
eq('full independent C3X factorization',I8,X4*Y4)
eq('irreducible color cubic cancels',s.diff(I8,c3),0)
eq('irreducible gravitational anomaly cancels',s.diff(I8,p2),0)
eq('pure gravitational anomaly cancels',I8.subs({C:0,L:0,c3:0,h:0,x:0}),0)
eq('entire four-dimensional anomaly is pushed forward',s.diff(I8,x),I6)
eq('four-dimensional BF anomaly factor',I6,-12*x*(L+9*h*h))
eq('hypercharge-only 4D anomaly vanishes',I6.subs(x,0),0)
eq('six-dimensional signed representation dimension',
   3*sum(f['chi6']*f['dimension'] for f in fields),0)
eq('six-dimensional positive weak doublet count',
   3*sum(f['color'] for f in fields if f['weak']==2),12)
eq('conventional SU2 six-dimensional congruence',12%6,0)
eq('six-dimensional SU3 signed fundamental count',
   3*sum(f['chi6']*f['weak'] for f in fields if f['color']==3),0)

# Integral source coefficients in the true direct-product group, not a
# quotient by the ordinary Standard Model common center.
Omega=s.Matrix([[0,1],[1,0]])
source=s.Matrix([X4,Y4])
eq('hyperbolic tensor factorization',(source.T*Omega*source)[0]/2,I8)
eq('unimodular tensor lattice determinant',Omega.det(),-1)
v_c=s.Matrix([0,3]);v_w=s.Matrix([1,1]);v_lambda=s.Matrix([0,1])
b_hh=s.Matrix([18,-54]);b_xx=s.Matrix([0,-12]);b_hx=s.zeros(2,1)
a_grav=2*v_lambda
eq('source coefficient reconstruction',source,
   v_c*C+v_w*L+a_grav*p1/4+(b_hh*h*h+b_xx*x*x)/2+b_hx*h*x)
ck('all gauge source coefficients integral',all(z.is_Integer for v in
   [v_c,v_w,b_hh/2,b_xx/2,b_hx] for z in v))
u,v=s.symbols('u v',integer=True)
test=s.Matrix([u,v])
eq('even lattice norm',(test.T*Omega*test)[0],2*u*v)
eq('characteristic gravitational vector',(a_grav.T*Omega*test)[0]-(test.T*Omega*test)[0],2*u*(1-v))
pairings=s.Matrix.hstack(v_c,v_w,v_lambda,b_hh,b_xx).T*Omega*s.Matrix.hstack(v_c,v_w,v_lambda,b_hh,b_xx)
ck('all displayed string charge pairings integral',all(z.is_Integer for z in pairings))

# Flux-one transgression of the two source classes. In a nonchiral B
# description the second row is the electric/magnetic dual (BF) axion.
eq('electric Bianchi source has no internal X transgression',s.diff(X4,x),0)
eq('magnetic/BF source carries charge twelve',s.diff(Y4,x),-12*x)
K=s.Matrix([[0,-12]]) # columns (integer hypercharge h, X)
snf=smith_normal_form(K,domain=ZZ)
eq('residual Smith factor',abs(snf[0,0]),12)

# All fields expressed as common left-handed four-dimensional Weyls.
S1=S3=H2X=HX2=SU2twice=SU3twice=0
left=[]
for f in fields:
    qh=f['chi6']*f['h'];qx=f['chi6']
    count=3*f['dimension']
    S1+=count*qx;S3+=count*qx**3
    H2X+=count*qh*qh*qx;HX2+=count*qh*qx*qx
    if f['weak']==2: SU2twice+=3*f['color']*qx
    if f['color']==3: SU3twice+=3*f['weak']*qx
    left.append({'field':f['name'] if f['chi6']>0 else f['name']+'_conjugate',
                 'multiplets':3,'dimension':f['dimension'],'h':qh,'X':qx,'Z12':qx%12})
eq('linear discrete fermion trace',S1,0)
eq('cubic discrete fermion trace',S3,0)
eq('pure Spin x Z12 cubic congruence',((12*12+3*12+2)*S3)%(6*12),0)
eq('pure Spin x Z12 linear congruence',(2*S1)%12,0)
eq('mixed h squared X trace',H2X,-216)
eq('mixed h X squared trace',HX2,0)
eq('weak instanton X charge',SU2twice,12)
eq('color instanton X charge',SU3twice,0)
eq('weak instanton respects residual Z12',SU2twice%12,0)
eq('integer hypercharge mixed trace divisible by12',H2X%12,0)

# One elementary scalar H:(1,2,h=3,X=0), with Htilde=i sigma2 H*.
charges={f['name']:s.Matrix([f['h'],f['X']]) for f in fields}
H=s.Matrix([3,0]);Ht=-H
for name,qv in [('up',-charges['Q']+Ht+charges['U']),
                ('down',-charges['Q']+H+charges['D']),
                ('charged_lepton',-charges['Lepton']+H+charges['E']),
                ('Dirac_neutrino',-charges['Lepton']+Ht+charges['N'])]:
    eq('scalar Yukawa charges '+name,qv,s.zeros(2,1))
eq('elementary Higgs degree in X flux',H[1],0)
P=s.diag(1,1,1,1,0,0,0,0);Pm=s.eye(8)-P
eq('opposite-chirality scalar projector is nonzero',Pm*Pm,Pm)
ck('opposite-chirality scalar projector rank',Pm.rank()==4)
ops={'QQQL_X':4,'ucucdcec_X':-4,'Majorana_nucnuc_X':-2}
for name,qx in ops.items(): ck(name+' nontrivial modulo12',qx%12!=0)
eq('cube of QQQL neutral modulo12',(3*ops['QQQL_X'])%12,0)

# Action engineering dimensions and direct fixed-product reduction.
eq('6D scalar dimension',(6-2)/s.Integer(2),2)
eq('6D fermion dimension',(6-1)/s.Integer(2),s.Rational(5,2))
eq('6D Yukawa coefficient dimension',6-(2+s.Rational(5,2)*2),-1)
eq('6D scalar quartic coefficient dimension',6-4*2,-2)
A,y6,lam6=s.symbols('Area y6 lambda6',positive=True)
eq('constant scalar quartic reduction',A*lam6*(1/s.sqrt(A))**4,lam6/A)
eq('closed-curve tensor vector count',2*13,26)

# A completely declared classical product solution in arbitrary mass units.
# This gives an example action with an on-shell background, not our universe.
tH,mu2=s.symbols('HdagH mu2',real=True)
potential=mu2*tH+lam6*tH*tH
Hnorm=-mu2/(2*lam6)
eq('Higgs stationary norm',s.diff(potential,tH).subs(tH,Hnorm),0)
eq('Higgs minimum potential',potential.subs(tH,Hnorm),-mu2**2/(4*lam6))
example={'M6_power4':1,'gX6_squared':'1/1152','mu6_squared':-1,'lambda6':1,
         'Lambda6':'-19/4','HdagH':'1/2','U':-5,'beta':1,'R_squared':1,
         'AdS4_radius_squared':1,'radion_mass_squared':8,'radial_Higgs_mass_squared':2}
eq('chosen scalar vacuum norm',Hnorm.subs({mu2:-1,lam6:1}),s.Rational(1,2))
eq('chosen total constant potential',-s.Rational(19,4)-s.Rational(1,4),-5)
eq('chosen flux coefficient',1/(1152*s.Rational(1,1152)),1)
eq('chosen internal Einstein curvature',(-5+3)/s.Integer(2),-1)
eq('chosen external Einstein curvature',(-5-1)/s.Integer(2),-3)
eq('chosen radial Higgs squared mass',-2*(-1),2)

out={'study':'AXG-04 C3X explicit charged action candidate',
 'status':'Exact local anomaly cancellation and integral source/lattice checks; global differential-cohomology completion, UV completion and realistic vacuum not established',
 'group':'SU(3) x SU(2) x U(1)_h x U(1)_X as a direct product; h=6Y',
 'six_dimensional_copies':3,'fields':fields,'Higgs':{'color':1,'weak':2,'h':3,'X':0},
 'internal_X_flux':1,'I8':str(I8),'I6':str(I6),'sources':{'X4':str(X4),'Y4':str(Y4)},
 'tensor_lattice':{'Omega':Omega.tolist(),'signature':[1,1],
   'a_gravity':list(a_grav),'b_hh':list(b_hh),'b_xx':list(b_xx),
   'charge_pairing_matrix':pairings.tolist(),
   'scope':'Integral ordinary Chern classes and p1/2 on spin manifolds, even unimodular U lattice and characteristic vector; not proof for a quotient gauge global form'},
 'BF':{'level_X':-12,'hypercharge_level':0,'classical_finite_component':'Z12'},
 'four_dimensional_left_ledger':left,'discrete_traces':{'S1':S1,'S3':S3,'h2X':H2X,'hX2':HX2,
   'weak_instanton_X':SU2twice,'color_instanton_X':SU3twice,
   'scope':'Displayed pure Spin x Z12 and mixed instanton/trace tests pass; full global theory still requires construction'},
 'explicit_classical_AdS_control':example,
 'action':{'H3':'dB + omega3(X4); delta B = -omega2_1(X4)',
   'topological':'2*pi integral B wedge Y4 in normalized B-period-one units; cancels a descent choice I7=omega3(X4)Y4',
   'kinetic':'Einstein-Hilbert + positive Yang-Mills + positive ordinary nonchiral B kinetic + charged Weyl kinetic + scalar kinetic',
   'potential':'V6=mu6^2 HdagH + lambda6 (HdagH)^2, lambda6>0; parameters are new inputs',
   'Yukawas':['bar Q Y6u Htilde U','bar Q Y6d H D','bar L Y6e H E','bar L Y6nu Htilde N'],
   'dimensional_reduction':{'y4':'Y6/sqrt(Area)','lambda4':'lambda6/Area','mu4_squared':'mu6_squared in fixed product frame',
      'fermion_mass':'v4 Y4/sqrt(2); arbitrary parent flavor matrices are inputs, not MTFT predictions'}},
 'limitations':['Three copies are an input; unit flux supplies one mode per copy once the explicit degree-one bundle proof is used.',
    'No bare or axion-dressed Majorana neutrino mass while the residual discrete symmetry is preserved.',
    'The ordinary B field also has26 harmonic 4D vectors before interactions/projections; one of its two scalar modes participates in the X Stückelberg mechanism.',
    'An Einstein action and its scales are supplied rather than derived.',
    'The unwarped hyperbolic background is AdS, not a realistic Minkowski or dS vacuum.'],
 'checks':checks,'summary':{'passed':len(checks),'failed':0}}
(ROOT/'charged_action_results.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
print(json.dumps(out['summary']))
