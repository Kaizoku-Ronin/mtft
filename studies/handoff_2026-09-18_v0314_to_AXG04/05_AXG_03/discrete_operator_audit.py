#!/usr/bin/env python3
"""AXG-03: integer axion dressings and the Spin x Z3 fermion anomaly.

The congruence test is Hsieh arXiv:1808.02881, equation (1.2). It diagnoses
the fermion sector alone; no global GS/topological cancellation is assumed.
"""
from pathlib import Path
import importlib.util
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('axg03_smflux',ROOT/'input/smflux_v0314.py')
sm=importlib.util.module_from_spec(spec);spec.loader.exec_module(sm)
names=('c','L','a','b','d')
N=[sm.M1_STACKS['N'][k] for k in names]
m=[sm.M1_STACKS['m'][k] for k in names]
Y=s.Matrix([s.Rational(sm.M1_STACKS['Y'][k]) for k in names])
BL=s.Matrix([s.Rational(sm.M1_STACKS['B-L'][k]) for k in names])
K=s.Matrix([[0,-2,1,1,0],[3,-4,0,0,1]])
k1=K[0,:];kd=K[1,:]-k1
Kh=s.Matrix.vstack(kd,3*k1)
g=s.Matrix([0,0,1,0,1])  # generator exp(2pi i g/3)
e=[s.eye(5)[:,i] for i in range(5)]
charges={'Q':e[0]-e[1],'u_c':e[2]-e[0],'d_c':e[3]-e[0],
         'L':e[4]-e[1],'e_c':e[3]-e[4],'nu_c':e[2]-e[4],
         'Hu':e[1]-e[2],'Hd':e[1]-e[3]}
checks=[]
def ck(name,ok):
    assert bool(ok),name
    checks.append(name)
def eq(name,a,b):
    if isinstance(a,s.MatrixBase):
        ck(name,(a-b).applyfunc(s.simplify)==s.zeros(*a.shape))
    else: ck(name,s.simplify(a-b)==0)
eq('new charge matrix',Kh,s.Matrix([[3,-2,-1,-1,1],[0,-6,3,3,0]]))
eq('same continuous kernel',Kh*s.Matrix.hstack(s.ones(5,1),Y,BL),s.zeros(2,3))
snf=smith_normal_form(Kh,domain=ZZ)
ck('Smith factors 1 and 3',[abs(snf[i,i]) for i in range(2)]==[1,3])
eq('order-three generator fixes both axion circles',Kh*g,s.Matrix([0,3]))
eq('generator lies in native-axion kernel',(kd*g)[0],0)
eq('tensor transgression flux',(Kh*s.Matrix(m)),s.Matrix([0,36]))

ops={
 'up_Yukawa':charges['Q']+charges['u_c']+charges['Hu'],
 'down_Yukawa':charges['Q']+charges['d_c']+charges['Hd'],
 'lepton_Yukawa':charges['L']+charges['e_c']+charges['Hd'],
 'neutrino_Dirac':charges['L']+charges['nu_c']+charges['Hu'],
 'Hu_Hd':charges['Hu']+charges['Hd'],
 'QQQL':3*charges['Q']+charges['L'],
 'uc_uc_dc_ec':2*charges['u_c']+charges['d_c']+charges['e_c'],
 'nu_c_nu_c':2*charges['nu_c'],
}
outops={}
for name,q in ops.items():
    zz=int((q.T*g)[0])%3
    if name in ['up_Yukawa','down_Yukawa','lepton_Yukawa','neutrino_Dirac']:
        eq(name+' still neutral',q,s.zeros(5,1))
    outops[name]={'charge':[int(x) for x in q],'Z3':zz}
for name,expected in [('Hu_Hd',2),('QQQL',1),('uc_uc_dc_ec',1)]:
    eq(name+' nontrivial Z3',outops[name]['Z3'],expected)
    solution=Kh.T.gauss_jordan_solve(ops[name])[0]
    ck(name+' no integer axion dressing',not all(v.is_Integer for v in solution))
    outops[name]['rational_dressing']=[str(v) for v in solution]
for name in ['Hu_Hd','QQQL','uc_uc_dc_ec']:
    solution=Kh.T.gauss_jordan_solve(3*ops[name])[0]
    ck('cube of '+name+' admits integer dressing',all(v.is_Integer for v in solution))
    outops[name]['cube_integer_dressing']=[str(v) for v in solution]
eq('leading pure holomorphic Higgs product cube',3*ops['Hu_Hd'],Kh.T*s.Matrix([0,-1]))
eq('Majorana B-L remains nonzero',(ops['nu_c_nu_c'].T*BL)[0],2)
ck('Majorana outside axion row span',Kh.T.row_join(ops['nu_c_nu_c']).rank()>Kh.rank())

# Enumerate ALL net chiral bifundamentals, including extra weak doublets.
ledger=[]
S1=S3=0
for i in range(5):
    for j in range(i+1,5):
        index=m[i]-m[j]
        if index==0: continue
        orientation=1 if index>0 else -1
        q=orientation*(e[i]-e[j])
        charge=int((q.T*g)[0])
        dim=N[i]*N[j]
        count=abs(index)*dim
        S1+=count*charge
        S3+=count*charge**3
        ledger.append({'sector':names[i]+names[j],'orientation':orientation,
                       'multiplicity':abs(index),'representation_dimension':dim,
                       'signed_integer_charge':charge,'Z3':charge%3,
                       'S1':count*charge,'S3':count*charge**3})
eq('linear discrete trace',S1,24)
eq('cubic discrete trace',S3,24)
qdict={names[i]:sm.Fr(int(g[i])) for i in range(5)}
eq('cubic trace matches package independent ledger',
   s.Rational(sm.abelian_anomaly_polynomial(sm.M1_STACKS['N'],sm.M1_STACKS['m'],qdict)),S3)
n=3
cubic_res=((n*n+3*n+2)*S3)%(6*n)
linear_res=(2*S1)%n
eq('Hsieh cubic congruence residue',cubic_res,12)
eq('equivalent cubic mod9 residue',S3%9,6)
eq('linear congruence vanishes',linear_res,0)
phase=s.frac(s.Rational(n*n+3*n+2,6*n)*S3)
eq('nonzero anomaly phase fraction',phase,s.Rational(2,3))
ck('fermion sector does not pass Spin x Z3 test',cubic_res!=0)
# Charge representatives may be changed by 3; the mod9 cubic condition is unchanged.
S3res=sum(row['multiplicity']*row['representation_dimension']*row['Z3']**3 for row in ledger)
eq('residue representatives give same mod9 anomaly',S3res%9,S3%9)
eq('inverse generator conjugates anomaly',(-S3)%9,3)

# An isolated repair control: add one positive-chirality 6D Weyl fermion in
# (det U_L)^(-1) z_a z_b. This is new matter, not a derived or complete parent.
qnew=k1.T
eq('repair character is an actual weak determinant power',qnew[1],-N[1])
eq('repair character is hypercharge neutral',(qnew.T*Y)[0],0)
eq('repair character is B-L neutral',(qnew.T*BL)[0],0)
eq('repair character is common-phase neutral',sum(qnew),0)
new_index=(qnew.T*s.Matrix(m))[0]
eq('repair fermion has twelve net modes',new_index,12)
new_discrete=int((qnew.T*g)[0])%3
eq('repair modes have Z3 charge one',new_discrete,1)
new_S1=S1+new_index*new_discrete
new_S3=S3+new_index*new_discrete**3
eq('repair passes cubic congruence',((n*n+3*n+2)*new_S3)%(6*n),0)
eq('repair passes linear congruence',(2*new_S1)%n,0)
eq('repair leaves gauge-inertia witness unchanged',
   (qnew.T*s.Matrix([1,0,0,0,-3]))[0],0)
repair={'added_6D_field':'one positive-chirality complex Weyl in (det U_L)^(-1) z_a z_b',
        'charge':list(qnew),'index':new_index,'Z3':new_discrete,
        'new_S1':new_S1,'new_S3':new_S3,'Spin_x_Z3_congruences':'pass',
        'scope':'Isolated repair of this discrete fermion test by new charged matter. Full 6D anomalies, spectrum, interactions and vacuum must be recomputed. This field does not change the restricted gauge-inertia obstruction.'}

# Classical free-field Hodge count, separate from any interacting flux spectrum.
genus=13
h=[1,2*genus,1]
eq('closed genus13 first Betti number',h[1],26)
eq('Euler characteristic',sum((-1)**i*h[i] for i in range(3)),-24)
mode_counts={'nonchiral_B2':{'4D_vectors':h[1],'4D_real_scalars':h[0]+h[2]},
             'one_chiral_B2':{'4D_vectors':genus,'4D_real_scalars':1}}
eq('self-dual vector count halves free pair',2*mode_counts['one_chiral_B2']['4D_vectors'],
   mode_counts['nonchiral_B2']['4D_vectors'])

result={
 'investigation':'AXG-03 discrete operators and fermion anomaly',
 'class':'EXACT charge algebra and conditional global-fermion anomaly test',
 'Khat':Kh.tolist(),'generator_cocharacter':list(g),'Smith_invariants':[1,3],
 'field_Z3_charges':{name:int((q.T*g)[0])%3 for name,q in charges.items()},
 'operators':outops,'fermion_ledger':ledger,
 'isolated_discrete_repair_control':repair,
 'discrete_anomaly':{'S1':S1,'S3':S3,'S3_mod9':S3%9,'S1_mod3':S1%3,
    'Hsieh_cubic_mod18':cubic_res,'phase_fraction':str(phase),
    'verdict':'nonzero fermion-sector Spin x Z3 anomaly',
    'scope':'does not include an unspecified global GS/Wu/topological or inflow sector; ordinary neutral scalars do not cancel this fermion anomaly automatically',
    'source':'https://arxiv.org/pdf/1808.02881, equation (1.2)'},
 'mode_counts':mode_counts,
 'mode_scope':'free ordinary tensor harmonic reduction on the closed curve; flux couplings, constraints, projections or masses must be derived before calling these interacting massless modes',
 'notes':['This is a control with a changed integer lattice, not AXG01 with a harmless rescaling.',
          'The Z3 selection rules are conditional on a quantum-consistent completion.',
          'No B_ij, proton lifetime, or Higgs vacuum is predicted.',
          'Index-zero vector-like matter does not alter these net chiral anomaly sums.'],
 'checks':checks,'summary':{'passed':len(checks),'failed':0}}
(ROOT/'discrete_operator_results.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
print(json.dumps(result['summary']))
