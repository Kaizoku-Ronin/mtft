#!/usr/bin/env python3
"""AXG-03: necessary tensor factorization tests for the cd-flipped lift.

No spectrum is declared complete. Neutral Weyl/tensor additions are controls.
Only ordinary invariant four-form Green-Schwarz products and one native scalar
shift k_delta are allowed. All arithmetic is exact over Q.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'input' / 'axg02_parent_anomaly_results.json'
data = json.loads(SOURCE.read_text())
I8 = s.sympify(data['zero_index_flip']['I8'])
fc, fl, fa, fb, fd = s.symbols('f_c f_L f_a f_b f_d')
C, L, p, q, c3 = s.symbols('c2_color c2_weak p1 p2 c3_color')
x, y, z, d, t, n = s.symbols('x y z d t n')
checks = []

def check(label, condition):
    assert bool(condition), label
    checks.append(label)

def eq(label, a, b):
    if isinstance(a, s.MatrixBase):
        check(label, (a-b).applyfunc(s.simplify) == s.zeros(*a.shape))
    else:
        check(label, s.expand(a-b) == 0)

def mat(M):
    return [[str(v) for v in row] for row in M.tolist()]

def gram(poly, independent):
    """For polynomial quadratic in algebraically independent inputs."""
    return s.hessian(poly, independent) / 2

# Rebuild the input from the specified representations, so the saved AXG-02
# polynomial is a cross-check rather than an opaque computational dependency.
f = [fc, fl, fa, fb, fd]
N = [3, 2, 1, 1, 1]
ch = [[3, 0, -C, c3/2, C*C/12], [2, 0, -L, 0, L*L/12],
      [1, 0, 0, 0, 0], [1, 0, 0, 0, 0], [1, 0, 0, 0, 0]]
rebuilt = 0
signed_dimension = 0
for i in range(5):
    for j in range(i+1, 5):
        coeff = {}
        for degree in (2, 4):
            coeff[degree] = sum(
                ch[i][a]*(-1)**b*ch[j][b]*(f[i]-f[j])**(degree-a-b)
                /s.factorial(degree-a-b)
                for a in range(degree+1) for b in range(degree-a+1))
        sign = -1 if (i, j) == (0, 4) else 1
        rebuilt += sign*(coeff[4]-p*coeff[2]/24
                         +N[i]*N[j]*(7*p*p-4*q)/5760)
        signed_dimension += sign*N[i]*N[j]
eq('rebuild entire AXG02 polynomial from Chern characters', rebuilt, I8)
eq('rebuild signed chiral representation dimension', signed_dimension, 18)

# Quotient by the scalar-cancellable ideal, and discard matter-invisible common
# phase by taking f_L=0. The other relative curvature coordinates are x,y,z.
delta = 3*fc - 2*fl - fa - fb + fd
sub = {fc:x, fl:0, fa:y, fb:z, fd:y+z-3*x}
eq('restriction lies in kernel of native shift', delta.subs(sub), 0)
Q8 = s.expand(I8.subs(sub))
eq('six-form color term is killed by quotient', s.diff(Q8,c3), 0)
eq('residual irreducible gravitational anomaly', s.diff(Q8,q), -s.Rational(1,80))

# W and S have explicitly fixed anomaly signs, rather than names assigned from
# a separate convention. S=L_8/8 is one signed real chiral two-form anomaly.
W = (7*p*p-4*q)/s.Integer(5760)
S = (7*q-p*p)/s.Integer(360)
eq('cd-flipped net Weyl gravitational dimension', Q8.subs({C:0,L:0,x:0,y:0,z:0}),18*W)
balanced = s.expand(Q8-18*W)
check('18 negative neutral Weyls remove both pure gravity terms',
      not balanced.has(q) and s.expand(balanced).coeff(p,2)==0)
with_controls = s.expand(Q8+n*W+t*S)
eq('tensor-aware p2 equation', s.diff(with_controls,q), (28*t-18-n)/1440)
tensor_balanced = s.expand(with_controls.subs(n,28*t-18))
eq('tensor-aware p2 cancellation', s.diff(tensor_balanced,q), 0)
eq('tensor-aware remaining p1 square', tensor_balanced.coeff(p,2),t/32)
eq('tensor family is one-parameter pure-gravity shift',tensor_balanced,balanced+t*p*p/32)

v3=s.Matrix([C,L,p])
A=gram(balanced.subs({x:0,y:0,z:0}),v3)
A_expected=s.Matrix([[s.Rational(1,4),s.Rational(1,2),s.Rational(1,16)],
                     [s.Rational(1,2),s.Rational(1,2),s.Rational(1,8)],
                     [s.Rational(1,16),s.Rational(1,8),0]])
eq('fixed pure SU3 SU2 gravity block',A,A_expected)
LA,DA=A.LDLdecomposition(hermitian=False)
eq('pure-block exact LDL reconstruction',A,LA*DA*LA.T)
eq('pure-block LDL pivots',DA,s.diag(s.Rational(1,4),-s.Rational(1,2),-s.Rational(1,64)))
eq('pure-block determinant',A.det(),s.Rational(1,512))
check('one ordinary nonchiral tensor product cannot have rank three',A.rank()==3)

# Independent forms on the 1D Abelian restriction y=z=0 are C,L,p,X=x^2.
X=s.symbols('X')
restricted=s.Poly(s.expand(balanced.subs({y:0,z:0})+d*p*p),x)
R4=sum(coef*X**(power[0]//2) for power,coef in restricted.terms())
check('one-direction restriction has only even x powers',all(k[0]%2==0 for k,_ in restricted.terms()))
v4=s.Matrix([C,L,p,X])
M=gram(R4,v4)
eq('rank-four restriction reconstructs polynomial',(v4.T*M*v4)[0],R4)
eq('rank-four determinant',M.det(),s.Rational(9,32)*(1-48*d))
LM,DM=M.subs(d,0).LDLdecomposition(hermitian=False)
eq('neutral-control exact LDL reconstruction',M.subs(d,0),LM*DM*LM.T)
eq('neutral-control inertia certificate pivots',DM,s.diag(s.Rational(1,4),-s.Rational(1,2),-s.Rational(1,64),144))
gauge_indices=[0,1,3]
G=M.extract(gauge_indices,gauge_indices)
LG,DG=G.LDLdecomposition(hermitian=False)
eq('gravity-independent gauge block LDL reconstruction',G,LG*DG*LG.T)
eq('gauge block requires two positive directions',DG,s.diag(s.Rational(1,4),-s.Rational(1,2),108))
eq('gauge block determinant',G.det(),-s.Rational(27,2))
schur=M[2:,2:]-M[2:,:2]*M[:2,:2].inv()*M[:2,2:]
eq('Schur complement exact expression',schur,s.Matrix([[d-s.Rational(1,64),-s.Rational(3,4)],[-s.Rational(3,4),108]]))
eq('Schur determinant threshold',schur.det(),108*d-s.Rational(9,4))
eq('threshold d value',s.solve(schur.det(),d)[0],s.Rational(1,48))
eq('signed tensor threshold',s.solve(schur.det().subs(d,t/32),t)[0],s.Rational(2,3))
check('singular threshold has rank three',M.subs(d,s.Rational(1,48)).rank()==3)

# Full Abelian quartic restriction and a nonunique formal Gram presentation.
# Six quadratic monomials have 21 symmetric Gram entries but only 15 quartic
# coefficients; their six syzygies make the full Gram inertia representation-
# dependent. The above 1D restriction has no such ambiguity.
ell=s.Matrix([s.diff(balanced,v).subs({C:0,L:0,p:0})/2 for v in v3])
quartic=s.expand(balanced.subs({C:0,L:0,p:0}))
residual=s.expand(quartic-(ell.T*A.inv()*ell)[0])
eq('full square completion identity',balanced,((v3+A.inv()*ell).T*A*(v3+A.inv()*ell))[0]+residual)
eq('positive Schur residual witness',residual.subs({x:1,y:0,z:0}),144)
quadratics=[x*x,y*y,z*z,x*y,x*z,y*z]
ab_gram=s.zeros(6)
for monomial,coefficient in s.Poly(quartic,x,y,z).terms():
    target=x**monomial[0]*y**monomial[1]*z**monomial[2]
    pairs=[(i,j) for i in range(6) for j in range(i,6)
           if s.expand(quadratics[i]*quadratics[j]-target)==0]
    denominator=sum(1 if i==j else 2 for i,j in pairs)
    check('quartic monomial has ordinary four-form factorization '+str(target),denominator>0)
    for i,j in pairs:
        ab_gram[i,j]=ab_gram[j,i]=coefficient/denominator
B=s.zeros(3,6)
for i in range(3):
    for j,monomial in enumerate(quadratics):
        B[i,j]=s.Poly(ell[i],x,y,z).coeff_monomial(monomial)
full=A.row_join(B).col_join(B.T.row_join(ab_gram))
forms=s.Matrix([C,L,p]+quadratics)
eq('explicit full formal four-form Gram presentation',(forms.T*full*forms)[0],balanced)
coeff_map=[]
mon4=[x**i*y**j*z**(4-i-j) for i in range(5) for j in range(5-i)]
for mon in mon4:
    coeff_map.append([s.Poly((1 if i==j else 2)*quadratics[i]*quadratics[j],x,y,z).coeff_monomial(mon)
                      for i in range(6) for j in range(i,6)])
eq('quartic coefficient map rank',s.Matrix(coeff_map).rank(),15)
eq('quartic Gram ambiguity dimension',21-s.Matrix(coeff_map).rank(),6)

result={
 'investigation':'AXG-03 tensor-factor necessary conditions',
 'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
 'class':'EXACT conditional algebra; not an anomaly-free parent construction',
 'quotient_substitution':{str(k):str(v) for k,v in sub.items()},
 'quotient_I8':str(Q8),
 'neutral_control':{'extra_signed_complex_Weyls':-18,'I8':str(balanced)},
 'tensor_backreaction':{'W':str(W),'S':str(S),'n_definition':'extra signed neutral complex Weyl count',
  't_definition':'number of real chiral tensors with anomaly +S minus number with anomaly -S',
  'p2_cancellation':'18+n=28t','p1_square_coefficient':'d=t/32',
  'scope':'no gravitinos or additional charged chiral fields; such fields alter these equations'},
 'pure_block':{'basis':['c2_color','c2_weak','p1'],'matrix':mat(A),'LDL_L':mat(LA),'LDL_D':mat(DA),'inertia':[1,2,0]},
 'one_direction_block':{'basis':['c2_color','c2_weak','p1','x^2'],'matrix':mat(M),
   'neutral_control_LDL_L':mat(LM),'neutral_control_LDL_D':mat(DM),
   'neutral_control_inertia':[2,2,0],'determinant':str(s.factor(M.det())),
   'Schur_complement':mat(schur),
   'inertia_cases':{'d<1/48':[2,2,0],'d=1/48':[2,1,1],'d>1/48':[3,1,0]},
   'tensor_integer_cases':{'t<=0':[2,2,0],'t>=1':[3,1,0]},
   'ordinary_nonchiral_product_bound':'one product X4*Y4 has rank <=2 and at most one eigenvalue of either sign'},
 'fixed_gauge_block':{'basis':['c2_color','c2_weak','x^2'],'matrix':mat(G),'LDL_L':mat(LG),'LDL_D':mat(DG),'inertia':[2,1,0]},
 'full_formal_factorization':{'basis':[str(v) for v in forms],'matrix':mat(full),
   'quartic':str(quartic),'square_completion_ell':[str(v) for v in ell],
   'square_completion_residual':str(residual),'positive_residual_witness':{'x':1,'y':0,'z':0,'value':144},
   'Gram_ambiguity_dimension':6,'scope':'one rational algebraic presentation; no kinetic, lattice or global consistency claim'},
 'conclusions':[
  'Scalar cancellation by k_delta cannot change this quotient polynomial.',
  'Irreducible p2 must be cancelled by actual chiral spectrum before ordinary tensor factorization.',
  'The neutral-18 control requires at least two tensor-pairing directions of each sign and at least rank four.',
  'The fixed gauge restriction always requires at least two positive directions in the displayed I8 convention.',
  'If imposing I8=positive_constant*Omega(X4,X4) with Omega signature (1,T), neutral changes alone cannot suffice.',
  'With the opposite overall tensor sign, d>=1/48 is necessary on this restriction, not sufficient for the full polynomial.',
  'The conventional (1,0) hyperino convention is positive W here; Park-Taylor use its overall negative and -Omega.',
  'Additional charged gauginos or charged matter alter the gauge block; no theorem here excludes those altered spectra.',
  'Formal indefinite factorization does not establish positive tensor/gauge kinetic energy, a quantized self-dual charge lattice, global anomalies, supersymmetric multiplets or a stabilized vacuum.'
 ],
 'references':[{'url':'https://arxiv.org/abs/1110.5916','sections':'2.2-2.4, 2.7',
                'use':'generalized scalar versus ordinary tensor cancellation; (1,T) pairing and separate kinetic-positivity requirements'},
               {'url':'https://arxiv.org/abs/1110.4639','sections':'2.2, equation 2.2',
                'use':'real self-dual field local anomaly is L/8; global anomaly is a separate condition'}],
 'checks':checks,'summary':{'passed':len(checks),'failed':0}
}
(ROOT/'tensor_factor_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result['summary']))
