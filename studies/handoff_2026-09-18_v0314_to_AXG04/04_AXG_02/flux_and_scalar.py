#!/usr/bin/env python3
"""AXG-02: flux compatibility and conditional two-field potential, exact."""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parent
checks=[]
def ck(name,condition):
    assert bool(condition),name
    checks.append(name)
def eq(name,lhs,rhs):
    if isinstance(lhs,s.MatrixBase):
        ck(name,(lhs-rhs).applyfunc(s.simplify)==s.zeros(*lhs.shape))
    else: ck(name,s.simplify(lhs-rhs)==0)

K=s.Matrix([[0,-2,1,1,0],[3,-4,0,0,1]])
m=s.Matrix([0,-3,3,3,0])
U=s.Matrix([[-1,1],[1,0]])
Knew=U*K
eq('AXG01 shifts have nonzero internal flux',K*m,s.Matrix([12,12]))
eq('change of axion basis is unimodular',U.det(),-1)
eq('new shifts isolate flux-compatible direction',Knew*m,s.Matrix([0,12]))
eq('primitive compatible shift',Knew[0,:],s.Matrix([[3,-2,-1,-1,1]]))
p,q=s.symbols('p q')
eq('general row combination flux',((p*K[0,:]+q*K[1,:])*m)[0],12*(p+q))
ck('compatible subspace rank one',s.Matrix([list(K*m)]).rank()==1)

# Direct-sum product-group adjoint has no charged bifundamental scalar.
N=(3,2,1,1,1)
adjoint_weights=[]
for i,ni in enumerate(N):
    for a in range(ni):
        for b in range(ni):
            central=s.zeros(5,1)
            adjoint_weights.append(central)
ck('16 product-adjoint generators',len(adjoint_weights)==16)
ck('all product-adjoint central charges vanish',all(v==s.zeros(5,1) for v in adjoint_weights))
Hu=s.Matrix([0,1,-1,0,0]); Hd=s.Matrix([0,1,0,-1,0])
ck('Hu absent from product-adjoint charge weights',Hu not in adjoint_weights)
ck('Hd absent from product-adjoint charge weights',Hd not in adjoint_weights)

# Dimensional scaling of ds6^2=rho^-2 ds4_E^2+rho^2 dsX0^2.
eq('six-dimensional measure power',s.Rational(4*(-2)+2*2,2),-2)
eq('four-dimensional Einstein coefficient power',-2+2,0)
eq('internal curvature potential power',-2-2,-4)
eq('internal gauge flux potential power',-2-4,-6)
eq('bulk potential power',-2,-2)
eq('neutral scalar four-dimensional kinetic power',-2+2,0)
eq('common-coupling flux weight',sum(N[i]*m[i]**2 for i in range(5)),36)

# Positive bulk potential no-go at each fixed neutral scalar.
rho,a,b,c=s.symbols('rho a b c',positive=True)
sigma,alpha,beta=s.symbols('sigma alpha beta',real=True)
Vplus=a/rho**4+b*s.exp(alpha*sigma)/rho**6+c*s.exp(beta*sigma)/rho**2
eq('positive potential radial derivative',s.diff(Vplus,rho),
   -4*a/rho**5-6*b*s.exp(alpha*sigma)/rho**7-2*c*s.exp(beta*sigma)/rho**3)
ck('positive potential derivative strictly negative',(-s.diff(Vplus,rho)).is_positive)

# Negative exponential bulk potential: classify the complete two-field Hessian.
t=s.symbols('t',real=True)
V=a*s.exp(-4*t)+b*s.exp(alpha*sigma-6*t)-c*s.exp(beta*sigma-2*t)
A,B,C=s.symbols('A B C',positive=True)
H=s.hessian(V,(t,sigma))
# Use ratios to replace amplitudes in derivatives exactly.
replace={a:A*s.exp(4*t),b:B*s.exp(-alpha*sigma+6*t),c:C*s.exp(-beta*sigma+2*t)}
grad=s.Matrix([s.diff(V,t),s.diff(V,sigma)]).subs(replace,simultaneous=True).applyfunc(s.simplify)
Habc=H.subs(replace,simultaneous=True).applyfunc(s.simplify)
eq('two-field stationarity equations',grad,s.Matrix([-4*A-6*B+2*C,alpha*B-beta*C]))
eq('two-field Hessian before stationarity',Habc,s.Matrix([[16*A+36*B-4*C,-6*alpha*B+2*beta*C],
                                                       [-6*alpha*B+2*beta*C,alpha**2*B-beta**2*C]]))
r=s.symbols('r',positive=True)
betanonzero=s.symbols('beta_nonzero',nonzero=True,real=True)
stationary={A:(r-3)*B/2,C:r*B,alpha:r*betanonzero,beta:betanonzero}
eq('stationary ratio solution',grad.subs(stationary,simultaneous=True),s.zeros(2,1))
Hs=Habc.subs(stationary,simultaneous=True).applyfunc(s.simplify)
eq('stationary Hessian matrix',Hs,
   B*s.Matrix([[4*(r+3),-4*r*betanonzero],[-4*r*betanonzero,r*(r-1)*betanonzero**2]]))
eq('stationary Hessian determinant',Hs.det(),4*r*(r-3)*(r+1)*betanonzero**2*B**2)
eq('negative vacuum energy',(A+B-C).subs(stationary,simultaneous=True),-(r+1)*B/2)
eps=s.symbols('eps',positive=True)
ck('positive leading Hessian minor for r>3',Hs[0,0].subs(r,3+eps).is_positive)
ck('positive determinant for r>3',s.factor(Hs.det().subs(r,3+eps)).is_positive)
eq('unique log-coordinate solution determinant',s.Matrix([[-2,alpha],[2,beta]]).det(),-2*(alpha+beta))
# At beta=0, alpha!=0 the scalar derivative is nonzero; alpha=beta=0 is flat.
eq('beta zero scalar derivative',grad[1].subs(beta,0),alpha*B)
eq('both exponents zero scalar derivative',grad[1].subs({alpha:0,beta:0}),0)
eq('both exponents zero scalar Hessian',Habc[1,1].subs({alpha:0,beta:0}),0)

# Explicit rational control, deliberately chosen as a check rather than model fit.
control={a:1,b:2,c:8,alpha:4,beta:1,t:0,sigma:0}
eq('control gradient zero',s.Matrix([s.diff(V,t),s.diff(V,sigma)]).subs(control),s.zeros(2,1))
eq('control Hessian',H.subs(control),s.Matrix([[56,-32],[-32,24]]))
eq('control Hessian determinant',H.subs(control).det(),320)
eq('control energy',V.subs(control),-5)
eq('control eigenvalue minus',(H.subs(control)-(40-16*s.sqrt(5))*s.eye(2)).det(),0)
eq('control eigenvalue plus',(H.subs(control)-(40+16*s.sqrt(5))*s.eye(2)).det(),0)

result={
 'investigation':'AXG-02 flux and scalar audit','class':'EXACT under stated assumptions',
 'flux':{'K':K.tolist(),'m':list(m),'K_m':list(K*m),'U':U.tolist(),'U_K':Knew.tolist(),
         'UK_m':list(Knew*m),'compatible_primitive_row':list(Knew[0,:]),
         'proof':'For smooth periodic native 6D scalars theta, Dtheta is a global 1-form; dDtheta=K F, so Stokes gives K m=0 on the closed internal curve.',
         'scope':'fixed-modulus scalar, no defects or zeros; necessary topological condition only, not a vacuum or holonomy check',
         'hybrid_candidate':'one native 6D zero-form along k2-k1; second 4D axion would need another origin, e.g. reduction of a tensor; no such tensor sector derived'},
 'Higgs_representation':{'product_adjoint_dimension':16,'charged_bifundamental_scalars_from_product_gauge_components':False,
                         'scope':'elementary dimensional reduction of the stated product-group gauge fields; no claim excluding added charged scalars, composites, or a separately derived unified group'},
 'potential':{'ansatz':'a*rho^-4 + b*exp(alpha*sigma)*rho^-6 +/- c*exp(beta*sigma)*rho^-2; a,b,c>0',
              'positive_bulk':'no finite radial stationary point for any real sigma',
              'negative_bulk_beta_nonzero':'stationary point exists iff r=alpha/beta>3; unique in (log rho,sigma)',
              'amplitude_relations':'A=(r-3)B/2, C=rB; A=a*rho^-4, B=b exp(alpha sigma)rho^-6, C=c exp(beta sigma)rho^-2',
              'Hessian':str(Hs),'Hessian_det':'4*r*(r-3)*(r+1)*beta^2*B^2',
              'energy':'-(r+1)*B/2 < 0',
              'special_cases':'beta=0,alpha!=0: no scalar stationary point; alpha=beta=0: scalar remains flat',
              'scope':'two-field stationary minimum with positive scalar kinetic metric; conditional AdS when other equations satisfied; no parent-derived exponent or coefficient'},
 'control':{'parameters':{'a':1,'b':2,'c':8,'alpha':4,'beta':1},'rho':1,'sigma':0,
            'Hessian':[[56,-32],[-32,24]],'determinant':320,'V':-5,
            'label':'arbitrary-unit mathematical positive control, not MTFT prediction'},
 'checks':checks,'summary':{'passed':len(checks),'failed':0}}
(ROOT/'flux_and_scalar_results.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
print(json.dumps(result['summary']))
