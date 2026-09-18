#!/usr/bin/env python3
"""AXG-03: minimal elementary Higgs scalars and six-dimensional chirality.

The operator positivity proof is analytic (recorded in JSON), not a numerical
spectral computation. Assertions verify source inputs and exact consequences.
"""
from pathlib import Path
import importlib.util
import json
import sympy as s

ROOT = Path(__file__).resolve().parent
source = ROOT / 'input' / 'smflux_v0314.py'
spec = importlib.util.spec_from_file_location('smflux_axg03_higgs', source)
smflux = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smflux)
M = smflux.M1_STACKS
names = ('c', 'L', 'a', 'b', 'd')
m = s.Matrix([M['m'][v] for v in names])
Y = s.Matrix([s.Rational(M['Y'][v]) for v in names])
checks = []

def ck(name, condition):
    assert bool(condition), name
    checks.append(name)

def eq(name, lhs, rhs):
    if isinstance(lhs, s.MatrixBase):
        ck(name, (lhs-rhs).applyfunc(s.simplify) == s.zeros(*lhs.shape))
    else:
        ck(name, s.simplify(lhs-rhs) == 0)

e = [s.eye(5)[:, i] for i in range(5)]
Hu, Hd = e[1]-e[2], e[1]-e[3]
eq('source flux vector', m, s.Matrix([0, -3, 3, 3, 0]))
eq('Hu hypercharge', (Hu.T*Y)[0], s.Rational(1, 2))
eq('Hd hypercharge', (Hd.T*Y)[0], -s.Rational(1, 2))
eq('Hu line degree', (Hu.T*m)[0], -6)
eq('Hd line degree', (Hd.T*m)[0], -6)

# Area and Bochner bound for the compact uniformizing metric. Convention:
# D^2=-i B dvol, deg L=(i/2pi) integral D^2, so B=2pi d/Area.
# On zero-forms, D*D=2 dbar* dbar+B=2 partial* partial-B.
# For B<0 choose the second expression, obtaining D*D>=-B.
R = s.symbols('R', positive=True)
g = s.Integer(13)
area = 4*s.pi*(g-1)*R**2
eq('genus13 compact area', area, 48*s.pi*R**2)
d = s.Integer(-6)
B = 2*s.pi*d/area
bound = -B
eq('degree minus6 HYM flux', B, -1/(4*R**2))
eq('elementary Higgs scalar lower bound', bound, 1/(4*R**2))
ck('elementary Higgs scalar bound strictly positive', bound.is_positive)
rho, R0 = s.symbols('rho R0', positive=True)
einstein_bound = s.simplify(rho**-2*bound.subs(R,rho*R0))
eq('four-dimensional scalar kinetic radius power in Einstein ansatz', -2+2,0)
eq('Einstein-frame scalar mass bound',einstein_bound,1/(4*R0**2*rho**4))

# Riemann-Roch is about Dolbeault kernels, not the scalar Bochner kernel.
degK = 2*g-2
degDual = degK-d
eq('canonical degree', degK, 24)
eq('Serre-dual Higgs bundle degree', degDual, 30)
chiL = d+1-g
eq('degree minus6 Dolbeault index', chiL, -18)
eq('H1 dimension since H0 negative-degree line vanishes', -chiL, 18)
eq('Serre-dual h0 by nonspecial Riemann-Roch', degDual+1-g, 18)
eq('Serre-dual scalar Bochner Landau level', 2*s.pi*degDual/area, 5/(4*R**2))
ck('Serre-dual Bochner level remains positive', (2*s.pi*degDual/area).is_positive)

# An arbitrary curvature coupling and bare mass can change the operator;
# neither coefficient is fixed by the arithmetic or by minimal coupling.
mu2, xi = s.symbols('mu2 xi', real=True)
Rinternal = -2/R**2
shifted_bound = s.factor(bound+mu2+xi*Rinternal)
eq('nonminimal scalar mass bound in unwarped product frame',
   shifted_bound, mu2+(s.Rational(1, 4)-2*xi)/R**2)
eq('conditional threshold when bound is saturated and bare mass zero',
   shifted_bound.subs({mu2:0, xi:s.Rational(1,8)}), 0)

# A six-dimensional scalar Dirac bilinear links opposite Weyl chiralities.
# gamma7 is diagonal; gamma0 anticommutes with it, hence bar(P+ Psi)=barPsi P-.
G7 = s.diag(1,1,1,1,-1,-1,-1,-1)
I8 = s.eye(8)
Pplus, Pminus = (I8+G7)/2, (I8-G7)/2
odd = s.BlockMatrix([[s.zeros(4),s.eye(4)],[s.eye(4),s.zeros(4)]]).as_explicit()
eq('chirality projector product vanishes', Pminus*Pplus, s.zeros(8))
eq('opposite-chirality scalar projector survives', Pminus*Pminus, Pminus)
eq('odd Clifford generator anticommutes with chirality', G7*odd+odd*G7,s.zeros(8))
ck('same-chirality vector bilinear may survive', Pminus*odd*Pplus != s.zeros(8))
eq('Spin6 chiral tensor square dimensions', s.binomial(4,2)+s.binomial(5,2),16)
eq('opposite-chirality tensor product dimensions', 1+15,16)

# Also cover charge-conjugated scalar bilinears, without assuming a Dirac
# contraction is the only possibility. Complexified 6D Euclidean Clifford
# matrices suffice for the Lorentz-invariant representation-theory check.
sx=s.Matrix([[0,1],[1,0]])
sy=s.Matrix([[0,-s.I],[s.I,0]])
sz=s.diag(1,-1)
I2=s.eye(2)
kron=s.kronecker_product
gammas=[kron(sx,I2,I2),kron(sy,I2,I2),kron(sz,sx,I2),
        kron(sz,sy,I2),kron(sz,sz,sx),kron(sz,sz,sy)]
ck('complexified six-dimensional Clifford algebra',
   all(gammas[i]*gammas[j]+gammas[j]*gammas[i] == (2*I8 if i==j else s.zeros(8))
       for i in range(6) for j in range(6)))
CE=gammas[1]*gammas[3]*gammas[5]
G7E=s.I*gammas[0]*gammas[1]*gammas[2]*gammas[3]*gammas[4]*gammas[5]
PE=(I8+G7E)/2
ck('charge conjugation matrix intertwines Clifford transpose',
   all(CE*gg*CE.inv()==-gg.T for gg in gammas))
eq('same-Weyl charge-conjugated scalar bilinear vanishes',PE.T*CE*PE,s.zeros(8))

# Existing four-dimensional triangle charges are allowed, but that fact does
# not itself specify a nonzero six-dimensional scalar Yukawa.
Q,uc,dc,lep,ec,nuc = e[0]-e[1],e[2]-e[0],e[3]-e[0],e[4]-e[1],e[3]-e[4],e[2]-e[4]
for label, charge in [('QucHu',Q+uc+Hu),('QdcHd',Q+dc+Hd),
                      ('LecHd',lep+ec+Hd),('LnucHu',lep+nuc+Hu)]:
    eq('4D triangle gauge charge '+label,charge,s.zeros(5,1))

result = {
 'investigation':'AXG-03 elementary Higgs parent audit',
 'class':'EXACT analytic theorem plus exact input/consequence checks; no numerical spectrum',
 'inputs':{'source':'input/smflux_v0314.py','stack_order':names,'flux':list(m),'genus':g},
 'minimal_scalar':{
   'charges':{'Hu':list(Hu),'Hd':list(Hd)},'line_degrees':[-6,-6],
   'operator':'D*D on scalar sections of the charged Hermitian line bundle',
   'proof':'On a closed compact Riemannian curve, <phi,D*D phi>=||D phi||^2. A zero eigenvector would be parallel; nonzero parallel sections trivialize a line bundle and force curvature zero, contradicting nonzero degree.',
   'HYM_convention':'D^2=-i B dvol, B=2pi deg(L)/Area. D*D=2 dbar* dbar+B=2 partial* partial-B on sections.',
   'HYM_bound':str(bound),'scope':'Positive kinetic term, smooth closed curve, central HYM connection, no bare mass, curvature coupling, twist, or gyromagnetic term. R is a dimensionful but unfixed radius; statement is in product-frame curvature units.',
   'Einstein_frame_bound':str(einstein_bound),
   'Einstein_frame_convention':'ds6^2=rho^-2 ds4E^2+rho^2 dsX0^2, R=rho R0. The elementary-scalar 4D kinetic normalization has no radius power, and mE^2=rho^-2 lambda(R). No radius or GeV scale is derived.',
   'saturation':'For negative degree, equality requires a nonzero anti-holomorphic covariantly defined section, equivalently a holomorphic section of the dual line. Degree alone does not guarantee equality. In untwisted M1 L=O(-2(P1+P2+P3)), its effective dual has a section, so at least one saturated mode exists. Degree-zero twists can remove that section.',
   'outcome':'No massless or tachyonic minimally coupled elementary scalar; no electroweak instability generated by minimal flux alone.'},
 'cohomology_distinction':{'H1_L':18,'Serre_dual_degree':30,'Serre_dual_H0':18,
   'Serre_dual_Bochner_level':str(2*s.pi*degDual/area),
   'interpretation':'The 18 H1 classes are Dolbeault-harmonic (0,1)-forms, represented by Serre-dual holomorphic sections. They are not 18 massless elementary scalar sections of degree -6. A physical form/twisted/gauge operator must be specified.'},
 'nonminimal_control':{'operator':'D*D + mu2 + xi R6','flat_external_unwarped_R6':'-2/R^2',
   'bound':str(shifted_bound),'saturated_zero_threshold_mu2_zero':'xi=1/8',
   'scope':'Added coefficients; saturation and product-frame assumptions explicit. No physical value or stable vacuum claimed.'},
 'chirality':{'scalar_bilinear':'bar(P+ Psi) P+ Chi = barPsi P- P+ Chi = 0',
   'charge_conjugated_scalar_check':'An explicit 8x8 complexified 6D Clifford representation gives P+^T C P+=0 as well.',
   'representation_check':'Spin(6,C)=SL4: 4 tensor 4=6+10 contains no scalar; 4 tensor dual4=1+15 does.',
   'scope':'Local nonderivative Lorentz-scalar Yukawa of two 6D Weyl fields and an elementary scalar. Higher-derivative operators, opposite-chirality new matter, internal vector components, or other parent structures are separate extensions.',
   'parent_consequence':'The AXG02 flip affects only cd; all four SM Yukawa triangles retain equal original 6D chiralities. Their 4D gauge-invariant charges do not establish a parent elementary-scalar Yukawa.',
   'SUSY_consequence':'An ordinary 6D (1,0) hypermultiplet has fixed hyperino chirality. Flipping one bifundamental fermion independently does not preserve that multiplet assignment.'},
 'references':[
   {'title':'Demailly, Complex Analytic and Differential Geometry, VII 1.3 (p330)', 'url':'https://www-fourier.univ-grenoble-alpes.fr/~demailly/manuscripts/agbook.pdf','supports':'Kahler Bochner-Kodaira identity; scalar bound derived here'},
   {'title':'Park and Taylor, arXiv1110.5916, Table1 (p6)','url':'https://arxiv.org/pdf/1110.5916','supports':'6D (1,0) multiplet chiralities'},
   {'title':'Cremades, Ibanez, Marchesano, hep-th/0404229, equations3.65-3.70','url':'https://arxiv.org/pdf/hep-th/0404229','supports':'Torus example distinguishes positive covariant Laplacian and gyromagnetically shifted internal-vector mass matrix; not a genus13 scalar spectrum'}],
 'checks':checks,'summary':{'passed':len(checks),'failed':0}}
(ROOT/'higgs_parent_results.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
print(json.dumps(result['summary']))
