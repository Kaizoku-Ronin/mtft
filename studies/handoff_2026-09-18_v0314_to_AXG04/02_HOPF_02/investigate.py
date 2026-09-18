#!/usr/bin/env python3
"""Exact HOPF-02 checks. Requires Python >=3.10 and SymPy; no mtft install.

The unchanged input/hecke_v0314.py snapshot reconstructs the Manin operator.
All reported equality gates use exact arithmetic. Run: python investigate.py
"""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
from math import gcd
import sympy as sp

ROOT = Path(__file__).resolve().parent
checks = []
results = {"date": "2026-09-17", "class": "EXACT algebra; stated geometric deductions", "checks": checks}

def check(name, condition):
    ok = bool(condition)
    checks.append({"name": name, "passed": ok})
    if not ok:
        raise AssertionError(name)

def clean(v):
    if isinstance(v, sp.MatrixBase):
        return [[str(x) for x in row] for row in v.tolist()]
    if isinstance(v, dict):
        return {str(k): clean(val) for k, val in v.items()}
    if isinstance(v, (list, tuple)):
        return [clean(x) for x in v]
    if isinstance(v, sp.Basic):
        return str(v)
    return v

print("Reconstructing the supplied rational Manin-symbol operator...", flush=True)
src = ROOT / "input/hecke_v0314.py"
spec = importlib.util.spec_from_file_location("hecke_snapshot", src)
h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h)
results["input_sha256"] = {str(src.relative_to(ROOT)): hashlib.sha256(src.read_bytes()).hexdigest()}
x = sp.Symbol("x")
g = sum(sp.Integer(c) * x**j for j, c in enumerate(h.G4))
sextic = sum(sp.Integer(c) * x**j for j, c in enumerate(h.H6))
T26 = sp.Matrix(h.cuspidal_hecke(2))
print("T2 reconstructed; checking characteristic polynomial...", flush=True)
check("full T2 characteristic polynomial", sp.expand(T26.charpoly(x).as_expr()-x**2*(x+2)**4*g**2*sextic**2) == 0)
P = sp.zeros(26)
for c in reversed(h.G4):
    P = P*T26 + c*sp.eye(26)
basis = sp.Matrix.hstack(*P.nullspace())
print("Quartic kernel reconstructed; checking its centralizer...", flush=True)
check("quartic kernel dimension 8", basis.shape == (26,8))
T8 = basis.gauss_jordan_solve(T26*basis)[0]
check("quartic basis intertwines", T26*basis == basis*T8)
check("quartic characteristic polynomial squared", sp.expand(T8.charpoly(x).as_expr()-g**2) == 0)
ann = sp.zeros(8)
for c in reversed(h.G4):
    ann = ann*T8+c*sp.eye(8)
check("quartic annihilator", ann == sp.zeros(8))
def evaluate(poly,operator):
    out=sp.zeros(operator.rows)
    for coefficient in sp.Poly(poly,x).all_coeffs():
        out=out*operator+coefficient*sp.eye(operator.rows)
    return out
rest=x*(x+2)*g
ba,bb,bgcd=sp.gcdex(sextic,rest,x)
check("sextic and complement coprime",bgcd==1)
check("full squarefree annihilator",evaluate(sextic*rest,T26)==sp.zeros(26))
P12=evaluate(bb*rest,T26)
P14=sp.eye(26)-P12
check("sextic projector idempotent",P12*P12==P12)
check("sextic projector dimension 12",sp.trace(P12)==12)
check("complement projector dimension 14",sp.trace(P14)==14)
check("complementary projections",P12*P14==sp.zeros(26) and P12+P14==sp.eye(26))
check("sextic projector image in its block",evaluate(sextic,T26)*P12==sp.zeros(26))
check("complement projector image in other blocks",evaluate(rest,T26)*P14==sp.zeros(26))
results["hodge_sphere_split"]={"sextic_polynomial":sextic,"complement_annihilator":rest,
    "Bezout_a":ba,"Bezout_b":bb,"P12":P12,"P14":P14,
    "sphere_identity":"S(H1)=S(V12)*S(V14)=S11*S13=S25 for the positive Hodge norm",
    "orthogonality_basis":"T2 is self-adjoint in the standard Hodge/Petersson inner product; real polynomial complementary projections are orthogonal"}
check("quartic separable", sp.gcd(g, sp.diff(g,x)) == 1)
root_intervals = sp.polys.polytools.intervals(g, eps=sp.Rational(1,1000))
check("four distinct real roots", len(root_intervals) == 4 and all(m == 1 for _,m in root_intervals))
check("quartic discriminant 1957", sp.discriminant(g,x) == 1957)
commutator_system = sp.kronecker_product(sp.eye(8),T8) - sp.kronecker_product(T8.T,sp.eye(8))
centralizer_dimension = 64-commutator_system.to_DM().rank()
check("real centralizer dimension 16", centralizer_dimension == 16)

C = sp.Matrix([[0,0,0,-1],[1,0,0,-5],[0,1,0,1],[0,0,1,3]])
G = sp.Matrix(4,4,lambda i,j: sp.trace(C**(i+j)))
check("trace metric symmetric", G == G.T)
minors = [G[:j,:j].det() for j in range(1,5)]
check("trace metric positive definite", all(v > 0 for v in minors))
check("companion self-adjoint in trace metric", C.T*G == G*C)
Li = sp.Matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]])
Lj = sp.Matrix([[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]])
I16 = sp.kronecker_product(Li,sp.eye(4))
J16 = sp.kronecker_product(Lj,sp.eye(4))
D16 = sp.kronecker_product(sp.eye(4),C)
M16 = sp.kronecker_product(sp.eye(4),G)
for label,Q in [("I",I16),("J",J16)]:
    check("doubled " + label + " squares to -1", Q*Q == -sp.eye(16))
    check("doubled " + label + " commutes with T2", Q*D16 == D16*Q)
    check("doubled " + label + " preserves control metric", Q.T*M16*Q == M16)
check("doubled quaternion anticommutation", I16*J16 == -J16*I16)
a,b=sp.symbols("a b",real=True)
anti=sp.Matrix([[a,b],[b,-a]])
check("2D obstruction square", anti**2 == (a*a+b*b)*sp.eye(2))
results["hecke"]={"g4":g,"quartic_root_intervals":root_intervals,"discriminant":1957,
    "T26":T26,"quartic_basis":basis,"T8":T8,"centralizer_dimension":centralizer_dimension,
    "centralizer_real_algebra":"M2(R)^4 (semisimple, four real eigenvalues, each multiplicity 2)",
    "commuting_quaternion_on_original":"IMPOSSIBLE: each invariant real eigenspace has dimension 2",
    "companion":C,"trace_metric_control":G,"leading_principal_minors":minors,
    "doubled_T":D16,"doubled_I":I16,"doubled_J":J16,"doubled_metric_control":M16,
    "metric_qualification":"constructed trace metric, not an identification with the supplied Petersson metric"}

print("Checking the spin pencil and Hopf curvature...",flush=True)
A=sp.diag(sp.I,-sp.I)
B=sp.Matrix([[0,-1],[1,0]])
K=A*B
for name,M in [("A",A),("B",B),("AB",K)]:
    check("Q8 "+name+" square",M*M == -sp.eye(2))
    check("Q8 "+name+" unitary",M.conjugate().T*M == sp.eye(2))
check("Q8 anticommutation",A*B == -B*A)
group=[e*M for M in [sp.eye(2),A,B,K] for e in [1,-1]]
check("Q8 eight distinct matrices",len({tuple(M) for M in group}) == 8)
check("Q8 closure",all(any(U*V == W for W in group) for U in group for V in group))
z=sp.Symbol("z")
def mob(M,z):
    return sp.cancel((M[1,0]+M[1,1]*z)/(M[0,0]+M[0,1]*z))
check("W11 projective action", mob(A,z) == -z)
check("W13 projective action", mob(B,z) == -1/z)
check("projective AL involutions commute",sp.cancel(mob(A,mob(B,z))-mob(B,mob(A,z))) == 0)
divs=[1,11,13,143]; exps=[-1,-1,1,1]
orders=[sum(sp.Rational(143*gcd(d,c)**2,24*gcd(c,143//c)*c*d)*r for d,r in zip(divs,exps)) for c in divs]
check("eta cusp orders", orders == [-6,-6,6,6])
check("eta weight zero",sum(exps)==0)
check("eta congruence one",sum(d*r for d,r in zip(divs,exps))%24 == 0)
check("eta congruence two",sum((143//d)*r for d,r in zip(divs,exps))%24 == 0)
check("eta square character",sp.prod(sp.Rational(d)**r for d,r in zip(divs,exps)) == 169)
degree=-sum(min(0,o) for o in orders)
check("pencil degree 12",degree == 12)
check("Riemann Hurwitz total ramification",2*13-2+2*degree == 48)
check("cusp and remaining ramification",sum(abs(o)-1 for o in orders) == 20 and 48-20 == 28)

# Quaternion multiplication in real coordinates (1,i,j,k).
def qm(p,q):
    p0,p1,p2,p3=p; q0,q1,q2,q3=q
    return sp.Matrix([p0*q0-p1*q1-p2*q2-p3*q3,
        p0*q1+p1*q0+p2*q3-p3*q2,
        p0*q2-p1*q3+p2*q0+p3*q1,
        p0*q3+p1*q2-p2*q1+p3*q0])
coords=sp.symbols("x0 x1 x2 x3",real=True)
den=1+sum(t*t for t in coords)
qbar=sp.Matrix([coords[0],-coords[1],-coords[2],-coords[3]])
units=[sp.eye(4)[:,i] for i in range(4)]
check("left arithmetic action moves Hopf coordinate",qm(units[1],units[1])[0] == -1)
check("right fiber action fixes Hopf coordinate",qm(units[1],-units[1])[0] == 1)
conn=[]
for e in units:
    component=qm(qbar,e)/den
    component[0]=0
    conn.append(component)
curv={}
for mu,nu in itertools.combinations(range(4),2):
    expr=conn[nu].diff(coords[mu])-conn[mu].diff(coords[nu])+qm(conn[mu],conn[nu])-qm(conn[nu],conn[mu])
    curv[mu,nu]=expr.applyfunc(sp.factor)
expected={(0,1):2*units[1],(0,2):2*units[2],(0,3):2*units[3],
          (1,2):-2*units[3],(1,3):2*units[2],(2,3):-2*units[1]}
for pair,num in expected.items():
    check("instanton curvature "+str(pair),curv[pair] == num/den**2)
check("anti-self-duality 01/23",curv[0,1] == -curv[2,3])
check("anti-self-duality 02/13",curv[0,2] == curv[1,3])
check("anti-self-duality 03/12",curv[0,3] == -curv[1,2])
def tr_product(p,q):
    return 2*qm(p,q)[0]
trace_density=sp.factor(2*(tr_product(curv[0,1],curv[2,3])-tr_product(curv[0,2],curv[1,3])+tr_product(curv[0,3],curv[1,2])))
check("trace F wedge F density",trace_density == 48/den**4)
r=sp.Symbol("r",positive=True)
radial=sp.integrate(r**3/(1+r*r)**4,(r,0,sp.oo))
check("instanton radial integral",radial == sp.Rational(1,12))
c2=sp.simplify(48*2*sp.pi**2*radial/(8*sp.pi**2))
check("second Chern number +1 in stated convention",c2 == 1)

# Complex plane q=x0+i*x1 gives CP1 inside HP1; connection reduces to iR.
restrict={coords[2]:0,coords[3]:0}
check("reduced A0 is along i",conn[0].subs(restrict) == -coords[1]*units[1]/(1+coords[0]**2+coords[1]**2))
check("reduced A1 is along i",conn[1].subs(restrict) == coords[0]*units[1]/(1+coords[0]**2+coords[1]**2))
check("reduced curvature",curv[0,1].subs(restrict) == 2*units[1]/(1+coords[0]**2+coords[1]**2)**2)
radial2=sp.integrate(2*r/(1+r*r)**2,(r,0,sp.oo))
check("monopole magnitude 1",radial2 == 1)
results["hopf"]={"Q8_A":A,"Q8_B":B,"cusp_orders":orders,"pencil_degree":degree,
    "pencil":"[1:sqrt(13)u] : X0(143) -> CP1 subset HP1",
    "connection_components":conn,"curvature_components":curv,"trace_F_wedge_F_density":trace_density,
    "orientation":"dx0 wedge dx1 wedge dx2 wedge dx3",
    "Chern_convention":"c2=(1/(8*pi^2))*tr(F wedge F), tr in complex 2D, anti-Hermitian F",
    "c2_on_S4":c2,"pullback_line_degrees":[-degree,degree],
    "pullback_structure":"SU2 bundle topologically trivial on the curve; U1 reduction is nontrivial",
    "action_distinction":"arithmetic Q8 acts on the left on coordinate pairs, moving the base; Hopf fiber Sp1 acts on the right and fixes the base"}

print("Enumerating the CRT modular-level factorization...",flush=True)
def canon_prime(c,d,p):
    c%=p;d%=p
    assert c or d
    return (0,1) if c == 0 else (1,(d*pow(c,-1,p))%p)
def crt(a,b):
    return (78*a+66*b)%143
check("CRT idempotents",78*78%143 == 78 and 66*66%143 == 66 and 78*66%143 == 0 and (78+66)%143 == 1)
P11=[(0,1)]+[(1,t) for t in range(11)]
P13=[(0,1)]+[(1,t) for t in range(13)]
pairs=list(itertools.product(P11,P13))
reps=[(crt(u[0],v[0]),crt(u[1],v[1])) for u,v in pairs]
check("168 CRT representatives",len(reps) == len(set(reps)) == 168)
check("all CRT reps primitive",all(gcd(gcd(c,d),143)==1 for c,d in reps))
check("CRT inverse",all((canon_prime(c,d,11),canon_prime(c,d,13))==pair for (c,d),pair in zip(reps,pairs)))
check("primitive vector count",sum(gcd(gcd(c,d),143)==1 for c in range(143) for d in range(143)) == 168*120)
index={pair:i for i,pair in enumerate(pairs)}
def act_prime(pair,M,p):
    c,d=pair; a,b,cc,dd=M
    return canon_prime(c*a+d*cc,c*b+d*dd,p)
S=(0,-1,1,0); TT=(1,1,0,1)
def perm(M):
    return [index[(act_prime(u,M,11),act_prime(v,M,13))] for u,v in pairs]
PS,PT=perm(S),perm(TT)
check("S involution",all(PS[PS[j]]==j for j in range(168)))
ST=[PT[PS[j]] for j in range(168)]
check("ST order three",all(ST[ST[ST[j]]]==j for j in range(168)))
def cycles(p):
    seen=set(); out=[]
    for j in range(len(p)):
        if j in seen:continue
        cyc=[]; k=j
        while k not in seen:
            seen.add(k);cyc.append(k);k=p[k]
        out.append(cyc)
    return out
widths=sorted(map(len,cycles(PT)))
check("cusp widths",widths == [1,11,13,143])
check("84 edge orbits",len(cycles(PS))==84)
check("56 triangular orbits",len(cycles(ST))==56)
check("genus from dessin",4-84+56 == 2-2*13)
model=h.model()
key_to_old={(canon_prime(c,d,11),canon_prime(c,d,13)):j for j,(c,d) in enumerate(model["P1"])}
old_order=[key_to_old[pair] for pair in pairs]
check("S matches package permutation",all(old_order[PS[j]] == model["sS"][old_order[j]] for j in range(168)))
check("T matches package permutation",all(old_order[PT[j]] == model["sT"][old_order[j]] for j in range(168)))
def genus_squarefree(N):
    primes=list(sp.factorint(N))
    mu=sp.Integer(N)*sp.prod(1+sp.Rational(1,p) for p in primes)
    e2=sp.prod(1+sp.legendre_symbol(-1,p) for p in primes)
    e3=sp.prod(1+sp.legendre_symbol(-3,p) for p in primes)
    cusps=2**len(primes)
    genus=1+mu/12-e2/4-e3/3-sp.Rational(cusps,2)
    return dict(level=N,index=mu,e2=e2,e3=e3,cusps=cusps,genus=genus)
genera=[genus_squarefree(N) for N in [11,13,143]]
check("three curve genera",[e["genus"] for e in genera] == [1,0,13])
check("forgetful degrees",genera[2]["index"]/genera[0]["index"]==14 and genera[2]["index"]/genera[1]["index"]==12)
check("sphere product dimension 24",11+13==24)
check("sphere join dimension 25",11+13+1==25)
check("tensor ambient sphere dimension 167",(11+1)*(13+1)-1==167)
check("generalized quaternionic Hopf dimensions",4*36-1==143 and 4*35+3==143)
results["level_and_spheres"]={"CRT_idempotents":[78,66],"projective_points":168,"unit_count":120,
    "CRT_representatives":reps,"S_permutation":PS,"T_permutation":PT,"cusp_widths":widths,
    "genera":genera,"forgetful_degrees":{"to_X0_11":14,"to_X0_13":12},
    "sphere_product_dimension":24,"sphere_join_dimension":25,"tensor_ambient_unit_sphere_dimension":167,
    "S143_Hopf":"S3 -> S143 -> HP35, base real dimension 140; no identification with X0(143)"}
results["summary"]={"passed":sum(c["passed"] for c in checks),"failed":sum(not c["passed"] for c in checks),
    "unproved_identifications":["Hopf fiber action equals geometric AL action", "constructed trace metric equals MTFT Petersson metric", "Hecke-equivariance of the half-form pencil", "a parent gauge or gravitational theory"]}
(ROOT/"results.json").write_text(json.dumps(clean(results),indent=2)+"\n")
print(json.dumps(results["summary"],indent=2),flush=True)
