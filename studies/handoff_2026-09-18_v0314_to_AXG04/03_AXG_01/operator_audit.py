#!/usr/bin/env python3
"""Exact charge tests for axion-dressed operators; no coefficient is predicted."""
from pathlib import Path
import json
import sympy as s

root=Path(__file__).resolve().parent
K=s.Matrix([[0,-2,1,1,0],[3,-4,0,0,1]])
Y=s.Matrix([s.Rational(1,6),0,-s.Rational(1,2),s.Rational(1,2),-s.Rational(1,2)])
BL=s.Matrix([s.Rational(1,3),0,0,0,-1])
e=[s.eye(5)[:,i] for i in range(5)]
q={"Q":e[0]-e[1],"u_c":e[2]-e[0],"d_c":e[3]-e[0],
   "L":e[4]-e[1],"e_c":e[3]-e[4],"nu_c":e[2]-e[4],
   "Hu":e[1]-e[2],"Hd":e[1]-e[3]}
checks=[]
def check(name,test):
    ok=bool(test);checks.append({"name":name,"passed":ok})
    assert ok,name
def vlist(v):return [str(x) for x in v]

ops={"up_Yukawa":q["Q"]+q["u_c"]+q["Hu"],
     "down_Yukawa":q["Q"]+q["d_c"]+q["Hd"],
     "lepton_Yukawa":q["L"]+q["e_c"]+q["Hd"],
     "neutrino_Dirac":q["L"]+q["nu_c"]+q["Hu"],
     "Hu_Hd":q["Hu"]+q["Hd"],
     "QQQL":3*q["Q"]+q["L"],
     "uc_uc_dc_ec":2*q["u_c"]+q["d_c"]+q["e_c"],
     "nu_c_nu_c":2*q["nu_c"]}
dress={"Hu_Hd":s.Matrix([-1,0]),"QQQL":s.Matrix([0,1]),"uc_uc_dc_ec":s.Matrix([2,-1])}
for name in ["up_Yukawa","down_Yukawa","lepton_Yukawa","neutrino_Dirac"]:
    check(name+" neutral",ops[name]==s.zeros(5,1))
for name,n in dress.items():
    # e^{i n.a} transforms with charge -K^T n.
    check(name+" axion dressing",ops[name]-K.T*n==s.zeros(5,1))
    check(name+" hypercharge zero",(ops[name].T*Y)[0]==0)
    check(name+" B-L zero",(ops[name].T*BL)[0]==0)
check("Hu hypercharge",(q["Hu"].T*Y)[0]==s.Rational(1,2))
check("Hd hypercharge",(q["Hd"].T*Y)[0]==-s.Rational(1,2))
check("Majorana hypercharge zero",(ops["nu_c_nu_c"].T*Y)[0]==0)
check("Majorana B-L is 2",(ops["nu_c_nu_c"].T*BL)[0]==2)
check("axions cannot cancel B-L",K*BL==s.zeros(2,1))
check("Majorana outside axion charge span",K.T.row_join(ops["nu_c_nu_c"]).rank()>K.rank())
# The homomorphisms on U(3)xU(2)xU(1)^3 use determinant characters.
check("color center divisibility",all(K[i,0]%3==0 for i in range(2)))
check("weak center divisibility",all(K[i,1]%2==0 for i in range(2)))
# Constant B-L flux conventions and a primitive gauge lattice are not new scalars.
# Gauge-invariant dressed Higgs bilinear has scalar dimension 2, coefficient mass^2.
result={"date":"2026-09-17","class":"EXACT charge algebra; no Wilson coefficient prediction",
    "conventions":"a -> a-K lambda; a has period 2pi; e^{i n.a} O has charge q_O-K^T n",
    "stack_order":["c","L","a","b","d"],"K":[vlist(K.row(i)) for i in range(2)],
    "field_charges":{k:vlist(v) for k,v in q.items()},
    "operators":{k:{"charge":vlist(v),"Y":str((v.T*Y)[0]),"B_minus_L":str((v.T*BL)[0]),
                   "integer_axion_exponent":vlist(dress[k]) if k in dress else None} for k,v in ops.items()},
    "characters":["(det U_L)^(-1) z_a z_b","(det U_c) (det U_L)^(-2) z_d"],
    "interpretation":{"Hu_Hd":"B_ij e^{-i a1} Hu_i dot Hd_j + h.c. allowed; B_ij unspecified",
        "QQQL":"e^{i a2} QQQL / Lambda^2 allowed as four-fermion operator if a coefficient is generated",
        "uc_uc_dc_ec":"e^{i(2a1-a2)} uc uc dc ec / Lambda^2 allowed if generated",
        "nu_c_nu_c":"not made invariant by these axions while B-L unbroken",
        "supersymmetric_variant":"if chiral T=s+i a exists, mu_ij e^{-T1} Hu_i dot Hd_j is gauge invariant; SUSY/Kahler/holomorphic coefficients not derived"},
    "checks":checks,"summary":{"passed":len(checks),"failed":0}}
(root/'operator_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result['summary']))
