#!/usr/bin/env python3
"""Exact polynomial reductions, with no specialization of the parameter c."""
from fractions import Fraction as Q
from pathlib import Path
import json

# Sparse bivariate polynomials: (degree in x, degree in c) -> rational.
def clean(p):return {k:v for k,v in p.items() if v}
def add(a,b):
    z=dict(a)
    for k,v in b.items():z[k]=z.get(k,0)+v
    return clean(z)
def scale(a,s):return clean({k:v*s for k,v in a.items()})
def mul(a,b):
    z={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():z[i+k,j+l]=z.get((i+k,j+l),0)+v*w
    return clean(z)
def dx(a):return {(i-1,j):i*v for (i,j),v in a.items() if i}
def exact_term(r,p):return add(mul(dx(r),p),scale(mul(r,dx(p)),Q(-1,2)))

def verify():
    P={(4,0):1,(2,1):2,(0,2):1,(0,1):1}
    N={(2,0):-1,(0,1):-1,(0,0):Q(-1,2)}
    denom0={(0,2):4,(0,1):4}
    R0={(3,0):-1,(1,1):-3,(1,0):-1}
    coeff0={(0,1):-1,(0,0):-1,(2,0):1}
    assert mul(denom0,N)==add(mul(coeff0,P),exact_term(R0,P))
    denom2={(0,1):4}
    R2={(1,1):1,(3,0):-1}
    coeff2={(0,1):-1,(2,0):1}
    assert mul(denom2,mul({(2,0):1},N))==add(mul(coeff2,P),exact_term(R2,P))
    return {'status':'EXACT','domain':'Q[c,x] after clearing denominators; c != 0,-1',
        'parameter_specializations_used':False,
        'verified_identities':2,
        'connection':[['-1/(4c)','1/(4c(c+1))'],['-1/4','1/(4c)']],
        'exact_differential_R0':'-x(x^2+3c+1)/(4c(c+1))',
        'exact_differential_R2':'x/4-x^3/(4c)',
        'scalar_equation':'I"+(1/c+1/(c+1)) I\'+(4c-1)/(16c^2(c+1)) I=0'}

if __name__=='__main__':
    result=verify()
    (Path(__file__).resolve().parent/'elliptic_exact_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
