"""Exact bridge between the v6 period-Manin basis and ``mtft.hecke``.

The period computation and the promoted Hecke engine chose different free
bases of the same 29-dimensional relative Manin-symbol quotient.  This module
reconstructs their change of basis from the *path endpoints*, using exact
continued fractions and the live Hecke model.  No frozen floating matrix is
used.

Conventions
-----------
R (29x29): columns are period generators expressed in mtft.hecke relative
coordinates.  det R = +1.

C (26x26): if B_P is the period cuspidal K basis and B_H is the promoted
``hecke.model()['K']`` basis, then B_H = B_P C.  det C = +1.

P = S^{-1} C maps promoted Hecke coordinates h to symplectic coordinates s:
    s = P h.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from functools import lru_cache
import json

from mtft import hecke as H
from .core import data_path, symplectic_change


def _basis_record():
    return json.loads(data_path("X0_143_period_basis_v6.json").read_text())


def _matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), Fr(0))
             for j in range(len(B[0]))] for i in range(len(A))]


def _transpose(A):
    return [list(r) for r in zip(*A)]


def _inverse(A):
    n = len(A)
    M = [[Fr(x) for x in row] + [Fr(1 if i == j else 0) for j in range(n)]
         for i, row in enumerate(A)]
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c]), None)
        if p is None:
            raise ValueError("singular matrix")
        M[c], M[p] = M[p], M[c]
        q = M[c][c]
        M[c] = [x / q for x in M[c]]
        for i in range(n):
            if i != c and M[i][c]:
                q = M[i][c]
                M[i] = [x - q*y for x, y in zip(M[i], M[c])]
    return [r[n:] for r in M]


def _pivots(A):
    M = [[Fr(x) for x in row] for row in A]
    nr, nc = len(M), len(M[0]); piv=[]; r=0
    for c in range(nc):
        p = next((i for i in range(r, nr) if M[i][c]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        q=M[r][c]; M[r]=[x/q for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c]:
                q=M[i][c]; M[i]=[x-q*y for x,y in zip(M[i],M[r])]
        piv.append(c);r+=1
        if r==nr: break
    return piv


def _det_bareiss(A):
    M=[[int(x) for x in r] for r in A]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        p=next((i for i in range(k,n) if M[i][k]),None)
        if p is None:return 0
        if p!=k:M[k],M[p]=M[p],M[k];sign=-sign
        pivot=M[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*pivot-M[i][k]*M[k][j])//prev
        prev=pivot
    return sign*M[-1][-1]


def _frac(pair):
    if pair is None:
        return None
    return Fr(int(pair[0]), int(pair[1]))


def _convergents(x: Fr):
    p,q=x.numerator,x.denominator; out=[]
    pm2,pm1,qm2,qm1=0,1,1,0; a,b=p,q
    while True:
        ak=a//b; pk,qk=ak*pm1+pm2,ak*qm1+qm2;out.append((pk,qk))
        if (pk,qk)==(p,q): return out
        a,b=b,a-ak*b;pm2,pm1,qm2,qm1=pm1,pk,qm1,qk


def _qcoords(m, v):
    ncols=len(m["cols"]); E=m["E"]
    x=[sum((Fr(m["Binv"][i][k])*v[k] for k in range(E)),Fr(0))
       for i in range(E)]
    return x[ncols:]


def _msym_edge(m,c,d):
    i=m["idx"][m["canon"](c % m["N"], d % m["N"])]
    return m["eid"][i],m["esign"][i]


def _zero_to(m, x: Fr):
    v=[Fr(0)]*m["E"]; cc=_convergents(x)
    for k in range(1,len(cc)):
        qk,qkm1=cc[k][1],cc[k-1][1]
        d=qkm1 if k%2 else -qkm1
        e,s=_msym_edge(m,qk,d);v[e]+=Fr(s)
    return v


def _zero_inf(m):
    v=[Fr(0)]*m["E"];e,s=_msym_edge(m,0,1);v[e]+=Fr(s);return v


def _inf_to(m,x:Fr):
    a=_zero_to(m,x);b=_zero_inf(m)
    return _qcoords(m,[a[i]-b[i] for i in range(m["E"])])


def _path_coords(m,a,b):
    if a is None:return _inf_to(m,b)
    if b is None:return [-x for x in _inf_to(m,a)]
    x,y=_inf_to(m,a),_inf_to(m,b)
    return [y[i]-x[i] for i in range(m["nq"])]


@lru_cache(maxsize=1)
def relative_basis_change():
    """R, 29x29 integer unimodular: period relative basis -> Hecke basis."""
    m=H.model(); rec=_basis_record(); cols=[]
    for p in rec["generator_paths"]:
        cols.append(_path_coords(m,_frac(p["start"]),_frac(p["end"])))
    R=[[cols[j][i] for j in range(29)] for i in range(29)]
    if any(x.denominator != 1 for r in R for x in r):
        raise ArithmeticError("period relative change is not integral")
    out=tuple(tuple(int(x) for x in r) for r in R)
    if _det_bareiss(out) != 1:
        raise ArithmeticError("period relative change is not unimodular")
    return out


@lru_cache(maxsize=1)
def cuspidal_basis_change():
    """C, 26x26 integer unimodular with B_H = B_period C."""
    m=H.model(); rec=_basis_record()
    R=[[Fr(x) for x in r] for r in relative_basis_change()]
    Kp=[[Fr(x) for x in r] for r in rec["K_period_29x26"]]
    PK=_matmul(R,Kp)
    KH=_transpose([[Fr(x) for x in r] for r in m["K"]])
    # independent rows of PK are pivots of PK^T
    rows=_pivots(_transpose(PK))
    if len(rows)!=26: raise ArithmeticError("period cusp basis rank != 26")
    A=[PK[i] for i in rows]; B=[KH[i] for i in rows]
    C=_matmul(_inverse(A),B)
    if _matmul(PK,C) != KH:
        raise ArithmeticError("period/Hecke cusp bridge residual nonzero")
    if any(x.denominator != 1 for r in C for x in r):
        raise ArithmeticError("period/Hecke cusp bridge is not integral")
    out=tuple(tuple(int(x) for x in r) for r in C)
    if _det_bareiss(out)!=1:
        raise ArithmeticError("period/Hecke cusp bridge is not unimodular")
    return out


@lru_cache(maxsize=1)
def hecke_to_symplectic_change():
    """P=S^-1 C, exact integer unimodular; symplectic coords s=P h."""
    S=[[Fr(x) for x in r] for r in symplectic_change()]
    C=[[Fr(x) for x in r] for r in cuspidal_basis_change()]
    P=_matmul(_inverse(S),C)
    if any(x.denominator!=1 for r in P for x in r):
        raise ArithmeticError("Hecke->symplectic map unexpectedly non-integral")
    return tuple(tuple(int(x) for x in r) for r in P)


__all__=["relative_basis_change","cuspidal_basis_change","hecke_to_symplectic_change"]
