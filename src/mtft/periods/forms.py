"""Native q-expansion evaluation and Bergman density for X_0(143).

The 13 raw rows are pinned to the v6 period ordering:
  f1 | four f2 embeddings | six f3 embeddings | 11a(q) | 11a(q^13).
The nontrivial newform coefficients come from exact power-basis a_n data in
``mtft.codifferent``; only the real embeddings are numerical.
"""
from __future__ import annotations

from functools import lru_cache
import json, math
import mpmath as mp

from mtft.codifferent import ORBITS, eigen_an
from .core import GENUS, data_path, omega_symplectic, riemann_matrix

CURVE_143A1=(0,-1,1,-1,-2)
CURVE_11A1=(0,-1,1,-10,-20)


def _basis_record():
    return json.loads(data_path("X0_143_period_basis_v6.json").read_text())


def _isprime(n):
    return n>=2 and all(n%d for d in range(2,int(n**0.5)+1))


def _curve_ap(curve,p):
    a1,a2,a3,a4,a6=curve; N=1
    for x in range(p):
        rhs=(x**3+a2*x*x+a4*x+a6)%p
        for y in range(p):
            if (y*y+a1*x*y+a3*y-rhs)%p==0:N+=1
    return p+1-N


def _elliptic_an(curve, conductor, nmax):
    ap={p:_curve_ap(curve,p) for p in range(2,nmax+1) if _isprime(p)}
    pp={}
    for p,a in ap.items():
        vals=[1,a];q=p*p
        r=2
        while q<=nmax:
            vals.append(a*vals[-1] if conductor%p==0
                        else a*vals[-1]-p*vals[-2])
            q*=p;r+=1
        pp[p]=vals
    out=[0]*(nmax+1);out[1]=1
    for n in range(2,nmax+1):
        x=n;v=1
        for p in pp:
            if x%p:continue
            e=0
            while x%p==0:x//=p;e+=1
            v*=pp[p][e]
        out[n]=v
    return out


def _peval(coords,x):
    s=mp.mpf('0')
    for c in reversed(coords):s=s*x+c
    return s


def _ordered_roots(orbit,dps):
    info=ORBITS[orbit]; poly=list(reversed(info["poly_low"]))
    roots=mp.polyroots(poly,maxsteps=500,error=False)
    roots=[mp.re(r) if abs(mp.im(r))<mp.mpf(10)**(-(dps//2)) else r for r in roots]
    a2=eigen_an(orbit)[2]
    targets=_basis_record()["a2_per_row"]
    targets=targets[1:5] if orbit=="f2" else targets[5:11]
    targets=[mp.mpf(t) for t in targets]
    vals=[_peval(a2,r) for r in roots]
    ordered=[]; used=set()
    for t in targets:
        j=min((j for j in range(len(roots)) if j not in used),key=lambda j:abs(vals[j]-t))
        used.add(j);ordered.append(roots[j])
    return ordered


@lru_cache(maxsize=8)
def _raw_tuple(nmax=140,dps=60):
    with mp.workdps(dps):
        f1=_elliptic_an(CURVE_143A1,143,nmax)
        old=_elliptic_an(CURVE_11A1,11,nmax)
        rows=[[mp.mpf(x) for x in f1]]
        for orbit in ("f2","f3"):
            an=eigen_an(orbit); roots=_ordered_roots(orbit,dps)
            for root in roots:
                rows.append([_peval(an[n],root) for n in range(nmax+1)])
        rows.append([mp.mpf(x) for x in old])
        rows.append([mp.mpf(0) if n%13 else mp.mpf(old[n//13])
                     for n in range(nmax+1)])
        return tuple(tuple(+x for x in row) for row in rows)


def raw_qexpansions(nmax=140,dps=60):
    """13 x (nmax+1) raw q-expansion coefficients in v6 row order."""
    return mp.matrix(_raw_tuple(nmax,dps))


def q_tail_bound(y,nmax=140):
    r=mp.e**(-2*mp.pi*mp.mpf(y)); M=mp.mpf(nmax)
    rho=r*((M+2)/(M+1))**mp.mpf('1.5')
    if rho>=1:return mp.inf
    return (M+1)**mp.mpf('1.5')*r**(M+1)/(1-rho)


def raw_form_values(z,nmax=140,dps=60):
    """Values of the 13 raw holomorphic forms at z via q-series."""
    with mp.workdps(dps):
        z=mp.mpc(z);q=mp.e**(2j*mp.pi*z);A=raw_qexpansions(nmax,dps)
        return mp.matrix([sum(A[i,n]*q**n for n in range(1,nmax+1))
                          for i in range(GENUS)])


def normalized_form_values(z,nmax=140,dps=60):
    """Alpha-normalized holomorphic differential values at z."""
    with mp.workdps(dps):
        O=omega_symplectic(dps); A=O[:,:GENUS]
        return +(A**-1 * raw_form_values(z,nmax,dps))


def bergman_density(z,nmax=140,dps=60):
    r"""Canonical Bergman scalar density w^* (Im tau)^-1 w.

    This is the coefficient density in the local coordinate dz, not a claim
    that the upper-half-plane coordinate is a global physical position.
    """
    with mp.workdps(dps):
        w=normalized_form_values(z,nmax,dps);t=riemann_matrix(dps)
        Y=mp.matrix([[mp.im(t[i,j]) for j in range(GENUS)] for i in range(GENUS)])
        return +mp.re((w.conjugate().T*(Y**-1)*w)[0])


__all__=["CURVE_143A1","CURVE_11A1","raw_qexpansions","q_tail_bound",
         "raw_form_values","normalized_form_values","bergman_density"]
