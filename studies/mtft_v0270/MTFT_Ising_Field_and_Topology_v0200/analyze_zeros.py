"""Exact polynomial structure, Sturm certificates for Lee-Yang zeros, Fisher diagnostics."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import time
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parent


def lee_yang():
    row=json.loads((ROOT/'joint_143.json').read_text())
    n,E,D=row['spins'],row['E'],row['counts']
    y,u=sp.symbols('y u')
    out=[]
    for a in (2,4,8):
        started=time.perf_counter()
        coeff=[sum(D[k][m]*a**(E-k) for k in range(E+1)) for m in range(n+1)]
        assert coeff==coeff[::-1]
        # P(y)/y^(n/2) = c_mid + sum_j c_(mid+j) C_j(u),
        # C_j(y+y^-1)=y^j+y^-j, with C_0=2 and C_1=u.
        mid=n//2
        prev,cur=sp.Poly(2,u),sp.Poly(u,u)
        reduced=sp.Poly(coeff[mid],u)
        for j in range(1,mid+1):
            reduced+=coeff[mid+j]*cur
            prev,cur=cur,sp.Poly(u,u)*cur-prev
        reduced=reduced.primitive()[1]
        field=sp.Poly.from_list(list(reversed(coeff)),y)
        substituted=sp.Poly(sp.expand(y**mid*reduced.as_expr().subs(u,y+1/y)),y)
        assert field.primitive()[1]==substituted.primitive()[1]
        intervals=sp.polys.polytools.intervals(reduced,eps=sp.Rational(1,10**18))
        assert sum(mult for _,mult in intervals)==mid
        assert all(-2<lo<=hi<2 and mult==1 for (lo,hi),mult in intervals)
        roots=[]
        for (lo,hi),mult in intervals:
            with mp.workdps(50):
                midpoint=mp.mpf(str(sp.N((lo+hi)/2,45)))
                theta=mp.acos(midpoint/2)
                roots.append({'u_interval':[str(lo),str(hi)],'multiplicity':mult,
                              'angle_radians':float(theta),'re_y':float(mp.cos(theta)),
                              'im_y_positive':float(mp.sin(theta))})
        out.append({'satisfied_weight':a,'unsatisfied_weight':1,'beta':math.log(a)/2,
                    'field_coefficients_ascending':[str(c) for c in coeff],
                    'reduced_u_coefficients_descending':[str(c) for c in reduced.all_coeffs()],
                    'gates':{'reciprocal_polynomial':True,'substitution_identity_exact':True,
                             'all_28_u_roots_real_simple_in_open_minus2_plus2':True,
                             'all_56_y_roots_on_unit_circle':True},
                    'positive_imaginary_roots':roots,'closest_angle':min(r['angle_radians'] for r in roots),
                    'seconds':time.perf_counter()-started})
        (ROOT/'lee_yang.json').write_text(json.dumps({'variable':'y=exp(2*eta), eta=beta*h. All field weights positive at real eta.',
          'certificate':'Exact reciprocal reduction and rational Sturm root isolation; complex root coordinates are numerical displays of algebraically certified unit-circle roots.',
          'rows':out},indent=2)+'\n')
        print('Lee-Yang',a,'angle edge',out[-1]['closest_angle'],'seconds',out[-1]['seconds'],flush=True)


def fisher():
    baseline=json.loads((ROOT/'baseline_results.json').read_text())
    row=next(r for r in baseline['levels'] if r['N']==143)
    x=sp.Symbol('x')
    D=sp.Poly.from_list(list(reversed(row['density_of_states'])),x)
    prefactor,factors=sp.factor_list(D)
    Q=next(f for f,m in factors if f.degree()==45)
    known=sp.Poly(2*(1+x)**28*(1+x*x),x)
    assert known*Q==D
    started=time.perf_counter()
    roots40=Q.nroots(n=40,maxsteps=500)
    roots65=Q.nroots(n=65,maxsteps=700)
    roots=[]
    coeff=[mp.mpf(int(c)) for c in Q.all_coeffs()]
    max_res=mp.mpf(0);max_shift=mp.mpf(0)
    with mp.workdps(80):
        coeff=[mp.mpf(int(c)) for c in Q.all_coeffs()]
        r40=[mp.mpc(str(sp.N(sp.re(z),42)),str(sp.N(sp.im(z),42))) for z in roots40]
        for z in roots65:
            r=mp.mpc(str(sp.N(sp.re(z),70)),str(sp.N(sp.im(z),70)))
            scale=mp.fsum(abs(c)*abs(r)**(len(coeff)-1-i) for i,c in enumerate(coeff))
            residual=abs(mp.polyval(coeff,r))/scale
            shift=min(abs(r-v) for v in r40)
            max_res=max(max_res,residual);max_shift=max(max_shift,shift)
            roots.append({'real':mp.nstr(r.real,65),'imag':mp.nstr(r.imag,65),
                          'normalized_residual':mp.nstr(residual,6)})
        assert max_res<mp.mpf('1e-55') and max_shift<mp.mpf('1e-30')
    out={'class':'Factorization EXACT; degree-45 roots DIAGNOSTIC with two precisions and normalized residual checks. No certified complex disks.',
         'variable':'x=exp(-2*beta), temperature variable; no real positive zero of the finite partition function.',
         'factorization':'2*(1+x)^28*(1+x^2)*Q45(x)',
         'Q45_coefficients_descending':[str(c) for c in Q.all_coeffs()],
         'known_roots':[{'real':-1,'imag':0,'multiplicity':28},{'real':0,'imag':1,'multiplicity':1},{'real':0,'imag':-1,'multiplicity':1}],
         'Q45_roots':roots,'max_normalized_residual':str(max_res),'max_40_vs_65_digit_root_shift':str(max_shift),
         'seconds':time.perf_counter()-started}
    (ROOT/'fisher_zeros.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Fisher roots',len(roots),'seconds',out['seconds'],'residual',str(max_res),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--which',choices=['lee','fisher','both'],default='both')
    a=p.parse_args()
    if a.which in ('lee','both'):lee_yang()
    if a.which in ('fisher','both'):fisher()
