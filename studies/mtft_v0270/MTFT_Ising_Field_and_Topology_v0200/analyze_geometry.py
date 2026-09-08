"""Fisher geometry of the exact (energy, magnetization) finite exponential family.

Natural coordinates are (beta, eta); sufficient statistics are
S=sum_edges spin_u*spin_v and M=sum_vertices spin_v. The metric is the full
sample Fisher metric Cov(S,M), without division by the number of spins.
Gaussian curvature K is used (scalar curvature is 2K).
"""
from __future__ import annotations
import json
from pathlib import Path
import itertools
import mpmath as mp
import numpy as np
from ising_extensions import joint_observables

ROOT=Path(__file__).resolve().parent


def metric_jets(points, degeneracy, beta, eta):
    logs=[beta*mp.mpf(x)+eta*mp.mpf(y) for x,y in points]
    mx=max(logs)
    weights=[mp.mpf(d)*mp.exp(v-mx) for d,v in zip(degeneracy,logs)]
    total=mp.fsum(weights)
    p=[v/total for v in weights]
    means=[mp.fsum(w*xy[a] for w,xy in zip(p,points)) for a in range(2)]
    centered=[(mp.mpf(x)-means[0],mp.mpf(y)-means[1]) for x,y in points]
    moment={}
    for order in (2,3,4):
        for j in range(order+1):
            moment[(order-j,j)]=mp.fsum(w*x**(order-j)*y**j for w,(x,y) in zip(p,centered))
    g=mp.matrix([[moment[(2-i-j,i+j)] for j in range(2)] for i in range(2)])
    third=lambda a,b,c:moment[(3-a-b-c,a+b+c)]
    fourth=lambda a,b,c,d:moment[(4-a-b-c-d,a+b+c+d)]-g[a,b]*g[c,d]-g[a,c]*g[b,d]-g[a,d]*g[b,c]
    return g,third,fourth


def three_curvatures(points, degeneracy, beta, eta, dps=60):
    from mtft.curvature import brioschi
    with mp.workdps(dps):
        beta,eta=mp.mpf(str(beta)),mp.mpf(str(eta))
        g,c3,c4=metric_jets(points,degeneracy,beta,eta)
        gi=g**-1;detg=mp.det(g)
        Kdet=-mp.det(mp.matrix([[g[0,0],g[0,1],g[1,1]],
                [c3(0,0,0),c3(0,0,1),c3(0,1,1)],
                [c3(0,0,1),c3(0,1,1),c3(1,1,1)]]))/(4*detg**2)
        Kbr=brioschi(g[0,0],g[0,1],g[1,1],c3(0,0,0),c3(0,0,1),c3(0,0,1),c3(0,1,1),c3(0,1,1),c3(1,1,1))
        Gamma=lambda a,b,c:mp.fsum(gi[a,d]*c3(d,b,c) for d in range(2))/2
        dg_inverse=[-gi*mp.matrix([[c3(i,j,l) for j in range(2)] for i in range(2)])*gi for l in range(2)]
        dGamma=lambda l,a,b,c:mp.fsum(dg_inverse[l][a,d]*c3(d,b,c)+gi[a,d]*c4(d,b,c,l) for d in range(2))/2
        def R(a):
            return dGamma(0,a,1,1)-dGamma(1,a,0,1)+mp.fsum(
                Gamma(a,0,m)*Gamma(m,1,1)-Gamma(a,1,m)*Gamma(m,0,1) for m in range(2))
        Kriemann=mp.fsum(g[0,a]*R(a) for a in range(2))/detg
        error=max(abs(Kdet-Kbr),abs(Kdet-Kriemann))
        assert error<mp.mpf('1e-45')
        return {'beta':str(beta),'eta':str(eta),'dps':dps,
                'K_determinant':mp.nstr(Kdet,50),'K_Brioschi':mp.nstr(Kbr,50),
                'K_Christoffel_Riemann':mp.nstr(Kriemann,50),
                'max_absolute_disagreement':mp.nstr(error,8),'fisher_determinant':mp.nstr(detg,35)}


def main():
    r=json.loads((ROOT/'joint_143.json').read_text())
    n,E,D=r['spins'],r['E'],r['counts']
    points=[];deg=[]
    for k,row in enumerate(D):
        for m,d in enumerate(row):
            if d:points.append((E-2*k,2*m-n));deg.append(d)
    categorical=three_curvatures([(0,0),(1,0),(0,1)],[1,1,1],0,0)
    assert abs(mp.mpf(categorical['K_determinant'])-mp.mpf('0.25'))<mp.mpf('1e-40')
    checks=[]
    for beta,eta in [(.3,0),(.64,0),(.64,.03),(1.,.08)]:
        c=three_curvatures(points,deg,beta,eta)
        double=joint_observables(D,n,E,beta,eta,higher=True)['gaussian_curvature']
        c['double_absolute_error']=abs(float(c['K_determinant'])-double)
        assert c['double_absolute_error']<1e-9
        checks.append(c)
        print('geometry',beta,eta,c['K_determinant'],c['max_absolute_disagreement'],flush=True)
    beta_grid=np.linspace(.05,1.,96)
    eta_grid=np.linspace(-.08,.08,81)
    names=['magnetization_per_spin','absolute_magnetization_per_spin','susceptibility_per_spin',
           'interaction_energy_per_spin','gaussian_curvature','fisher_determinant']
    grids={name:np.zeros((len(beta_grid),len(eta_grid))) for name in names}
    for i,beta in enumerate(beta_grid):
        for j,eta in enumerate(eta_grid):
            o=joint_observables(D,n,E,beta,eta,higher=True)
            for name in names:grids[name][i,j]=o[name]
    assert np.all(grids['fisher_determinant']>0)
    assert np.all(np.isfinite(grids['gaussian_curvature']))
    m_error=float(np.max(abs(grids['magnetization_per_spin']+grids['magnetization_per_spin'][:,::-1])))
    K_error=float(np.max(abs(grids['gaussian_curvature']-grids['gaussian_curvature'][:,::-1])))
    assert m_error<1e-10 and K_error<1e-9
    np.savez_compressed(ROOT/'field_geometry_grid.npz',beta=beta_grid,eta=eta_grid,**grids)
    summary={'coordinates':'Natural parameters beta and eta=beta*h; full Fisher metric, not divided by n.',
             'curvature_convention':'Gaussian K; scalar R=2K; categorical simplex gate K=+1/4.',
             'claim_class':'Numerical/DIAGNOSTIC finite-ensemble geometry; selected points cross-checked at 60 digits via three curvature assemblies.',
             'categorical_gate':categorical,'cross_checks':checks,
             'grid_shape':[len(beta_grid),len(eta_grid)],'beta_domain':[float(beta_grid[0]),float(beta_grid[-1])],
             'eta_domain':[float(eta_grid[0]),float(eta_grid[-1])],
             'grid_K_range':[float(grids['gaussian_curvature'].min()),float(grids['gaussian_curvature'].max())],
             'spin_flip_symmetry_errors':{'magnetization':m_error,'curvature':K_error}}
    (ROOT/'geometry_checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('geometry grid',summary['grid_K_range'],flush=True)


if __name__=='__main__':main()
