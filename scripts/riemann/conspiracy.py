"""Multi-zero conspiracy test for L2: do several off-line quadruplets cancel?
Case A: two quadruplets, same beta=0.75, at gamma_1 and gamma_2      -> slope should stay -0.25
Case B: four quadruplets, beta=0.75, at gamma_1..gamma_4             -> slope should stay -0.25
Case C: mixed beta=0.75 (gamma_1) + beta=0.9 (gamma_3)               -> largest beta wins: slope -0.4
Case D: adversarial: beta=0.75 at gamma_1 with COEFFICIENT NEGATED on one member
        (breaks the canonical coefficient formula -- tests whether sign flips can cancel;
         not realizable by an actual L-function, pure stress test)
"""
import mpmath as mp, json, time, math
t0=time.time(); mp.mp.dps=30
gammas=[mp.mpf(s) for s in json.load(open("zeros.json"))]
def G(w): return -mp.zeta(w,derivative=1)/mp.zeta(w)
MAIN=[(2*G(2),mp.mpf(3)),(mp.zeta(0),mp.mpf(2))]
def coeff(rho): return -mp.gamma(rho+1)*mp.zeta(rho-1)
def kappa(X,zs):   # zs = list of (rho, coefficient)
    S=DS=D2S=mp.mpf(0)
    for c,a in MAIN:
        T=c*X**(-a); S+=T; DS+=-a*T; D2S+=a*a*T
    for rho,cf in zs:
        a=rho+1; T=cf*X**(-a)
        S+=2*mp.re(T); DS+=2*mp.re(-a*T); D2S+=2*mp.re(a*a*T)
    return (D2S*S-DS*DS)/(S*S)
def slopeD(zs, note):
    ys=[mp.mpf(10)**(mp.mpf(-7)+mp.mpf(5.2)*i/90) for i in range(91)]
    bins={}
    for y in ys:
        X=2*mp.pi*y
        dk=kappa(X,zs)-kappa(X,[])
        b=math.floor(mp.log10(y)*2)/2
        bins.setdefault(b,[]).append(abs(dk))
    pts=sorted((b,mp.log10(mp.sqrt(sum(v*v for v in vals)/len(vals)))) for b,vals in bins.items())
    n=len(pts); sx=sum(p[0] for p in pts); sy=sum(p[1] for p in pts)
    sxx=sum(p[0]*p[0] for p in pts); sxy=sum(p[0]*p[1] for p in pts)
    sl=(n*sxy-sx*sy)/(n*sxx-sx*sx) - mp.mpf(3)/2
    print(f"  {note}: slope(log|D|) = {mp.nstr(sl,4)}")
    return sl

half=mp.mpf(1)/2
on=[(mp.mpc(half,g),coeff(mp.mpc(half,g))) for g in gammas[:12]]
def quad(beta,g): return [(mp.mpc(beta,g),coeff(mp.mpc(beta,g))),(mp.mpc(1-beta,g),coeff(mp.mpc(1-beta,g)))]

b=mp.mpf("0.75")
print("multi-zero conspiracy tests (prediction: slope = 1/2 - max beta):")
slopeD(quad(b,gammas[0])+quad(b,gammas[1])+on[2:], "A: 2x beta=0.75            (pred -0.25)")
slopeD(quad(b,gammas[0])+quad(b,gammas[1])+quad(b,gammas[2])+quad(b,gammas[3])+on[4:], "B: 4x beta=0.75            (pred -0.25)")
slopeD(quad(b,gammas[0])+on[1:2]+quad(mp.mpf("0.9"),gammas[2])+on[3:], "C: 0.75@g1 + 0.90@g3       (pred -0.40)")
zsD=[(mp.mpc(b,gammas[0]),coeff(mp.mpc(b,gammas[0]))),(mp.mpc(1-b,gammas[0]),-coeff(mp.mpc(1-b,gammas[0])))]+on[1:]
slopeD(zsD, "D: adversarial sign flip   (pred -0.25 still: same-gamma partners cannot cancel across X)")
print(f"({time.time()-t0:.1f}s)")
