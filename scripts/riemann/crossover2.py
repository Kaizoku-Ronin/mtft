"""Stable deep scan.  Delta kappa = D[P/Q], P = DZ*M - Z*DM, DP = D2Z*M - Z*D2M,
Q = M(M+Z), DQ = DM*(M+Z) + M*(DM+DZ);  Delta kappa = (DP*Q - P*DQ)/Q^2.
Only O(1)-relative cancellations -> works at any depth with modest dps."""
import mpmath as mp, json, time, math
t0=time.time(); mp.mp.dps=60
gammas=[mp.mpf(s) for s in json.load(open("zeros.json"))]
def G(w): return -mp.zeta(w,derivative=1)/mp.zeta(w)
MAIN=[(2*G(2),mp.mpf(3)),(mp.zeta(0),mp.mpf(2))]
def coeff(rho): return -mp.gamma(rho+1)*mp.zeta(rho-1)

def sums(X, zs):
    M=DM=D2M=mp.mpf(0)
    for c,a in MAIN:
        T=c*X**(-a); M+=T; DM+=-a*T; D2M+=a*a*T
    Z=DZ=D2Z=mp.mpf(0)
    for rho,cf in zs:
        a=rho+1; T=cf*X**(-a)
        Z+=2*mp.re(T); DZ+=2*mp.re(-a*T); D2Z+=2*mp.re(a*a*T)
    return M,DM,D2M,Z,DZ,D2Z

def dkappa(X, zs):
    M,DM,D2M,Z,DZ,D2Z = sums(X,zs)
    P  = DZ*M - Z*DM
    DP = D2Z*M - Z*D2M
    Q  = M*(M+Z)
    DQ = DM*(M+Z) + M*(DM+DZ)
    return (DP*Q - P*DQ)/(Q*Q)

def slope_window(zs, lo, hi, note):
    ys=[mp.mpf(10)**(mp.mpf(lo)+mp.mpf(hi-lo)*i/90) for i in range(91)]
    bins={}
    for y in ys:
        X=2*mp.pi*y
        dk=dkappa(X,zs)
        b=math.floor(float(mp.log10(y))*2)/2
        bins.setdefault(b,[]).append(abs(dk))
    pts=sorted((b,mp.log10(mp.sqrt(sum(v*v for v in vals)/len(vals)))) for b,vals in bins.items())
    n=len(pts); sx=sum(p[0] for p in pts); sy=sum(p[1] for p in pts)
    sxx=sum(p[0]*p[0] for p in pts); sxy=sum(p[0]*p[1] for p in pts)
    sl=(n*sxy-sx*sy)/(n*sxx-sx*sx) - mp.mpf(3)/2
    print(f"  y in [1e{lo}, 1e{hi}]  {note}: slope(log|D|) = {mp.nstr(sl,4)}")

half=mp.mpf(1)/2
on=[(mp.mpc(half,g),coeff(mp.mpc(half,g))) for g in gammas[:12]]
def quad(beta,g): return [(mp.mpc(beta,g),coeff(mp.mpc(beta,g))),(mp.mpc(1-beta,g),coeff(mp.mpc(1-beta,g)))]
zsC=quad(mp.mpf("0.75"),gammas[0])+on[1:2]+quad(mp.mpf("0.9"),gammas[2])+on[3:]

# sanity: stable formula reproduces the B2 single-quadruplet slope
zsB=quad(mp.mpf("0.75"),gammas[0])+on[1:]
print("sanity (stable formula vs stage2): single 0.75 quadruplet, [1e-7,1.6e-2]:")
slope_window(zsB,-7,mp.log10(mp.mpf("0.016")), "pred -0.25")

print("\nCase C crossover scan (pred -0.25 above ~1e-50, -0.40 below):")
for lo,hi,note in [(-7,-2,"shallow "),(-30,-20,"mid     "),(-48,-40,"near    "),(-70,-58,"deep    "),(-95,-80,"deepest ")]:
    slope_window(zsC,lo,hi,note)
print(f"({time.time()-t0:.1f}s)")
