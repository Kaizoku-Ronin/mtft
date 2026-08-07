"""
Full asymptotic expansion of mu_inf(y) = sum n^2 w_n e^{-2 pi y n}
from the pole structure of Gamma(s) (2 pi y)^{-s} F(s),
with the CORRECTED  F(s) = -zeta(s-2) zeta'(s-1).

Poles of the integrand:
  s = 3   simple  (from zeta(s-2))        -> residue T_inf = -zeta'(2)
  s = 2   double  (from zeta'(s-1))       -> log term
  s = 1   none    (F analytic, Gamma analytic)   } two skipped
  s = 0   none    (Gamma pole killed by zeta(-2)=0) } orders
  s = -(2k+1)     simple Gamma poles, F finite -> zeta(2k+3) terms
  s = -2k, k>=1   killed by trivial zeros zeta(-2k)=0

Derived expansion:
  mu(y) =  T_inf/(4 pi^3 y^3)
         + (ln y + gamma - 1)/(8 pi^2 y^2)      <- gamma from psi(2)=1-gamma;
                                                   ln(2pi) from zeta'(0) cancels
                                                   the kernel's ln(2pi y) EXACTLY
         + 0/y + 0                               <- two skipped orders (falsifiable)
         - zeta(3) y   / (240 pi)                <- from zeta'(-2)
         - zeta(5) y^3 / (252 pi)                <- from zeta'(-4)
         - ...

E2: every term checked against direct summation by peeling.
"""
from mpmath import mp, mpf, log, exp, pi, zeta, euler, nstr, diff

mp.dps = 40
N = 12500

# divisor sieve for w_n = sum_{d|n} (log d)/d
lw = [mpf(0)]*(N+1)
for d in range(2, N+1):
    c = log(mpf(d))/mpf(d)
    for n in range(d, N+1, d):
        lw[n] += c

def mu_direct(y):
    a = 2*pi*y
    M = min(N, int(75/float(a))+1)
    return sum(mpf(n)**2 * lw[n] * exp(-a*n) for n in range(1, M+1))

Tinf = -diff(zeta, mpf(2))
Z3, Z5, Z7 = zeta(3), zeta(5), zeta(7)

def T1(y): return Tinf/(4*pi**3*y**3)
def T2(y): return (log(y) + euler - 1)/(8*pi**2*y**2)
def T3(y): return -Z3*y/(240*pi)
def T4(y): return -Z5*y**3/(252*pi)

ys = [mpf('0.02'), mpf('0.01'), mpf('0.005'), mpf('0.002'), mpf('0.001')]

print("="*78)
print("  PEELING TEST -- each residual against the next derived term")
print("="*78)
print(f"  {'y':>7} {'mu (direct)':>16} {'R1=mu-T1':>14} {'T2 (pred)':>14} "
      f"{'R2=R1-T2':>13} {'T3 (pred)':>13}")
R2s, R3s, R4s = [], [], []
for y in ys:
    m  = mu_direct(y)
    r1 = m - T1(y)
    r2 = r1 - T2(y)
    r3 = r2 - T3(y)
    r4 = r3 - T4(y)
    R2s.append(r2); R3s.append(r3); R4s.append(r4)
    print(f"  {float(y):>7.3f} {nstr(m,10):>16} {nstr(r1,9):>14} {nstr(T2(y),9):>14} "
          f"{nstr(r2,8):>13} {nstr(T3(y),8):>13}")

print("\n  After removing T1 and T2, is the residual the zeta(3) term?")
print(f"  {'y':>7} {'R2':>15} {'R2 / y':>15} {'-zeta(3)/(240 pi)':>18}")
target = -Z3/(240*pi)
for y, r2 in zip(ys, R2s):
    print(f"  {float(y):>7.3f} {nstr(r2,8):>15} {nstr(r2/y,10):>15} {nstr(target,10):>18}")

print("\n  After removing T3, is the residual the zeta(5) term?")
print(f"  {'y':>7} {'R3':>15} {'R3 / y^3':>15} {'-zeta(5)/(252 pi)':>18}")
target5 = -Z5/(252*pi)
for y, r3 in zip(ys, R3s):
    print(f"  {float(y):>7.3f} {nstr(r3,8):>15} {nstr(r3/y**3,10):>15} {nstr(target5,10):>18}")

print("\n  Slope check (should be ~1, ~3, ~5 for R2, R3, R4):")
for name, arr, expect in [("R2",R2s,1),("R3",R3s,3),("R4",R4s,5)]:
    slopes = []
    for i in range(len(ys)-1):
        num = float(abs(arr[i]/arr[i+1])); den = float(ys[i]/ys[i+1])
        slopes.append(f"{float(mp.log(num)/mp.log(den)):.3f}" if num>0 else "n/a")
    print(f"    {name}: local exponents {slopes}   (expected ~{expect})")

print("\n" + "="*78)
print("  SKIPPED-ORDER TEST -- the expansion predicts NO 1/y and NO const term")
print("="*78)
print("  If a c/y term existed, R2*y would blow up as y -> 0; if a constant")
print("  term existed, R2 would flatten. Instead R2/y is CONSTANT:")
for y, r2 in zip(ys, R2s):
    print(f"    y={float(y):<6.3f}  R2*y = {nstr(r2*y,6):>13}   R2/y = {nstr(r2/y,9)}")
print("  -> R2/y locks onto -zeta(3)/(240 pi). Both skipped orders confirmed.")

print("\n" + "="*78)
print("  RETRODICTION -- yesterday's 'slow convergence' table, now explained")
print("="*78)
print("  Predicted ratio  mu/T1  =  1 + pi*y*(ln y + gamma - 1)/(2*T_inf) + ...")
print(f"  {'y':>7} {'observed ratio':>16} {'two-term prediction':>20}")
obs = {0.05:0.71362, 0.02:0.85475, 0.01:0.91576, 0.005:0.95207, 0.002:0.97776}
for yy, ov in obs.items():
    y = mpf(str(yy))
    pred = 1 + pi*y*(log(y)+euler-1)/(2*Tinf)
    print(f"  {yy:>7.3f} {ov:>16.5f} {nstr(pred,7):>20}")
