"""
NULL MODEL TEST for the AG Dissection constant scan.

E2 rule: check the quantity by a route that does not pass through
the formula being checked.

Method: run the IDENTICAL search machinery against DECOY targets --
numbers with no theoretical meaning. If the machinery produces just as
many "discoveries" for meaningless targets, the machinery is the source
of the hits, not the mathematics.
"""
import math, random
from itertools import combinations

random.seed(20260406)

# ---- exact same constant list used in the session scan ----
C = {
 'pi':3.14159265358979324,'e':2.71828182845904524,'gamma':0.57721566490153286,
 'delta':4.66920160910299067,'alpha_F':2.50290787509589282,
 'zeta2':math.pi**2/6,'zeta3':1.20205690315959429,'zeta4':math.pi**4/90,
 'zeta5':1.03692775514336993,'Tinf':0.93754825431584375,'M':0.26149721284764278,
 'Omega':0.56714329040978387,'Catalan':0.91596559417721902,
 'ln2':math.log(2),'ln3':math.log(3),'phi':(1+math.sqrt(5))/2,
 'sqrt2':math.sqrt(2),'sqrt3':math.sqrt(3),'Khinchin':2.68545200106530645,
 'Glaisher':1.28242712910062264,'LandauRam':0.76422365358922066,
 'Artin':0.37395581361920229,'TwinC2':0.66016181584686957,
 'Plastic':1.32471795724474603,'Viswanath':1.13198824879045735,
 'Levy':1.18656911041562545,
}
names = list(C); vals = [C[k] for k in names]; n = len(names)

def all_expressions():
    """Every 2-constant expression the session scan actually evaluated."""
    out = []
    for i, j in combinations(range(n), 2):
        a, b = vals[i], vals[j]
        out += [a/b, b/a, a+b, a-b, b-a, a*b,
                a*a+b, a+b*b, a*a*b, a*b*b, a*a-b]
    return [v for v in out if math.isfinite(v)]

def triple_expressions():
    out = []
    idx = list(range(min(10, n)))
    for i, j, k in combinations(idx, 3):
        a, b, c = vals[i], vals[j], vals[k]
        out += [a*b*c, a*b+c, a*b-c, a+b+c, (a+b)*c, a/(b*c), (a*b)/c]
    return [v for v in out if math.isfinite(v)]

EXPR = all_expressions() + triple_expressions()
print(f"Expressions evaluated by the scan machinery: {len(EXPR):,}")

def count_hits(targets, tol):
    """How many expressions land within `tol` relative error of some target."""
    h = 0
    for t in targets:
        if t == 0: continue
        lo, hi = abs(t)*(1-tol), abs(t)*(1+tol)
        h += sum(1 for v in EXPR if lo <= abs(v) <= hi)
    return h

# ---- the REAL targets used in the session ----
REAL = {
 'y_c':0.18213038440927,'3/13':3/13,'2/3':2/3,'12':12.0,'2':2.0,
 'gamma':0.57721566490153286,'Tinf':0.93754825431584375,
 'e^gamma':math.exp(0.57721566490153286),'1/2':0.5,'1':1.0,
 '13':13.0,'11':11.0,'alpha_inv':137.035999084,
}

# ---- DECOY targets: same numeric ranges, zero theoretical meaning ----
DECOY = {
 'd1':0.19372,'d2':0.25714,'d3':0.71236,'d4':11.0,'d5':3.0,
 'd6':0.61934,'d7':0.88171,'d8':1.69422,'d9':0.47,'d10':1.13,
 'd11':14.0,'d12':9.0,'d13':128.44,
}

print("\n" + "="*72)
print("  TEST 1 -- Real targets vs decoy targets, identical machinery")
print("="*72)
print(f"  {'tolerance':>10} {'REAL hits':>11} {'DECOY hits':>11}   verdict")
print(f"  {'-'*10} {'-'*11} {'-'*11}   {'-'*30}")
for tol in [1e-2, 5e-3, 1e-3, 5e-4, 1e-4]:
    r = count_hits(REAL.values(), tol)
    d = count_hits(DECOY.values(), tol)
    verdict = "indistinguishable" if 0.5 <= (r/max(d,1)) <= 2.0 else "DIFFERENT"
    print(f"  {tol:>10.0e} {r:>11,} {d:>11,}   {verdict}")

# ---- Test 2: distribution over many random target sets ----
print("\n" + "="*72)
print("  TEST 2 -- Distribution of hit counts over 400 random target sets")
print("="*72)
def random_targets(k=13):
    out = []
    for _ in range(k):
        r = random.random()
        if r < 0.55: out.append(random.uniform(0.15, 1.0))
        elif r < 0.8: out.append(random.uniform(1.0, 15.0))
        else: out.append(random.uniform(100.0, 150.0))
    return out

for tol in [1e-3, 1e-4]:
    counts = sorted(count_hits(random_targets(), tol) for _ in range(400))
    real_ct = count_hits(REAL.values(), tol)
    pct = 100.0 * sum(1 for c in counts if c >= real_ct) / len(counts)
    print(f"\n  tolerance {tol:.0e}")
    print(f"    random target sets: median {counts[200]}, "
          f"5th pct {counts[20]}, 95th pct {counts[380]}")
    print(f"    REAL target set:    {real_ct}")
    print(f"    -> {pct:.1f}% of random target sets score at least as high")
    print(f"    -> {'NOT significant' if pct > 5 else 'significant'}")

# ---- Test 3: the continued-fraction claim under Gauss-Kuzmin ----
print("\n" + "="*72)
print("  TEST 3 -- 'CF(Tinf - M) contains both 11 and 13' under Gauss-Kuzmin")
print("="*72)
def p_gk(k):  # P(partial quotient = k)
    return math.log2(1 + 1.0/(k*(k+2)))
p11, p13, T = p_gk(11), p_gk(13), 15
q11 = 1-(1-p11)**T; q13 = 1-(1-p13)**T
print(f"  P(a_i = 11) = {p11:.5f}   P(a_i = 13) = {p13:.5f}")
print(f"  In {T} terms: P(>=1 eleven) = {q11:.4f}, P(>=1 thirteen) = {q13:.4f}")
print(f"  P(a given constant shows BOTH) ~ {q11*q13:.4f}  ({q11*q13*100:.1f}%)")
n_checked = 11
print(f"  Constants whose CF was inspected: {n_checked}")
print(f"  P(at least one of {n_checked} shows both) = "
      f"{1-(1-q11*q13)**n_checked:.4f}  ({(1-(1-q11*q13)**n_checked)*100:.1f}%)")
print("  -> a ~1-in-7 outcome. Expected, not remarkable.")

# ---- Test 4: does delta = 2*phi/ln2 look like an identity-with-correction? ----
print("\n" + "="*72)
print("  TEST 4 -- Precision profile of the headline relations")
print("="*72)
DELTA=4.66920160910299067; PHI=(1+math.sqrt(5))/2; LN2=math.log(2)
E=math.e; GAMMA=0.57721566490153286; Z2=math.pi**2/6; S2=math.sqrt(2)
rel = [
 ("delta = 2*phi/ln2",            DELTA,        2*PHI/LN2),
 ("e*delta - ln2 = 12",           E*DELTA-LN2,  12.0),
 ("zeta(2) - sqrt2 = 3/13",       Z2-S2,        3/13),
 ("M*Tinf*e = 2/3",               0.26149721284764278*0.93754825431584375*E, 2/3),
 ("pi*Omega = e^gamma",           math.pi*0.56714329040978387, math.exp(GAMMA)),
]
print(f"  {'relation':<28} {'agrees to':>12}   {'first disagreeing digit':>24}")
for nm, a, b in rel:
    d = abs(a-b)/abs(b)
    digits = -math.log10(d)
    print(f"  {nm:<28} {digits:>9.1f} digits   position {int(digits)+1} of the decimal")
print("\n  Reference -- a relation that IS explained by a theorem:")
try:
    from mpmath import mp, mpf, exp, sqrt, pi as mpi, floor
    mp.dps = 40
    v = exp(mpi*sqrt(mpf(163))); nb = floor(v+mpf('0.5'))
    print(f"    e^(pi*sqrt(163)) vs nearest integer: agrees to "
          f"{float(-mp.log10(abs(v-nb)/v)):.1f} digits  (CM theory explains this)")
except Exception:
    print("    e^(pi*sqrt(163)) agrees to ~12 digits (CM theory explains this)")
print("\n  The session's relations agree to 3-4 digits and then stop.")
print("  A real identity-plus-correction keeps going once the correction is found.")
print("  An accident stops exactly where the search tolerance was set.")
