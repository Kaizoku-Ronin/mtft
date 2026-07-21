from sympy import symbols, groebner, Rational, expand, factor, QQ
import itertools, time
x, y, z = symbols('x y z')
F1 = (1+x*y)**3*z + y**2*(1+x*y)*(4+3*x*y)
F2 = y + 3*x*(1+x*y)**2*z + 3*x*y**2*(4+3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
u = 1+x*y

# 0) u-form identities and equivariance (exact asserts)
assert expand(F1 - (u**3*z + y**2*u + 3*y**2*u**2)) == 0
assert expand(F2 - (y + 3*x*y**2 + 9*x*y**2*u + 3*x*u**2*z)) == 0
assert expand(F3 - x*(2 - 3*x*y - x**2*z)) == 0
s = {x:-x, y:-y}
assert expand(F1.subs(s, simultaneous=True) - F1) == 0
assert expand(F2.subs(s, simultaneous=True) + F2) == 0
assert expand(F3.subs(s, simultaneous=True) + F3) == 0
print("u-form + Z/2 equivariance: verified exactly")

def fiber_count(target, note):
    t0=time.time()
    a,b,c = target
    G = groebner([F1-a, F2-b, F3-c], x, y, z, order='grevlex')
    if G.exprs == [1] or G.exprs == [Rational(1)]:
        print(f"  fiber over {note} {target}: EMPTY  ({time.time()-t0:.0f}s)"); return 0
    lms = [tuple(m.monoms()[0]) for m in [g.as_poly(x,y,z) for g in G.exprs]]
    lms = [g.as_poly(x,y,z).LM(order='grevlex').exponents for g in G.exprs]
    N = 30
    cnt = 0
    for e in itertools.product(range(N), repeat=3):
        if not any(all(e[i]>=l[i] for i in range(3)) for l in lms):
            cnt += 1
    print(f"  fiber over {note} {target}: {cnt} points (with multiplicity)  ({time.time()-t0:.0f}s)")
    return cnt

print("\n1) generic degree (fiber counts at random rational targets):")
fiber_count((Rational(3), Rational(1), Rational(2)), "generic")
fiber_count((Rational(-2), Rational(5,3), Rational(1)), "generic")

print("\n2) special targets:")
fiber_count((Rational(-1,4), Rational(0), Rational(0)), "collision")
fiber_count((Rational(7), Rational(2), Rational(0)), "wall Z=0")

print("\n3) candidate missed loci (aligned target roles (x,y,z)~(Z,Y,X)):")
# u-analog 1+Z*Y=0 ; Q-analog 2-3ZY-Z^2 X=0
fiber_count((Rational(0), Rational(1), Rational(-1)), "1+ZY=0 pt A")
fiber_count((Rational(5), Rational(1), Rational(-1)), "1+ZY=0 pt B")
fiber_count((Rational(2), Rational(0), Rational(1)), "Q_t=0 pt C")
