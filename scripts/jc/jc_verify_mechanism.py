from sympy import symbols, expand, factor, Rational, groebner, sqrt
import itertools, time
x,y,z,T,A,B,C = symbols('x y z T A B C')
F1 = (1+x*y)**3*z + y**2*(1+x*y)*(4+3*x*y)
F2 = y + 3*x*(1+x*y)**2*z + 3*x*y**2*(4+3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
p3 = 27*A**2*C**2 - 18*A*B*C + B**3*C + 16*A - B**2
p1 = 4 - 3*B*C
p0 = -2*C

# (a) tautological identity: P(x; F1,F2,F3) == 0
taut = expand((p3*T**3 + p1*T + p0).subs({T:x, A:F1, B:F2, C:F3}))
print("tautological identity P(x;F)=0:", taut == 0)

# (b) discriminant of the cubic: square or not -> monodromy
Delta = expand(-4*p3*p1**3 - 27*p3**2*p0**2)
print("Delta = -p3*(4 p1^3 + 27 p3 p0^2); factored:")
print("  ", factor(Delta))

# (c) pullback of the escape wall: factor p3(F1,F2,F3)
pb = factor(p3.subs({A:F1,B:F2,C:F3}, simultaneous=True))
print("p3 o F factored:")
print("  ", pb)

def fiber_count(target, note):
    a,b,c = target
    G = groebner([F1-a, F2-b, F3-c], x, y, z, order='grevlex')
    if list(G.exprs) == [1]:
        print(f"  fiber over {note} {target}: EMPTY"); return
    lms = [g.as_poly(x,y,z).LM(order='grevlex').exponents for g in G.exprs]
    cnt = sum(1 for e in itertools.product(range(30),repeat=3)
              if not any(all(e[i]>=l[i] for i in range(3)) for l in lms))
    print(f"  fiber over {note} {target}: {cnt} points")

print("\npredictions on the escape wall:")
fiber_count((0,0,2), "p3=0 axis pt (predict 1)")
fiber_count((Rational(1,3), 2, Rational(2,3)), "p3=p1=0 curve pt (predict EMPTY)")
# generic p3=0, p1!=0 point: pick B=1, solve p3=0 for A with C=1: 27A^2-18A+1+16A-1=27A^2-2A=A(27A-2)
fiber_count((Rational(2,27), 1, 1), "p3=0 generic (predict 1)")
fiber_count((0,1,1), "p3=0 (A=0,B=C=1: p3=1-1=0) (predict 1)")
