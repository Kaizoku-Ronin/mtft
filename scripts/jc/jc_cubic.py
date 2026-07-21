from sympy import symbols, resultant, factor, expand, Poly, Rational, discriminant, groebner
import time
t0=time.time()
x, y = symbols('x y')
A, B, C = symbols('A B C')
u = 1+x*y
# z eliminated: z = (A - y^2 u(1+3u))/u^3  and  z = (2x-3x^2 y - C)/x^3
E1 = expand(x**3*(A - y**2*u*(1+3*u)) - u**3*(2*x - 3*x**2*y - C))
E2 = expand(x**2*(y - B) + 3*x**3*y**2*(1+3*u) + 3*u**2*(2*x - 3*x**2*y - C))
R = resultant(Poly(E1, y), Poly(E2, y))
print(f"resultant computed ({time.time()-t0:.0f}s)")
Rf = factor(R)
print("factored resultant:")
print(Rf)
