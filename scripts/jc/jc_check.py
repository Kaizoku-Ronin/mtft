from sympy import symbols, Matrix, expand, Rational, groebner, solve, simplify
x, y, z = symbols('x y z')

F1 = (1+x*y)**3 * z + y**2 * (1+x*y) * (4+3*x*y)
F2 = y + 3*x*(1+x*y)**2 * z + 3*x*y**2 * (4+3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z

J = Matrix([[F.diff(v) for v in (x,y,z)] for F in (F1,F2,F3)])
det = expand(J.det())
print("Jacobian determinant (fully expanded):", det)
print("constant == -2 ?", det == -2)

pts = [(0, 0, Rational(-1,4)), (1, Rational(-3,2), Rational(13,2)), (-1, Rational(3,2), Rational(13,2))]
print("\nimages (exact rational arithmetic):")
for p in pts:
    img = tuple(expand(F.subs({x:p[0], y:p[1], z:p[2]})) for F in (F1,F2,F3))
    print(f"  F{p} = {img}")

print("\npoints pairwise distinct:", len(set(pts)) == 3)

# full fiber over (-1/4, 0, 0) via Groebner basis
print("\nfull fiber over (-1/4,0,0):")
gb = groebner([F1 + Rational(1,4), F2, F3], x, y, z, order='lex')
sols = solve(list(gb.exprs), [x, y, z], dict=True)
for s in sols:
    print("  ", s)
print("fiber size (over C, as found):", len(sols))
