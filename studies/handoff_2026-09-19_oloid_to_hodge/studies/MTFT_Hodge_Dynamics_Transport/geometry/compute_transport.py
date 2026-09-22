#!/usr/bin/env python3
"""Exact Gauss--Manin reduction over Q(c), using only the standard library.

For D_2: v^2=Q(u,c)=(u-c)((u^2+c)^2+c), calculate the connection
of eta_j=u^j du/v, j=0,...,3.  Polynomial identities certify every
entry, independently of numerical integration or a choice of cycles.
"""
from fractions import Fraction as F
from pathlib import Path
import json


def trim(a):
    a = list(map(F, a))
    while len(a) > 1 and not a[-1]:
        a.pop()
    return tuple(a or [F(0)])


def padd(a, b):
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def pneg(a):
    return tuple(-x for x in a)


def pmul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def pdiv(a, b):
    a, b = list(trim(a)), trim(b)
    if b == (F(0),):
        raise ZeroDivisionError
    q = [F(0)] * max(1, len(a)-len(b)+1)
    while trim(a) != (F(0),) and len(a) >= len(b):
        k = len(a)-len(b)
        t = a[-1]/b[-1]
        q[k] += t
        for j in range(len(b)):
            a[j+k] -= t*b[j]
        a = list(trim(a))
    return trim(q), trim(a)


def pgcd(a, b):
    while b != (F(0),):
        a, b = b, pdiv(a,b)[1]
    return tuple(x/a[-1] for x in a) if a != (F(0),) else (F(1),)


class Rat:
    """Rational functions in c, stored in reduced form over Q."""
    def __init__(self, n=0, d=1):
        if isinstance(n, Rat):
            self.n, self.d = n.n, n.d
            return
        n = trim(n if isinstance(n, (tuple,list)) else [n])
        d = trim(d if isinstance(d, (tuple,list)) else [d])
        if d == (F(0),):
            raise ZeroDivisionError
        g = pgcd(n,d)
        n, rn = pdiv(n,g)
        d, rd = pdiv(d,g)
        assert not any(rn) and not any(rd)
        scale = d[-1]
        self.n = tuple(x/scale for x in n)
        self.d = tuple(x/scale for x in d)
    def __add__(self, b):
        b = Rat(b)
        return Rat(padd(pmul(self.n,b.d),pmul(b.n,self.d)), pmul(self.d,b.d))
    __radd__ = __add__
    def __neg__(self):
        return Rat(pneg(self.n),self.d)
    def __sub__(self,b):
        return self + -Rat(b)
    def __rsub__(self,b):
        return Rat(b) + -self
    def __mul__(self,b):
        b = Rat(b)
        return Rat(pmul(self.n,b.n),pmul(self.d,b.d))
    __rmul__ = __mul__
    def __truediv__(self,b):
        b = Rat(b)
        return Rat(pmul(self.n,b.d),pmul(self.d,b.n))
    def __rtruediv__(self,b):
        return Rat(b)/self
    def __pow__(self,n):
        if n < 0:
            return (Rat(1)/self)**(-n)
        ans = Rat(1)
        for _ in range(n):
            ans *= self
        return ans
    def __bool__(self):
        return any(self.n)
    def __eq__(self,b):
        b = Rat(b)
        return self.n == b.n and self.d == b.d
    def at(self,c):
        def ev(p):
            out = F(0)
            for a in reversed(p):
                out = out*c+a
            return out
        return ev(self.n)/ev(self.d)
    def encoded(self):
        return {"numerator_ascending": [str(x) for x in self.n],
                "denominator_ascending": [str(x) for x in self.d]}
    def derivative(self):
        dn = trim([i*x for i,x in enumerate(self.n)][1:])
        dd = trim([i*x for i,x in enumerate(self.d)][1:])
        return Rat(padd(pmul(dn,self.d),pneg(pmul(self.n,dd))),pmul(self.d,self.d))


def solve(a, b):
    a = [[Rat(x) for x in row]+[Rat(b[i])] for i,row in enumerate(a)]
    n = len(a)
    for k in range(n):
        candidates = [i for i in range(k,n) if a[i][k]]
        if not candidates:
            raise ValueError("singular reduction matrix")
        r = min(candidates,key=lambda i:len(a[i][k].n)+len(a[i][k].d))
        a[k],a[r] = a[r],a[k]
        pivot = a[k][k]
        a[k] = [x/pivot for x in a[k]]
        for i in range(k+1,n):
            t = a[i][k]
            if t:
                a[i] = [x-t*y for x,y in zip(a[i],a[k])]
    x = [Rat(0)]*n
    for i in reversed(range(n)):
        x[i] = a[i][-1]-sum(a[i][j]*x[j] for j in range(i+1,n))
    return x


def poly_u_mul(a,b):
    out = [Rat(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] += x*y
    return out


def poly_u_add(a,b):
    out = [(a[i] if i<len(a) else Rat(0))+(b[i] if i<len(b) else Rat(0))
           for i in range(max(len(a),len(b)))]
    while len(out)>1 and not out[-1]:
        out.pop()
    return out


def poly_u_derivative(a):
    return [i*x for i,x in enumerate(a)][1:] or [Rat(0)]


def poly_u_scale(a,s):
    return [s*x for x in a]


def poly_u_compose(a,b):
    out = [Rat(0)]
    for x in reversed(a):
        out = poly_u_add(poly_u_mul(out,b),[x])
    return out


def exact_primitive_numerator(r,q):
    return poly_u_add(poly_u_mul(poly_u_derivative(r),q),
                      poly_u_scale(poly_u_mul(r,poly_u_derivative(q)),Rat(-1)/2))


def zero_polynomial(a):
    return not any(a)


def main():
    c = Rat([0,1])
    # Q = u^5 - c u^4 + 2c u^3 - 2c^2 u^2 + (c^2+c)u - c^3-c^2.
    q = [-c**3-c**2,c**2+c,-2*c**2,2*c,-c,Rat(1)]
    qc = [-3*c**2-2*c,2*c+1,-4*c,Rat(2),Rat(-1),Rat(0)]
    qu = [(k+1)*q[k+1] for k in range(5)]
    # N = A Q + (R_u Q - R Q_u/2), with deg(A)<=3, deg(R)<=4.
    cols = []
    for k in range(4):
        cols.append([Rat(0)]*k+q+[Rat(0)]*(3-k))
    for k in range(5):
        col = [Rat(0)]*9
        if k:
            for i,z in enumerate(q):
                col[i+k-1] += k*z
        for i,z in enumerate(qu):
            col[i+k] -= z/2
        cols.append(col)
    matrix = [[cols[j][i] for j in range(9)] for i in range(9)]
    connection, certificates = [],[]
    for j in range(4):
        rhs = [Rat(0)]*9
        for i,z in enumerate(qc):
            rhs[i+j] = -z/2
        sol = solve(matrix,rhs)
        assert all(sum(matrix[i][k]*sol[k] for k in range(9)) == rhs[i]
                   for i in range(9)), "polynomial certificate failed"
        connection.append(sol[:4])
        certificates.append(sol[4:])
        print(f"Certified connection row {j}",flush=True)
    assert sum(connection[j][j] for j in range(4)) == 0
    c2 = [[-1/(4*c),1/(4*c*(c+1))],[-Rat(1)/4,1/(4*c)]]
    # Independently certify C2, using its residue-free basis {du/w,u^2 du/w}.
    p2 = [c*c+c,Rat(0),2*c,Rat(0),Rat(1)]
    p2c = [x.derivative() for x in p2]
    p2basis = [[Rat(1)],[Rat(0),Rat(0),Rat(1)]]
    c2cert = []
    for j,basis in enumerate(p2basis):
        rhs = poly_u_scale(poly_u_mul(basis,p2c),Rat(-1)/2)
        target = poly_u_add(poly_u_scale(p2basis[0],c2[j][0]),
                            poly_u_scale(p2basis[1],c2[j][1]))
        residual = poly_u_add(rhs,poly_u_scale(poly_u_mul(target,p2),Rat(-1)))
        columns = [exact_primitive_numerator([Rat(0)]*k+[Rat(1)],p2) for k in [1,3]]
        # Highest-degree and constant coefficients suffice to solve; verify all.
        def coeff(a,i): return a[i] if i<len(a) else Rat(0)
        small = [[coeff(col,i) for col in columns] for i in [0,6]]
        rr = solve(small,[coeff(residual,i) for i in [0,6]])
        primitive = [Rat(0),rr[0],Rat(0),rr[1]]
        assert zero_polynomial(poly_u_add(residual,poly_u_scale(exact_primitive_numerator(primitive,p2),Rat(-1))))
        c2cert.append(primitive)
    # Direct genus-3 verification includes parameter derivatives of the maps.
    u = [c,Rat(0),Rat(1)]
    p3 = poly_u_compose(p2,u)
    p3c = [x.derivative() for x in p3]
    c3basis = [poly_u_mul([Rat(0),Rat(2)],poly_u_compose(b,u)) for b in p2basis]
    c3basis += [poly_u_scale(poly_u_compose([Rat(0)]*j+[Rat(1)],u),Rat(2)) for j in range(4)]
    c3connection = [[Rat(0)]*6 for _ in range(6)]
    for i in range(2):
        c3connection[i][:2] = c2[i]
    for i in range(4):
        c3connection[i+2][2:] = connection[i]
    c3cert = []
    for j,b in enumerate(c3basis):
        rhs = poly_u_add(poly_u_mul([x.derivative() for x in b],p3),
                         poly_u_scale(poly_u_mul(b,p3c),Rat(-1)/2))
        target = [Rat(0)]
        for k in range(6):
            target = poly_u_add(target,poly_u_scale(c3basis[k],c3connection[j][k]))
        residual = poly_u_add(rhs,poly_u_scale(poly_u_mul(target,p3),Rat(-1)))
        if j<2:
            # A varying pullback adds d_u of the contraction with du/dc=1.
            primitive = poly_u_add(poly_u_compose(c2cert[j],u),poly_u_compose(p2basis[j],u))
        else:
            k = j-2
            numerator = poly_u_add(poly_u_compose(certificates[k],u),
                                    poly_u_compose([Rat(0)]*k+[Rat(1)],u))
            assert not numerator[0]
            primitive = numerator[1:]  # Divide by x because v=x*y.
        assert zero_polynomial(poly_u_add(residual,poly_u_scale(exact_primitive_numerator(primitive,p3),Rat(-1)))), f"C3 row {j} failed"
        c3cert.append(primitive)
    print("Certified inherited C2 and all six genus-3 pullback identities",flush=True)
    specializations = {}
    for point in [-2,-3,F(-3,2),1,2]:
        # Numerical values here are exact rationals, not floating-point results.
        specializations[str(point)] = {
            "D2_connection": [[str(x.at(point)) for x in row] for row in connection],
            "C2_connection": [[str(x.at(point)) for x in row] for row in c2],
            "certificate_R_coefficients": [[str(x.at(point)) for x in row] for row in certificates],
        }
    payload = {
        "family": "v^2=(u-c)*((u^2+c)^2+c)",
        "basis": [f"u^{j} du/v" for j in range(4)],
        "convention": "d_c eta_j = sum_k A[j][k] eta_k + d_u(R_j/v)",
        "connection": [[x.encoded() for x in row] for row in connection],
        "certificate_R_coefficients": [[x.encoded() for x in row] for row in certificates],
        "specializations": specializations,
        "trace_zero": True,
        "polynomial_identity_certificates_passed": 4,
        "coefficient_equations_per_identity": 9,
        "additional_C2_identities_passed": 2,
        "direct_C3_pullback_identities_passed": 6,
        "C3_exact_primitive_coefficients": [[x.encoded() for x in row] for row in c3cert],
        "genus3_split_basis": ["2x dx/y", "2x(x^2+c)^2 dx/y"] +
                              [f"2(x^2+c)^{j} dx/y" for j in range(4)],
        "genus3_connection": "block_diag(C2_connection,D2_connection)",
        "projector": "sigma^*(x^j dx/y)=(-1)^(j+1) x^j dx/y; P_+=(1+sigma^*)/2, P_-=(1-sigma^*)/2",
        "valid_parameters": "c*(c+1)*(c^3+2c^2+c+1) != 0",
    }
    dest = Path(__file__).resolve().parent / "transport_results.json"
    dest.write_text(json.dumps(payload,indent=2)+"\n")
    print(f"Wrote {dest}",flush=True)


if __name__ == "__main__":
    main()
