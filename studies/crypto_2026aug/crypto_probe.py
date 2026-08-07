"""
crypto_probe — what the v0.11.x structures actually offer cryptography.

Q1  MODULE TYPE (answered in run 1): every isotypic block of T2 on the
    29-dim Manin quotient appears with multiplicity exactly 2.  So
    H_1^new is a rank-2 module over the Hecke algebra: the Kyber shape.

Q2  GLUING.  T (x) Q = K1 x K2 x K3 splits.  Does the INTEGRAL Hecke
    order T split too?  If it does, an attacker projects onto the
    smallest factor and wins.  If it is GLUED at congruence primes,
    the projection is not integral and the module does not decompose
    over Z.  Measured: disc(Z[T2]|new) vs prod_i disc(Z[a2]|K_i).
    The ratio's primes are exactly the congruence primes.

Q3  EXPANSION.  Is the Hecke correspondence graph on P^1(Z/N)
    Ramanujan?  That is the CGL/LPS hash requirement.  Measured
    directly against the bound 2 sqrt(p).
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import sympy as sp
from x0143_particle_box import P1, tessellation, ModularSymbols, hecke_matrices

N = 143
x = sp.Symbol('x')

# the three new-part factors of charpoly(T2), read off in run 1
F1 = sp.Poly(x, x)                                          # a2 = 0
F2 = sp.Poly(x**4 - 3*x**3 - x**2 + 5*x + 1, x)
F3 = sp.Poly(x**6 - 10*x**4 + 2*x**3 + 24*x**2 - 7*x - 12, x)


def main():
    p1 = P1(N)
    tessellation(N)
    ms = ModularSymbols(p1)
    print(f"P^1(Z/{N}) = {len(p1.reps)} cosets   Manin dim {ms.dim}   "
          f"cuspidal {ms.cuspidal_dim}")

    # ---- Q2: congruence gluing of the integral Hecke order --------------
    print("\n" + "=" * 70)
    print("Q2  GLUING — does the integral Hecke order split over Z?")
    print("=" * 70)

    d1 = 1                                    # disc of Q
    d2 = int(sp.discriminant(F2.as_expr(), x))
    d3 = int(sp.discriminant(F3.as_expr(), x))
    print(f"  disc Z[a2] for K1 (deg 1) = {d1}")
    print(f"  disc Z[a2] for K2 (deg 4) = {d2}   = {sp.factorint(d2)}")
    print(f"  disc Z[a2] for K3 (deg 6) = {d3}   = {sp.factorint(d3)}")

    # New-part minimal polynomial of T2 and its discriminant.
    Pnew = (F1 * F2 * F3).as_expr()
    dnew = int(sp.discriminant(Pnew, x))
    print(f"\n  disc Z[T2] on the whole new part = {dnew}")
    print(f"    = {sp.factorint(abs(dnew))}")

    ratio = sp.Rational(abs(dnew), abs(d1 * d2 * d3))
    print(f"\n  disc(new) / prod_i disc(K_i) = {ratio}")
    print(f"    = {sp.factorint(int(ratio))}")
    print("\n  The primes above are the CONGRUENCE PRIMES: the resultants")
    print("  Res(F_i, F_j) measure mod-l congruences between the orbits.")
    for (a, na), (b, nb) in ((("F1", F1), ("F2", F2)),
                             (("F1", F1), ("F3", F3)),
                             (("F2", F2), ("F3", F3))):
        r = int(sp.resultant(na.as_expr(), nb.as_expr(), x))
        print(f"    Res({a},{b}) = {r:>12}   = {sp.factorint(abs(r))}")

    # ---- Q3: Ramanujan / expander test on the coset graph ---------------
    print("\n" + "=" * 70)
    print("Q3  EXPANSION — Hecke correspondence graph on P^1(Z/N)")
    print("=" * 70)
    print("   p | deg p+1 | 2sqrt(p) | 2nd |lam| (cuspidal) | Ramanujan?")
    print("  ---+---------+----------+----------------------+-----------")
    n = len(p1.reps)
    for p in (2, 3, 5, 7, 17, 19, 23, 29, 31, 37, 41, 43):
        if N % p == 0:
            continue
        # (p+1)-regular graph: x -> x.g for each Hecke matrix g
        A = np.zeros((n, n))
        for g in hecke_matrices(p):
            for i, xr in enumerate(p1.reps):
                j = p1.index[p1.act(xr, tuple(map(tuple, np.array(g))))]
                A[j, i] += 1.0
        ev = np.linalg.eigvals(A)
        mags = np.sort(np.abs(ev))[::-1]
        triv = p + 1.0
        nontriv = [m for m in mags if abs(m - triv) > 1e-7 * triv]
        lam2 = max(nontriv) if nontriv else 0.0
        bound = 2.0 * np.sqrt(p)
        print(f"  {p:>3}| {triv:>7.1f} | {bound:>8.4f} | {lam2:>20.6f} |"
              f" {'YES' if lam2 <= bound + 1e-7 else 'NO'}")

    print("\n  (deg = p+1 is the trivial/Eisenstein eigenvalue; the test is")
    print("   whether every OTHER eigenvalue sits inside the Ramanujan disc.)")


if __name__ == "__main__":
    main()
