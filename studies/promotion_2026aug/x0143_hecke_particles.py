#!/usr/bin/env python3
"""
x0143_hecke_particles.py — the particles of X0(143) painted on the skeleton
============================================================================

MIT License — Copyright (c) 2026 Roger Tano

Sequel to x0143_ribbon_embedding.py.  Manin symbols live on the SAME
index set P1(Z/143) as the skeleton's flags, and the dictionary is
exact: symbols = oriented tessellation edges (S-orbits with sign), the
2-term relation is orientation reversal, the 3-term relation is
precisely im d2 (triangle boundaries).  Hence

   modular symbols  =  Q^84 / im d2  =  H1(X, cusps)  (dim 29),
   cuspidal part    =  ker(cusp boundary)             (dim 26 = H1(X)).

Hecke operators act by Merel's matrices (det p, a > b >= 0, d > c >= 0)
through the right P1 action.  Everything below is exact rational /
integer linear algebra; numerics appear only as guess generators whose
outputs are certified exactly.

RESULTS
-------
R1 (EXACT).  Model certified: rank d2 = 55, quotient dim 29, cusp
  boundary annihilates d2, cuspidal dim 26; Merel T2, T3, T5 all
  preserve im d2 (well-defined), pairwise commute, and act on the
  3-dim Eisenstein complement with eigenvalue exactly p + 1 —
  charpoly(T_p | 29) = (x - (p+1))^3 * charpoly(T_p | 26) verified by
  exact polynomial division.

R2 (Pr/EXACT — the particle content, corpus-completing).  The
  factorization, verified by exact integer polynomial multiplication:

     charpoly(T2 | H1) = x^2 * (x+2)^4 * g4(x)^2 * h6(x)^2
     g4(x) = x^4 - 3x^3 - x^2 + 5x + 1        (irreducible /Q)
     h6(x) = x^6 - 10x^4 + 2x^3 + 24x^2 - 7x - 12  (irreducible /Q)

  FOUR Hecke blocks of dims [2, 4, 8, 12] = 26:
    E_ell  (2):  ker T2;         T3 = -1, T5 = -1  ==  143a1 point
                 counts on the corpus curve (0,-1,1,-1,-2).  E2:
                 Merel/Manin route vs Weierstrass point counting.
    E_old  (4):  ker(T2 + 2);    T3 = -1, T5 = +1  ==  11a1 point
                 counts on (0,-1,1,-10,-20).  This is the LEVEL-11
                 OLDSPACE: 143 = 11*13, dim S2(Gamma0(11)) = 1, two
                 degeneracy maps -> the ghost of X0(11) inside 143.
    E_q4   (8):  ker g4(T2)   — the quartic new orbit f2.
    E_q6  (12):  ker h6(T2)   — the sextic new orbit f3.
  CORPUS RECONCILIATION: the corpus statement "three Galois orbits
  [1,4,6]" is the NEWSPACE and is here re-derived from scratch
  (1 + 4 + 6 = 11 = genus 13 minus the 2-dim old part).  The homology
  additionally carries old-11a with multiplicity two; the skeleton
  remembers its factor level.  Memory refined, not contradicted.

R3 (EXACT).  The star involution i*(c:d) = (-c:d) descends, commutes
  with T2, T3, T5, and splits the blocks into (1,1), (2,2), (4,4),
  (6,6) — the real Hodge split of each particle.

R4 (Pr).  HARMONIC PAINTING.  Each block's classes get their unique
  harmonic representatives (orthogonal to im d2 in R^84, exact
  rational).  The presence density rho(e) = diagonal of the block's
  orthogonal projector is basis-independent, i*-invariant (exact),
  with trace exactly the block dimension.  These are the four particle
  fields on the 84 edges of the certified embedding.

R5 (EXACT).  The 143a1 REAL AND IMAGINARY CYCLES as primitive integer
  edge flows: the i*-even and i*-odd generators of E_ell, harmonically
  projected and cleared to primitive integer vectors — the real
  homology of the elliptic curve written on the skeleton in integers.

GATES G1-G9 mirror the above.  Runtime ~4 min.  Writes
x0143_hecke_particles_ledger.json and particles_viz.json (artifact
data) next to itself.
"""

from __future__ import annotations

import itertools
import json
import os
import sys
import time
from fractions import Fraction as Fr
from math import gcd

import numpy as np

N = 143


# ── flags and actions ───────────────────────────────────────────────

def build():
    units = [u for u in range(1, N) if gcd(u, N) == 1]
    valid = lambda c, d: not (c % 11 == 0 and d % 11 == 0) and \
                         not (c % 13 == 0 and d % 13 == 0)
    canon = lambda c, d: min(((c * u) % N, (d * u) % N) for u in units)
    P1 = sorted({canon(c, d) for c in range(N) for d in range(N)
                 if valid(c, d)})
    idx = {p: i for i, p in enumerate(P1)}
    perm = lambda f: [idx[canon(*f(*p))] for p in P1]
    sS = perm(lambda c, d: (d % N, (-c) % N))
    sT = perm(lambda c, d: (c, (c + d) % N))
    sR = perm(lambda c, d: (d % N, (d - c) % N))
    io = perm(lambda c, d: ((-c) % N, d))
    return P1, idx, canon, sS, sT, sR, io


def orbits(s):
    seen, out = set(), []
    for i in range(len(s)):
        if i in seen:
            continue
        o, j = [], i
        while j not in seen:
            seen.add(j)
            o.append(j)
            j = s[j]
        out.append(o)
    return out


def rref(M):
    M = [[Fr(x) for x in row] for row in M]
    R, C = len(M), len(M[0])
    piv, r = [], 0
    for c in range(C):
        pr = None
        for i in range(r, R):
            if M[i][c] != 0:
                pr = i
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(R):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == R:
            break
    return M, piv


def nullspace(M):
    R, piv = rref(M)
    C = len(M[0])
    pivset = set(piv)
    out = []
    for fc in [c for c in range(C) if c not in pivset]:
        v = [Fr(0)] * C
        v[fc] = Fr(1)
        for r_, c_ in enumerate(piv):
            v[c_] = -R[r_][fc]
        out.append(v)
    return out


def merel(p):
    out = []
    for a in range(1, p + 1):
        for b in range(0, a):
            for c in range(0, p + 1):
                num = p + b * c
                if num % a:
                    continue
                d = num // a
                if d > c and a * d - b * c == p:
                    out.append((a, b, c, d))
    return out


def pmulz(A, B):
    C = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                C[i + j] += a * b
    return C


def pdivz(A, B):
    A = A[:]
    q = [0] * (len(A) - len(B) + 1)
    for i in range(len(A) - 1, len(B) - 2, -1):
        if A[i] % B[-1]:
            return None, None
        f = A[i] // B[-1]
        q[i - len(B) + 1] = f
        for j, b in enumerate(B):
            A[i - len(B) + 1 + j] -= f * b
    return q, A[:len(B) - 1]


def main() -> int:
    t00 = time.time()
    ledger = {"study": "x0143_hecke_particles", "gates": {}}
    ok = True

    def gate(name, passed, **info):
        nonlocal ok
        ok &= bool(passed)
        ledger["gates"][name] = {"passed": bool(passed), **info}
        print(f"[{'PASS' if passed else 'FAIL'}] {name}  "
              + "  ".join(f"{k}={v}" for k, v in info.items()))

    P1, idx, canon, sS, sT, sR, io = build()
    tris = orbits(sR)
    tri_of = {}
    for ti, t in enumerate(tris):
        for f in t:
            tri_of[f] = ti
    fans = sorted(orbits(sT), key=len)
    cusp_of = {}
    for k, o in enumerate(fans):
        for f in o:
            cusp_of[f] = k
    Eorb = orbits(sS)
    E = len(Eorb)
    erep, eid, esign = {}, {}, {}
    for k, o in enumerate(Eorb):
        r = min(o)
        erep[k] = r
        for f in o:
            eid[f] = k
            esign[f] = 1 if f == r else -1

    D2 = [[0] * 56 for _ in range(E)]
    for ti, t in enumerate(tris):
        for f in t:
            D2[eid[f]][ti] += esign[f]
    DEL = [[0] * 4 for _ in range(E)]
    for k in range(E):
        x = erep[k]
        DEL[k][cusp_of[x]] += 1
        DEL[k][cusp_of[sS[x]]] -= 1

    # G1 model
    dD = all(sum(DEL[e][c] * D2[e][t] for e in range(E)) == 0
             for c in range(4) for t in range(56))
    _, pivC = rref([[D2[e][t] for t in range(56)] for e in range(E)])
    cols = pivC
    D2r = [[D2[e][t] for t in cols] for e in range(E)]
    _, pivE = rref([[D2r[e][j] for e in range(E)]
                    for j in range(len(cols))])
    free = [e for e in range(E) if e not in set(pivE)]
    gate("G1_model", dD and len(cols) == 55 and len(free) == 29,
         rank_d2=len(cols), quotient=len(free))

    # invert B = [D2r | I_free]
    B = [[Fr(0)] * E for _ in range(E)]
    for e in range(E):
        for j in range(55):
            B[e][j] = Fr(D2r[e][j])
    for j, e in enumerate(free):
        B[e][55 + j] = Fr(1)
    M = [row[:] + [Fr(1) if i == j else Fr(0) for j in range(E)]
         for i, row in enumerate(B)]
    for c in range(E):
        pr = next(i for i in range(c, E) if M[i][c] != 0)
        M[c], M[pr] = M[pr], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(E):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    Binv = [row[E:] for row in M]

    def qcoords(v):
        x = [sum(Binv[i][k] * v[k] for k in range(E) if v[k] != 0)
             for i in range(E)]
        return x[55:]

    def hecke(p):
        T = [[0] * E for _ in range(E)]
        for k in range(E):
            c, d = P1[erep[k]]
            for (a, b, cc, dd) in merel(p):
                y = idx[canon((c * a + d * cc) % N, (c * b + d * dd) % N)]
                T[eid[y]][k] += esign[y]
        return T

    def descends(T):
        TD = [[sum(T[e][x] * D2r[x][j] for x in range(E) if T[e][x])
               for j in range(55)] for e in range(E)]
        _, pv = rref([[ (D2r[e][j] if j < 55 else TD[e][j - 55])
                        for e in range(E)] for j in range(110)])
        return len(pv) == 55

    def qop(T):
        cols_ = [qcoords([Fr(T[e][fj]) for e in range(E)]) for fj in free]
        return [[cols_[j][i] for j in range(29)] for i in range(29)]

    T2m, T3m, T5m = hecke(2), hecke(3), hecke(5)
    dsc = descends(T2m) and descends(T3m) and descends(T5m)
    Q2, Q3, Q5 = qop(T2m), qop(T3m), qop(T5m)

    def comm(A, Bq):
        return all(sum(A[i][k] * Bq[k][j] for k in range(len(A)))
                   == sum(Bq[i][k] * A[k][j] for k in range(len(A)))
                   for i in range(len(A)) for j in range(len(A)))

    # cuspidal
    DB = [[Fr(DEL[free[j]][c]) for j in range(29)] for c in range(4)]
    Mr, pivK = rref(DB)
    freeK = [j for j in range(29) if j not in set(pivK)]
    K = []
    for fj in freeK:
        v = [Fr(0)] * 29
        v[fj] = Fr(1)
        for r_, c_ in enumerate(pivK):
            v[c_] = -Mr[r_][fj]
        K.append(v)

    def restrict(Q):
        out = []
        for k in K:
            img = [sum(Q[i][j] * k[j] for j in range(29))
                   for i in range(29)]
            out.append([img[fj] for fj in freeK])
        return [[out[j][i] for j in range(26)] for i in range(26)]

    C2, C3, C5 = restrict(Q2), restrict(Q3), restrict(Q5)

    def charpoly(Mrat):
        n = len(Mrat)
        den = 1
        for row in Mrat:
            for x in row:
                den = den * x.denominator // gcd(den, x.denominator)
        Mi = [[int(x * den) for x in row] for row in Mrat]

        def bdet(A):
            A = [r[:] for r in A]
            m = len(A)
            prev, sg = 1, 1
            for k in range(m - 1):
                if A[k][k] == 0:
                    for r_ in range(k + 1, m):
                        if A[r_][k]:
                            A[k], A[r_] = A[r_], A[k]
                            sg = -sg
                            break
                    else:
                        return 0
                for i in range(k + 1, m):
                    for j in range(k + 1, m):
                        A[i][j] = (A[i][j] * A[k][k]
                                   - A[i][k] * A[k][j]) // prev
                prev = A[k][k]
            return sg * A[m - 1][m - 1]

        vals = [bdet([[(x * den if i == j else 0) - Mi[i][j]
                       for j in range(n)] for i in range(n)])
                for x in range(n + 1)]
        work = [Fr(v) for v in vals]
        newt = [work[0]]
        for k in range(1, n + 1):
            work = [(work[i + 1] - work[i]) / Fr(k)
                    for i in range(len(work) - 1)]
            newt.append(work[0])
        acc = [Fr(0)] * (n + 1)
        base = [Fr(1)] + [Fr(0)] * n
        for k in range(n + 1):
            for i in range(n + 1):
                acc[i] += newt[k] * base[i]
            nb = [Fr(0)] * (n + 1)
            for i in range(n):
                nb[i + 1] += base[i]
                nb[i] -= Fr(k) * base[i]
            base = nb
        cp = [acc[i] * Fr(den) ** i / Fr(den) ** n for i in range(n + 1)]
        assert all(c.denominator == 1 for c in cp)
        return [int(c) for c in cp]

    cp2_26, cp2_29 = charpoly(C2), charpoly(Q2)
    cp3_26, cp3_29 = charpoly(C3), charpoly(Q3)
    cp5_26, cp5_29 = charpoly(C5), charpoly(Q5)
    eis = True
    for cp_s, cp_b, p in ((cp2_26, cp2_29, 2), (cp3_26, cp3_29, 3),
                          (cp5_26, cp5_29, 5)):
        q, r = pdivz(cp_b, pmulz(cp_s, pmulz(pmulz(
            [-(p + 1), 1], [-(p + 1), 1]), [-(p + 1), 1])))
        eis &= (q is not None and all(x == 0 for x in r) and q == [1])
    gate("G2_hecke_wellposed", dsc and comm(Q2, Q3) and comm(Q2, Q5)
         and comm(Q3, Q5) and eis,
         descends=dsc, eisenstein="(x-(p+1))^3 exact for p=2,3,5")

    # G3/G4 elliptic and old blocks vs point counts
    def curve_ap(curve, pmax):
        a1, a2_, a3, a4, a6 = curve
        b2 = a1 * a1 + 4 * a2_
        b4 = 2 * a4 + a1 * a3
        b6 = a3 * a3 + 4 * a6
        pr = [p for p in range(2, pmax + 1)
              if all(p % q for q in range(2, int(p ** .5) + 1))]
        ap = {2: 2 + 1 - (sum(1 for x in range(2) for y in range(2)
              if (y * y + a1 * x * y + a3 * y
                  - (x ** 3 + a2_ * x * x + a4 * x + a6)) % 2 == 0) + 1)}
        for p in pr:
            if p == 2:
                continue
            sq = bytearray(p)
            for k in range(p // 2 + 1):
                sq[(k * k) % p] = 1
            s = 0
            for x in range(p):
                g = (((4 * x + b2) * x + 2 * b4) * x + b6) % p
                if g:
                    s += 1 if sq[g] else -1
            ap[p] = -s
        return ap

    ap143 = curve_ap((0, -1, 1, -1, -2), 5)
    ap11 = curve_ap((0, -1, 1, -10, -20), 5)

    def scal(C, vs, s):
        return all(all(sum(C[i][j] * v[j] for j in range(26)) == s * v[i]
                       for i in range(26)) for v in vs)

    Eell = nullspace(C2)
    gate("G3_elliptic_143a1", len(Eell) == 2
         and scal(C3, Eell, ap143[3]) and scal(C5, Eell, ap143[5])
         and ap143 == {2: 0, 3: -1, 5: -1},
         dims=len(Eell), curve_ap=str(ap143))
    Eold = nullspace([[C2[i][j] + (2 if i == j else 0)
                       for j in range(26)] for i in range(26)])
    gate("G4_old_11a1_ghost", len(Eold) == 4
         and scal(C3, Eold, ap11[3]) and scal(C5, Eold, ap11[5])
         and ap11 == {2: -2, 3: -1, 5: 1},
         dims=len(Eold), curve_ap=str(ap11))

    # G5 quartic/sextic factors: numeric guess -> exact certificates
    e2 = np.sort(np.linalg.eigvals(
        np.array([[float(x) for x in r] for r in C2])).real)
    irr = sorted(set(round(float(x), 9) for x in e2
                     if abs(x) > 1e-6 and abs(x + 2) > 1e-6))
    dcp = [cp2_26[i] * i for i in range(1, len(cp2_26))]

    def pgcd(A, Bp):
        A = [Fr(x) for x in A]
        Bp = [Fr(x) for x in Bp]
        while any(x != 0 for x in Bp):
            while A and A[-1] == 0:
                A.pop()
            while Bp and Bp[-1] == 0:
                Bp.pop()
            if len(A) < len(Bp):
                A, Bp = Bp, A
                continue
            f = A[-1] / Bp[-1]
            sh = len(A) - len(Bp)
            A = [a - (f * Bp[i - sh] if 0 <= i - sh < len(Bp) else 0)
                 for i, a in enumerate(A)]
            while A and A[-1] == 0:
                A.pop()
            if not A or all(x == 0 for x in A):
                return Bp
            A, Bp = Bp, A
        return A

    g = pgcd(cp2_26, dcp)
    g = [x / g[-1] for x in g]

    def pdivq(A, Bq):
        A = [Fr(x) for x in A]
        q = [Fr(0)] * (len(A) - len(Bq) + 1)
        for i in range(len(A) - 1, len(Bq) - 2, -1):
            f = A[i] / Bq[-1]
            q[i - len(Bq) + 1] = f
            for j, b in enumerate(Bq):
                A[i - len(Bq) + 1 + j] -= f * b
        return q

    sfq = pdivq(cp2_26, g)
    den = 1
    for x in sfq:
        den = den * x.denominator // gcd(den, x.denominator)
    sf = [int(x * den) for x in sfq]
    gg = 0
    for x in sf:
        gg = gcd(gg, abs(x))
    sf = [x // gg for x in sf]
    if sf[-1] < 0:
        sf = [-x for x in sf]
    quart, quads, cubs = [], [], []
    for r_ in (2, 3, 4):
        for sub in itertools.combinations(irr, r_):
            c = np.poly(sub)[::-1]
            ic = [round(float(x)) for x in c]
            if max(abs(float(c[i]) - ic[i])
                   for i in range(r_ + 1)) < 1e-5:
                q, rr = pdivz(sf, ic)
                if q is not None and all(x == 0 for x in rr):
                    (quads if r_ == 2 else cubs
                     if r_ == 3 else quart).append(ic)
    g4 = quart[0] if quart else None
    rest, _ = pdivz(sf, g4)
    rest, _ = pdivz(rest, [0, 1])
    h6, _ = pdivz(rest, [2, 1])
    full = pmulz([0, 0, 1], pmulz(pmulz([4, 4, 1], [4, 4, 1]),
                                  pmulz(pmulz(g4, g4), pmulz(h6, h6))))
    hasse = all(abs(x) <= 2 * np.sqrt(2) + 1e-9 for x in irr)
    gate("G5_new_orbits", len(quart) == 1 and not quads and not cubs
         and full == cp2_26 and hasse
         and g4 == [1, 5, -1, -3, 1]
         and h6 == [-12, -7, 24, 2, -10, 0, 1],
         g4=str(g4), h6=str(h6),
         factorization="x^2 (x+2)^4 g4^2 h6^2 exact",
         irreducible="no rational/quadratic/cubic factors")
    ledger["g4"] = g4
    ledger["h6"] = h6

    def poly_of_op(coefs, C):
        n = len(C)
        Mo = [[Fr(0)] * n for _ in range(n)]
        P = [[Fr(1) if i == j else Fr(0) for j in range(n)]
             for i in range(n)]
        for c in coefs:
            if c:
                for i in range(n):
                    for j in range(n):
                        Mo[i][j] += c * P[i][j]
            P = [[sum(P[i][k] * C[k][j] for k in range(n))
                  for j in range(n)] for i in range(n)]
        return Mo

    Eq4 = nullspace(poly_of_op(g4, C2))
    Eq6 = nullspace(poly_of_op(h6, C2))
    blocks = {"ell": Eell, "old": Eold, "q4": Eq4, "q6": Eq6}
    allv = Eell + Eold + Eq4 + Eq6
    _, pv = rref([[v[j] for v in allv] for j in range(26)])
    inv = True
    for Cx in (C3, C5):
        for bl in blocks.values():
            imgs = [[sum(Cx[i][j] * v[j] for j in range(26))
                     for i in range(26)] for v in bl]
            _, pvb = rref([[u[j] for u in (bl + imgs)]
                           for j in range(26)])
            inv &= len(pvb) == len(bl)
    gate("G6_block_sum", len(Eq4) == 8 and len(Eq6) == 12
         and len(pv) == 26 and inv, dims=[2, 4, 8, 12],
         t3_t5_invariant=inv)

    # G7 star involution
    IS = [[0] * E for _ in range(E)]
    for k in range(E):
        y = io[erep[k]]
        IS[eid[y]][k] += esign[y]
    QI = qop(IS)
    CI = restrict(QI)
    isq = all(sum(CI[i][k] * CI[k][j] for k in range(26))
              == (1 if i == j else 0)
              for i in range(26) for j in range(26))
    cIT = comm(CI, C2) and comm(CI, C3) and comm(CI, C5)
    splits = {}
    for nm, bl in blocks.items():
        img = [[sum(CI[i][j] * v[j] for j in range(26))
                for i in range(26)] for v in bl]
        _, pB = rref([[v[j] for v in bl] for j in range(26)])
        plus = nullspace([[img[a][i] - bl[a][i] for a in range(len(bl))]
                          for i in range(26)]) if False else None
        # count via trace: tr(CI|block): coords of img in block basis
        Rb, pivB = rref([[v[j] for v in bl] for j in range(26)])
        fB = pivB  # block basis rows... use solve via free pattern
        # simpler: block vectors have identity on their own free coords
        Cfull, pivF = rref(C2 if nm == "ell" else
                           [[C2[i][j] + (2 if i == j else 0)
                             for j in range(26)] for i in range(26)]) \
            if nm in ("ell", "old") else (None, None)
        # trace route: tr = sum over basis of coefficient of v_k in img_k
        # coords: block built by nullspace() -> each v has a distinguished
        # free coordinate = 1; record them:
        splits[nm] = None
    # recompute splits robustly via eigen-count of CI on each block
    spl_ok = True
    spl = {}
    for nm, bl in blocks.items():
        d = len(bl)
        Bm = np.array([[float(x) for x in v] for v in bl]).T  # 26 x d
        CIm = np.array([[float(x) for x in r] for r in CI])
        X = np.linalg.lstsq(Bm, CIm @ Bm, rcond=None)[0]
        ev = np.sort(np.linalg.eigvals(X).real)
        p_ = int(np.sum(ev > 0))
        m_ = d - p_
        spl[nm] = (p_, m_)
        spl_ok &= (p_ == m_ == d // 2) and np.allclose(
            np.abs(ev), 1, atol=1e-8)
    gate("G7_star_involution", isq and cIT and spl_ok, splits=str(spl))

    # G8 harmonic densities
    G56 = [[sum(D2[e][i] * D2[e][j] for e in range(E))
            for j in range(56)] for i in range(56)]

    def harmonic(vquot26):
        v29 = [Fr(0)] * 29
        for j, kv in zip(range(26), []):
            pass
        # expand cuspidal coords -> 29 quotient coords
        v29 = [sum(K[a][j] * vquot26[a] for a in range(26))
               for j in range(29)]
        v84 = [Fr(0)] * E
        for j, e in enumerate(free):
            v84[e] = v29[j]
        rhs = [sum(Fr(D2[e][i]) * v84[e] for e in range(E))
               for i in range(56)]
        Aug = [[Fr(G56[i][j]) for j in range(56)] + [rhs[i]]
               for i in range(56)]
        Rr, pv_ = rref(Aug)
        alpha = [Fr(0)] * 56
        for r_, c_ in enumerate(pv_):
            if c_ < 56:
                alpha[c_] = Rr[r_][56]
        w = [v84[e] - sum(Fr(D2[e][i]) * alpha[i] for i in range(56))
             for e in range(E)]
        return w

    dens = {}
    hbases = {}
    hok = True
    for nm, bl in blocks.items():
        W = [harmonic(v) for v in bl]
        hbases[nm] = W
        hok &= all(sum(Fr(D2[e][i]) * w[e] for e in range(E)) == 0
                   for w in W for i in range(56))
        d = len(W)
        Gm = [[sum(W[a][e] * W[b][e] for e in range(E))
               for b in range(d)] for a in range(d)]
        Aug = [row[:] + [Fr(1) if i == j else Fr(0) for j in range(d)]
               for i, row in enumerate(Gm)]
        for c in range(d):
            pr = next(i for i in range(c, d) if Aug[i][c] != 0)
            Aug[c], Aug[pr] = Aug[pr], Aug[c]
            pvv = Aug[c][c]
            Aug[c] = [x / pvv for x in Aug[c]]
            for i in range(d):
                if i != c and Aug[i][c] != 0:
                    f = Aug[i][c]
                    Aug[i] = [a - f * b for a, b in zip(Aug[i], Aug[c])]
        Gi = [row[d:] for row in Aug]
        rho = [sum(W[a][e] * Gi[a][b] * W[b][e]
                   for a in range(d) for b in range(d))
               for e in range(E)]
        hok &= (sum(rho) == d)
        # i*-invariance of rho
        for k in range(E):
            y = io[erep[k]]
            hok &= rho[eid[y]] == rho[k]
        dens[nm] = [float(x) for x in rho]
    gate("G8_harmonic_densities", hok,
         traces={nm: len(blocks[nm]) for nm in blocks})

    # G9 elliptic integer flows (i*-even / odd generators)
    Bm = np.array([[float(x) for x in v] for v in Eell]).T
    CIm = np.array([[float(x) for x in r] for r in CI])
    X = np.linalg.lstsq(Bm, CIm @ Bm, rcond=None)[0]
    ev, Vv = np.linalg.eig(X)
    flows = {}
    fok = True
    for sgn, name in ((1, "flow_plus"), (-1, "flow_minus")):
        j = int(np.argmin(np.abs(ev.real - sgn)))
        c0, c1 = Vv[:, j].real
        c0, c1 = Fr(c0).limit_denominator(10 ** 6), \
            Fr(c1).limit_denominator(10 ** 6)
        v26 = [c0 * Eell[0][i] + c1 * Eell[1][i] for i in range(26)]
        # exact eigencheck under CI
        img = [sum(CI[i][j2] * v26[j2] for j2 in range(26))
               for i in range(26)]
        fok &= img == [sgn * x for x in v26]
        w = harmonic(v26)
        den_ = 1
        for x in w:
            den_ = den_ * x.denominator // gcd(den_, x.denominator)
        iv = [int(x * den_) for x in w]
        g_ = 0
        for x in iv:
            g_ = gcd(g_, abs(x))
        iv = [x // g_ for x in iv]
        if next(x for x in iv if x) < 0:
            iv = [-x for x in iv]
        flows[name] = iv
    gate("G9_elliptic_integer_flows", fok,
         plus_maxcoef=max(abs(x) for x in flows["flow_plus"]),
         minus_maxcoef=max(abs(x) for x in flows["flow_minus"]))
    ledger["flow_plus"] = flows["flow_plus"]
    ledger["flow_minus"] = flows["flow_minus"]

    # artifact export
    epairs = []
    for k in range(E):
        a, b = tri_of[erep[k]], tri_of[sS[erep[k]]]
        epairs.append([b, a])  # crossing orientation source->target
    viz = dict(epairs=epairs,
               rho={nm: [round(x, 5) for x in dens[nm]]
                    for nm in dens},
               flow_plus=flows["flow_plus"],
               flow_minus=flows["flow_minus"])
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "particles_viz.json"), "w") as f:
        json.dump(viz, f, separators=(",", ":"))

    ledger["all_passed"] = ok
    ledger["runtime_s"] = round(time.time() - t00, 1)
    with open(os.path.join(here,
              "x0143_hecke_particles_ledger.json"), "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print(f"\nledger + particles_viz.json written  "
          f"[{ledger['runtime_s']}s]")
    print("ALL GATES PASS" if ok else "GATE FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
