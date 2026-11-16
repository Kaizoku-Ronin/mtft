#!/usr/bin/env python3
"""
mp12_tg4_verify.py  --  Route A for MP12 (test T-G4 of the gravity map)

What replaces the Rademacher number at level 143?  For each closed line of
instrument 12.7 (a lift of a closed geodesic of the modular surface to
Gamma_0(143)\\PSL(2,R)) this script computes 30 integers:

    Psi                     the Rademacher number, the fibre coordinate
    phi_1, phi_A, phi_B     three windings round the cusps (eta quotients)
    2 + 4 + 8 + 12          periods on the four Hecke blocks of H_1(X_0(143))

Written from scratch: nothing is imported from the lab or from mtft.
Predictions P1 to P11 are in MP12_TG4_Plan.md and were written before the
first run; P12 and P13 were added there before their own runs.  The last
section (H1 to H5) is post hoc and is counted apart.  Route B:
mp12_tg4_mtft.py (the programme's package) and test_core_tg4.js (the
page's core).

Run:  python3 mp12_tg4_verify.py [mp12_numbers.json]      (about ten seconds)
Needs sympy, mpmath, numpy.

Outcome on 12 November 2026: 25 of 26 pre-registered checks, 6 of 6 post
hoc.  The one miss is the conjecture that closes P13, and it is a recorded
miss: the script exits with status 1 because of it.
"""
import sys, json, math, time, random
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
import sympy as sy

t_start = time.time()
random.seed(20261112)
N = 143
RESULTS, NUM = [], {}


def check(tag, name, cond, detail=""):
    RESULTS.append((tag, name, bool(cond), detail))
    print(("ok   " if cond else "FAIL ") + f"{tag:4s} {name}" + (f"   [{detail}]" if detail else ""), flush=True)


# ---------------------------------------------------------------- 2 x 2 integer matrices
def mm(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))


def inv(A):
    return ((A[1][1], -A[0][1]), (-A[1][0], A[0][0]))


def mpow(A, n):
    R_, B = ((1, 0), (0, 1)), A
    while n:
        if n & 1:
            R_ = mm(R_, B)
        B = mm(B, B)
        n >>= 1
    return R_


I2 = ((1, 0), (0, 1))
S = ((0, -1), (1, 0))
T = ((1, 1), (0, 1))
Ti = ((1, -1), (0, 1))
R = ((0, 1), (-1, -1))            # R = S^-1 T, R^3 = 1, T = S R
R2 = mm(R, R)
assert mm(S, R) == T and mm(R, R2) == I2

# ---------------------------------------------------------------- A: the 168 darts and the graph
units = [u for u in range(1, N) if math.gcd(u, N) == 1]
CAN = {}
for c in range(N):
    for d in range(N):
        if math.gcd(math.gcd(c, d), N) == 1:
            CAN[(c, d)] = min(((u * c) % N, (u * d) % N) for u in units)
pts = sorted(set(CAN.values()))
idx = {p: i for i, p in enumerate(pts)}
ND = len(pts)


def act(i, g):
    c, d = pts[i]
    return idx[CAN[((c * g[0][0] + d * g[1][0]) % N, (c * g[0][1] + d * g[1][1]) % N)]]


pS = [act(i, S) for i in range(ND)]
pR = [act(i, R) for i in range(ND)]
pR2 = [act(i, R2) for i in range(ND)]
pT = [act(i, T) for i in range(ND)]
x0 = idx[CAN[(0, 1)]]


def orbits(p):
    seen, out = [False] * len(p), []
    for i in range(len(p)):
        if not seen[i]:
            o, j = [], i
            while not seen[j]:
                seen[j] = True; o.append(j); j = p[j]
            out.append(o)
    return out


verts, edges, faces = orbits(pR), orbits(pS), orbits(pT)
vof = {x: k for k, o in enumerate(verts) for x in o}
eof = {x: k for k, o in enumerate(edges) for x in o}
sgn = {x: (1 if x == min(o) else -1) for o in edges for x in o}
fof = {x: k for k, o in enumerate(faces) for x in o}
fpos = {x: j for o in faces for j, x in enumerate(o)}
NV, NE, NF = len(verts), len(edges), len(faces)
loops = [k for k, o in enumerate(edges) if vof[o[0]] == vof[o[1]]]
# spanning tree by breadth-first search on the vertices
tree, seenv, queue = set(), {vof[x0]}, [vof[x0]]
while queue:
    v = queue.pop(0)
    for x in verts[v]:
        w = vof[pS[x]]
        if w not in seenv:
            seenv.add(w); tree.add(eof[x]); queue.append(w)
nontree = [k for k in range(NE) if k not in tree]
ntpos = {k: j for j, k in enumerate(nontree)}
B1 = len(nontree)
genus = (2 - (NV - NE + NF)) // 2
check("P1", "56 vertices, 84 edges (one loop), 4 faces of lengths 1, 11, 13, 143; connected", (NV, NE, NF) == (56, 84, 4) and len(loops) == 1 and sorted(len(o) for o in faces) == [1, 11, 13, 143] and len(seenv) == NV and all(len(o) == 3 for o in verts) and all(len(o) == 2 for o in edges), f"V, E, F = {NV}, {NE}, {NF}; loops {len(loops)}; faces {sorted(len(o) for o in faces)}")
check("P1", "b1 = 29 and genus 13: H_1 of the punctured curve is Z^29, of the 3-fold rank 30", B1 == 29 and genus == 13, f"b1 = {B1}, genus {genus}")
NUM["P1"] = {"V": NV, "E": NE, "F": NF, "loops": len(loops), "faces": sorted(len(o) for o in faces), "b1": B1, "genus": genus}


# ---------------------------------------------------------------- B: words, paths and the class map
def word(g):
    """g = T^q0 S T^q1 S ... T^qn up to sign, as tokens ('T', q) and ('S',)."""
    (a, b), (c, d) = g
    out = []
    while c != 0:
        q = a // c
        a, b = a - q * c, b - q * d
        out.append(('T', q)); out.append(('S',))
        a, b, c, d = c, d, -a, -b
    out.append(('T', b * d))
    return out


def trace_path(tokens, x, vec):
    """Move the dart x along the tokens, adding the signed edge crossings to vec; return the last dart."""
    for tok in tokens:
        if tok[0] == 'S':
            vec[eof[x]] += sgn[x]; x = pS[x]
        else:
            k = tok[1]
            if k == 0:
                continue
            o = faces[fof[x]]; w = len(o); p = fpos[x]
            full, rem = divmod(abs(k), w)
            if full:
                s_ = full if k > 0 else -full
                for y in o:
                    vec[eof[y]] += s_ * sgn[y]
            if k > 0:
                for j in range(rem):
                    y = o[(p + j) % w]; vec[eof[y]] += sgn[y]
                x = o[(p + k) % w]
            else:
                for j in range(1, rem + 1):
                    y = o[(p - j) % w]; vec[eof[y]] -= sgn[y]
                x = o[(p + k) % w]
    return x


def cls(g, start=None):
    vec = [0] * NE
    x = x0 if start is None else start
    end = trace_path(word(g), x, vec)
    if end != x:
        raise ValueError("not a closed path: the element is not in the group")
    return tuple(vec[k] for k in nontree), vec


def is_cycle(vec):
    """A closed path gives a cycle of the graph: every vertex has as many exits as entries."""
    bal = [0] * NV
    for k, o in enumerate(edges):
        tail = o[0] if sgn[o[0]] == 1 else o[1]
        bal[vof[tail]] -= vec[k]; bal[vof[pS[tail]]] += vec[k]
    return all(v == 0 for v in bal)


# coset representatives along the tree: g_x with x0 . g_x = x, crossing tree edges only
grep = {x0: I2}
queue = [x0]
while queue:
    x = queue.pop(0)
    for y, g in ((pR[x], R), (pR2[x], R2)) + (((pS[x], S),) if eof[x] in tree else ()):
        if y not in grep:
            grep[y] = mm(grep[x], g); queue.append(y)
assert len(grep) == ND and all(act(x0, grep[x]) == x for x in grep)
gens = []
for k in nontree:
    o = edges[k]; tail = o[0] if sgn[o[0]] == 1 else o[1]
    g = mm(mm(grep[tail], S), inv(grep[pS[tail]]))
    assert g[1][0] % N == 0 and g[0][0] * g[1][1] - g[0][1] * g[1][0] == 1
    gens.append(g)
unit_ok = all(cls(g)[0] == tuple(1 if j == i else 0 for j in range(B1)) and is_cycle(cls(g)[1]) for i, g in enumerate(gens))


def rand_gamma(length):
    g = I2
    for _ in range(length):
        h = random.choice(gens)
        g = mm(g, h if random.random() < 0.5 else inv(h))
    return g


add = lambda u, v: tuple(a + b for a, b in zip(u, v))
hom_ok = conj_ok = True
for _ in range(200):
    g, h = rand_gamma(random.randint(1, 6)), rand_gamma(random.randint(1, 6))
    hom_ok &= cls(mm(g, h))[0] == add(cls(g)[0], cls(h)[0])
    conj_ok &= cls(mm(mm(h, g), inv(h)))[0] == cls(g)[0]
check("P2", "the 29 tree generators have unit classes; cls is additive and conjugation-invariant (200 random pairs)", unit_ok and hom_ok and conj_ok)
face_cls, par = {}, {}
for o in faces:
    w = len(o); x = min(o)
    p = mm(mm(grep[x], mpow(T, w)), inv(grep[x]))
    vec = [0] * NE
    for y in o:
        vec[eof[y]] += sgn[y]
    face_cls[w] = tuple(vec[k] for k in nontree); par[w] = p
face_ok = all(cls(par[w])[0] == face_cls[w] for w in par) and all(sum(face_cls[w][j] for w in face_cls) == 0 for j in range(B1))
check("P2", "the loop round each cusp has the face cycle as its class, and the four face cycles add up to zero", face_ok)
NUM["P2"] = {"generators": [[list(r) for r in g] for g in gens], "face_cycles": {str(w): list(v) for w, v in face_cls.items()}}


# ---------------------------------------------------------------- C: Dedekind sums, Rademacher, eta quotients
def dedekind(h, k):
    """s(h, k) for k > 0 and gcd(h, k) = 1, by the reciprocity law."""
    h %= k
    s, sg = Fr(0), 1
    while h:
        s += sg * ((Fr(h, k) + Fr(k, h) + Fr(1, h * k)) / 12 - Fr(1, 4))
        sg = -sg
        h, k = k % h, h
    return s


def dedekind_slow(h, k):
    saw = lambda x: Fr(0) if x.denominator == 1 else x - math.floor(x) - Fr(1, 2)
    return sum(saw(Fr(r, k)) * saw(Fr(h * r, k)) for r in range(1, k))


ded_ok = all(dedekind(h, k) == dedekind_slow(h, k) for k in range(1, 60) for h in range(1, k + 1) if math.gcd(h, k) == 1)
sign = lambda v: (v > 0) - (v < 0)


def Phi(g):
    (a, b), (c, d) = g
    if c == 0:
        return Fr(b, d)
    return Fr(a + d, c) - 12 * sign(c) * dedekind(d, abs(c))


def Psi(g):
    (a, b), (c, d) = g
    return Phi(g) - 3 * sign(c * (a + d))


UNITS = {"phi1": {1: 1, 11: -1, 13: -1, 143: 1}, "phiA": {1: 12, 11: -12}, "phiB": {1: 2, 13: -2}}


def winding(g, r):
    (a, b), (c, d) = g
    if c < 0 or (c == 0 and d < 0):
        a, b, c, d = -a, -b, -c, -d
    tot = Fr(0)
    for dl, e in r.items():
        tot += e * (Fr(dl * b, d) if c == 0 else Fr(dl * (a + d), c) - 12 * dedekind(d, c // dl))
    return tot / 24


def phis(g):
    return tuple(winding(g, UNITS[k]) for k in ("phi1", "phiA", "phiB"))


# the units are functions on the group (eta by Euler's pentagonal series at 160 digits)
mp.mp.dps = 160


def eta(tau):
    tau = mp.mpc(tau)
    q = mp.exp(2j * mp.pi * tau)
    K = int(mp.sqrt(170 * mp.log(10) / (3 * mp.pi * tau.imag))) + 3
    s = mp.mpc(1)
    for k in range(1, K + 1):
        s += (-1) ** k * (q ** (k * (3 * k - 1) // 2) + q ** (k * (3 * k + 1) // 2))
    return mp.exp(mp.pi * 1j * tau / 12) * s


def uval(r, tau):
    v = mp.mpc(1)
    for dl, e in r.items():
        v *= eta(dl * tau) ** e
    return v


def small_gamma():
    cp = random.choice([-3, -2, -1, 1, 2, 3]); c = N * cp
    while True:
        d = random.randint(-60, 60)
        if d != 0 and math.gcd(c, d) == 1:
            break
    a = pow(d, -1, abs(c)); b = (a * d - 1) // c
    return ((a, b), (c, d))


winv = mp.mpf(0)
for _ in range(30):
    (a, b), (c, d) = small_gamma()
    tau = mp.mpf(-d) / c + mp.mpc(random.uniform(-0.5, 0.5), random.uniform(0.6, 1.6)) / abs(c)
    for r in UNITS.values():
        winv = max(winv, abs(uval(r, (a * tau + b) / (c * tau + d)) / uval(r, tau) - 1))
check("P3", "Dedekind sums by reciprocity = by definition (k < 60); the three units are unchanged by 30 random elements", ded_ok and winv < mp.mpf(10) ** -20, f"worst {float(winv):.1e}")
PHI_GEN = [phis(g) for g in gens]
int_ok = all(v.denominator == 1 for row in PHI_GEN for v in row)
lin_ok = True
for _ in range(300):
    g = rand_gamma(random.randint(1, 8)); c_ = cls(g)[0]
    lin_ok &= phis(g) == tuple(sum(c_[j] * PHI_GEN[j][i] for j in range(B1)) for i in range(3))
check("P3", "the windings are integers and integer linear forms of cls (fitted on the generators, 300 random elements)", int_ok and lin_ok)
cusp_tab = {w: tuple(int(v) for v in phis(par[w])) for w in (1, 143, 13, 11)}
want_tab = {1: (5, -5, -1), 143: (5, 65, 11), 13: (-5, -65, 1), 11: (-5, 5, -11)}
check("P3", "on the cusp loops of widths 1, 143, 13, 11: phi_1 = 5, 5, -5, -5; phi_A = -5, 65, -65, 5; phi_B = -1, 11, 1, -11", cusp_tab == want_tab, str(cusp_tab))

# the same integers by integrating (1/24) sum delta r_delta E2(delta z) dz from tau to gamma tau
NT = 6000
sig1 = np.zeros(NT + 1)
for dd in range(1, NT + 1):
    sig1[dd::dd] += dd


def E2np(z):
    q = np.exp(2j * np.pi * z)
    nt = min(NT, int(40 / (2 * np.pi * z.imag.min())) + 50)            # terms until |q|^n < e^-40
    qn = np.cumprod(np.repeat(q[:, None], nt, axis=1), axis=1)
    return 1 - 24 * np.sum(sig1[1:nt + 1][None, :] * qn, axis=1)


def winding_numeric(g, r):
    (a, b), (c, d) = g
    z0, z1 = complex(-d, 1) / c, complex(a, 1) / c
    panels = max(8, int(abs(z1 - z0) * abs(c) * 2))
    xg, wg = np.polynomial.legendre.leggauss(24)
    tot = 0j
    for p in range(panels):
        lo, hi = z0 + (z1 - z0) * p / panels, z0 + (z1 - z0) * (p + 1) / panels
        z = lo + (hi - lo) * (xg + 1) / 2
        f = sum(e * dl * E2np(dl * z) for dl, e in r.items())
        tot += np.sum(wg * f) * (hi - lo) / 2
    return tot / 24


SMALL = [((1, 0), (143, 1)), ((12, 1), (143, 12)), ((84, -47), (143, -80)), ((1, 0), (286, 1)), ((72, 1), (143, 2)), ((191, 2), (286, 3))]
assert all(g[0][0] * g[1][1] - g[0][1] * g[1][0] == 1 and g[1][0] % N == 0 for g in SMALL)
wnum = 0.0
for g in SMALL:
    for k_, r in UNITS.items():
        wnum = max(wnum, abs(winding_numeric(g, r) - float(winding(g, r))))
check("P3", "the integral of (1/24) sum delta r_delta E2(delta z) dz gives the same integers (six small elements)", wnum < 1e-8, f"worst {wnum:.1e}; exact values " + "; ".join(str(tuple(int(v) for v in phis(g))) for g in SMALL))
NUM["P3"] = {"phi_on_generators": [[int(v) for v in row] for row in PHI_GEN], "cusp_loops": {str(w): list(v) for w, v in cusp_tab.items()}, "numeric_err": wnum}


# ---------------------------------------------------------------- D: the closed lines of 12.7
def positive_word(A):
    """Conjugate A to a matrix with entries >= 0 and peel it into letters U = T, V = [[1,0],[1,1]]; return (#U, #V)."""
    frontier, seen = [A], {A}
    found = None
    for _ in range(12):
        nxt = []
        for B in frontier:
            for B_ in (B, ((-B[0][0], -B[0][1]), (-B[1][0], -B[1][1]))):
                if min(B_[0] + B_[1]) >= 0:
                    found = B_
            if found:
                break
            for h in (S, T, Ti):
                C = mm(mm(h, B), inv(h))
                if C not in seen:
                    seen.add(C); nxt.append(C)
        if found:
            break
        frontier = nxt
    (a, b), (c, d) = found
    nu = nv = 0
    while (a, b, c, d) != (1, 0, 0, 1):
        if a >= c and b >= d:
            a, b = a - c, b - d; nu += 1
        else:
            c, d = c - a, d - b; nv += 1
    return nu, nv


def chi(g):
    return sum(t[1] if t[0] == 'T' else 3 for t in word(g)) % 6


MOTIONS = {"cat": ((2, 1), (1, 1)), "s10": ((84, -47), (143, -80)), "pell": ((12, 1), (143, 12))}
PSI_A = {k: Psi(A) for k, A in MOTIONS.items()}
letters = {k: positive_word(A) for k, A in MOTIONS.items()}
check("P4", "Psi of the three motions is 0, -1, -21, and equals #U - #V of the positive word", [PSI_A[k] for k in ("cat", "s10", "pell")] == [0, -1, -21] and all(PSI_A[k] == letters[k][0] - letters[k][1] for k in MOTIONS), f"letters (#U, #V): {letters}")
LINES = []
for name, A in MOTIONS.items():
    pA = [act(i, A) for i in range(ND)]
    wA = word(A)
    for o in sorted(orbits(pA), key=lambda o: (len(o), min(o))):
        x, ell = min(o), len(o)
        vec = [0] * NE
        y = x
        for _ in range(ell):
            y = trace_path(wA, y, vec)
        assert y == x
        gam = mm(mm(grep[x], mpow(A, ell)), inv(grep[x]))
        assert gam[1][0] % N == 0
        c_direct = cls(gam)[0]
        LINES.append({"motion": name, "ell": ell, "dart": x, "cycle": frozenset(o), "gamma": gam, "cls": tuple(vec[k] for k in nontree), "cls_direct": c_direct,
                      "psi": Psi(gam), "phi": phis(gam), "chi": chi(gam)})
types = {name: sorted(l["ell"] for l in LINES if l["motion"] == name) for name in MOTIONS}
check("P4", "28 closed lines: cycles 7^4 35^4, 1^4 5^4 6^4 30^4, 1 11 13 143; tracing the cycle and decomposing the large matrix give the same class", types == {"cat": [7] * 4 + [35] * 4, "s10": [1] * 4 + [5] * 4 + [6] * 4 + [30] * 4, "pell": [1, 11, 13, 143]} and all(l["cls"] == l["cls_direct"] and is_cycle([0] * NE) for l in LINES), f"largest entry has {max(len(str(abs(v))) for l in LINES for row in l['gamma'] for v in row)} digits")
check("P4", "Psi of a lift is ell times Psi of the motion", all(l["psi"] == l["ell"] * PSI_A[l["motion"]] for l in LINES))
CHI_GEN = [chi(g) for g in gens]
chi_lin = lambda c_: sum(c_[j] * CHI_GEN[j] for j in range(B1)) % 6
cong_ok = all(int(l["psi"]) % 6 == l["chi"] == chi_lin(l["cls"]) for l in LINES)
check("P4", "Psi = chi (mod 6) on the 28 lines, and chi is a linear form of cls modulo 6", cong_ok, f"chi on the generators: {CHI_GEN}")
tr_ok = all(all(sum(l["cls"][j] for l in LINES if l["motion"] == name) == 0 for j in range(B1)) for name in MOTIONS)
fib = {name: int(sum(l["psi"] for l in LINES if l["motion"] == name)) for name in MOTIONS}
check("P5", "transfer: for each motion the classes of all its lifts add up to zero; the fibre numbers add up to 168 Psi(A)", tr_ok and all(fib[k] == 168 * PSI_A[k] for k in MOTIONS), f"sums of Psi: {fib}")
cat = [l for l in LINES if l["motion"] == "cat"]
pairs, self_paired = [], 0
for l in cat:
    partner = frozenset(pS[x] for x in l["cycle"])
    m_ = next(k for k in cat if k["cycle"] == partner)
    if m_ is l:
        self_paired += 1
    pairs.append((l["ell"], m_["cls"] == tuple(-v for v in l["cls"]), m_ is not l))
check("P6", "reciprocity: each cat lift and its partner through S have opposite classes; no lift is its own partner", all(p[1] and p[2] for p in pairs) and self_paired == 0, f"{len(pairs) // 2} pairs")
zero_lines = [(l["motion"], l["ell"]) for l in LINES if not any(l["cls"])]
NUM["P4"] = {"psi_motions": {k: int(v) for k, v in PSI_A.items()}, "letters": letters, "chi_generators": CHI_GEN}


# ---------------------------------------------------------------- E: Manin symbols, T_2 and the blocks
def merel(p):
    out = []
    for a in range(1, p + 1):
        for b in range(0, a):
            for c in range(0, p + 1):
                if (p + b * c) % a == 0:
                    d = (p + b * c) // a
                    if d > c and a * d - b * c == p:
                        out.append(((a, b), (c, d)))
    return out


tail = [o[0] if sgn[o[0]] == 1 else o[1] for o in edges]              # the dart whose Manin symbol is +E_k
pTS = [pS[pT[i]] for i in range(ND)]                                    # x -> x.(TS), of order 3
rel = []
for o in orbits(pTS):
    row = [0] * NE
    for x in o:
        row[eof[x]] += sgn[x]
    rel.append(row)
RelM = sy.Matrix(rel)
rref, piv = RelM.rref()
free = [k for k in range(NE) if k not in piv]
Q = sy.zeros(len(free), NE)                                             # quotient coordinates: coefficients on the free edges
fpos_ = {k: j for j, k in enumerate(free)}
for k in free:
    Q[fpos_[k], k] = 1
for i, p_ in enumerate(piv):
    for k in free:
        Q[fpos_[k], p_] = -rref[i, k]
DQ = sy.zeros(NF, len(free))                                            # the boundary: cusp of g(infinity) minus cusp of g(0)
Dfull = sy.zeros(NF, NE)
for k in range(NE):
    Dfull[fof[tail[k]], k] += 1; Dfull[fof[pS[tail[k]]], k] -= 1
for j, k in enumerate(free):
    DQ[:, j] = Dfull[:, k]
bd_ok = (Dfull * RelM.T).is_zero_matrix
Kb = DQ.nullspace()
Kmat = sy.Matrix.hstack(*Kb)                                            # 29 x 26, a basis of the cuspidal part


def hecke_q(p):
    Tm = sy.zeros(NE, NE)
    for k in range(NE):
        c, d = pts[tail[k]]
        for M in merel(p):
            u, v = (c * M[0][0] + d * M[1][0]) % N, (c * M[0][1] + d * M[1][1]) % N
            if math.gcd(math.gcd(u, v), N) == 1:
                y = idx[CAN[(u, v)]]
                Tm[eof[y], k] += sgn[y]
    well = (Q * Tm * RelM.T).is_zero_matrix
    J = sy.zeros(NE, len(free))
    for j, k in enumerate(free):
        J[k, j] = 1
    return Q * Tm * J, well


def to_cusp(Tq):                                                        # the matrix of an operator on the cuspidal part, in the basis Kb
    return (Kmat.T * Kmat).inv() * Kmat.T * Tq * Kmat


T2q, w2 = hecke_q(2)
T3q, w3 = hecke_q(3)
T2c, T3c = to_cusp(T2q), to_cusp(T3q)
x_ = sy.symbols('x')
cp2 = sy.factor_list(T2c.charpoly(x_).as_expr())
cpq = sy.factor(T2q.charpoly(x_).as_expr())
g4 = x_ ** 4 - 3 * x_ ** 3 - x_ ** 2 + 5 * x_ + 1
h6 = x_ ** 6 - 10 * x_ ** 4 + 2 * x_ ** 3 + 24 * x_ ** 2 - 7 * x_ - 12
fac = {sy.expand(f): m for f, m in cp2[1]}
check("P8", "Manin model: dimension 29, cuspidal 26; relations and T_2, T_3 respect the boundary and the relations; T_2 T_3 = T_3 T_2", len(free) == 29 and len(Kb) == 26 and bd_ok and w2 and w3 and (T2q * T3q - T3q * T2q).is_zero_matrix and (DQ * T2q * Kmat).is_zero_matrix)
check("P8", "charpoly of T_2 on the cuspidal part = x^2 (x + 2)^4 g4^2 h6^2; eigenvalue 3 three times on the rest", fac == {x_: 2, x_ + 2: 4, sy.expand(g4): 2, sy.expand(h6): 2} and sy.Poly(cpq, x_).rem(sy.Poly((x_ - 3) ** 3, x_)).is_zero, str(sy.factor(T2c.charpoly(x_).as_expr())))
polys = {"ell": x_, "old": x_ + 2, "q4": g4, "q6": h6}
BLK = {}
for name, f in polys.items():
    co = sy.Poly(f, x_).all_coeffs()[::-1]
    Mf = sy.zeros(26, 26); P_ = sy.eye(26)
    for cf in co:
        Mf += cf * P_; P_ = P_ * T2c
    BLK[name] = Mf.nullspace()
dims = {k: len(v) for k, v in BLK.items()}
Wc = sy.Matrix.hstack(*[v for name in ("ell", "old", "q4", "q6") for v in BLK[name]])      # 26 x 26: block bases in cuspidal coordinates
Winv = Wc.inv()
# a_p on the 2-dimensional block against the curve y^2 + y = x^3 - x^2 - x - 2
Bell = sy.Matrix.hstack(*BLK["ell"])


def ap_curve(p):
    if p == 2:
        return 2 + 1 - (1 + sum(1 for x in range(2) for y in range(2) if (y * y + y - (x ** 3 - x * x - x - 2)) % 2 == 0))
    s = 0
    for x in range(p):
        v = (4 * x ** 3 - 4 * x * x - 4 * x - 7) % p
        if v:
            s += 1 if pow(v, (p - 1) // 2, p) == 1 else -1
    return -s


ap_blk = {}
for p in (2, 3, 5, 7):
    Tp = to_cusp(hecke_q(p)[0])
    img = Tp * Bell
    lam = [sy.Rational(img[i, 0], Bell[i, 0]) for i in range(26) if Bell[i, 0] != 0][0]
    ap_blk[p] = lam if (img - lam * Bell).is_zero_matrix else None
check("P8", "block dimensions 2, 4, 8, 12; on the 2-dimensional block T_p = a_p of y^2 + y = x^3 - x^2 - x - 2 for p = 2, 3, 5, 7", dims == {"ell": 2, "old": 4, "q4": 8, "q6": 12} and all(ap_blk[p] == ap_curve(p) for p in (2, 3, 5, 7)) and [ap_curve(p) for p in (2, 3, 5, 7)] == [0, -1, -1, -2], f"block: {ap_blk}")


def manin(g):
    """The modular symbol {infinity, g(infinity)} in the 29 quotient coordinates."""
    (a, b), (c, d) = g
    if c < 0:
        a, c = -a, -c
    vec = [0] * NE
    if c != 0:
        pm, qm, p_, q_ = 1, 0, a // c, 1                                # convergents p_{-1}/q_{-1}, p_0/q_0
        num, den = a - (a // c) * c, c
        y = idx[CAN[((-1) % N, 0)]]; vec[eof[y]] += sgn[y]              # j = 0: the symbol (-q_0 : q_{-1}) = (-1 : 0)
        j = 0
        while num:
            j += 1
            aj = den // num
            num, den = den - aj * num, num
            pm, qm, p_, q_ = p_, q_, aj * p_ + pm, aj * q_ + qm
            y = idx[CAN[(((-1) ** (j - 1) * q_) % N, qm % N)]]; vec[eof[y]] += sgn[y]
        assert (p_, q_) == (a, c)
    return Q * sy.Matrix(vec)


Mgen = sy.Matrix.hstack(*[manin(g) for g in gens])                     # 29 x 29: the Manin class of each generator
cusp_ok = (DQ * Mgen).is_zero_matrix
face_zero = all((Mgen * sy.Matrix(face_cls[w])).is_zero_matrix for w in face_cls)
direct_ok = all((manin(l["gamma"]) - Mgen * sy.Matrix(l["cls"])).is_zero_matrix for l in LINES)
check("P8", "every class is cuspidal, the cusp loops map to zero, and the symbol of each lift from its large matrix = the linear image of its cls", cusp_ok and face_zero and direct_ok and Mgen.rank() == 26)
# coordinates in the block bases, then integer coordinates in a Z-basis of each block's projected lattice
KtK = (Kmat.T * Kmat).inv() * Kmat.T                                    # quotient coordinates -> cuspidal coordinates
ALPHA = (Winv * KtK * Mgen).T                                           # 29 x 26: row e = block coefficients of generator e


def lattice_basis(rows):
    A = [list(r) for r in rows]
    m, n = len(A), len(A[0])
    r = 0
    for c in range(n):
        while True:
            nz = [i for i in range(r, m) if A[i][c] != 0]
            if len(nz) <= 1:
                break
            piv_ = min(nz, key=lambda i: abs(A[i][c]))
            for i in nz:
                if i != piv_:
                    q = A[i][c] // A[piv_][c]
                    A[i] = [u - q * v for u, v in zip(A[i], A[piv_])]
        nz = [i for i in range(r, m) if A[i][c] != 0]
        if nz:
            A[r], A[nz[0]] = A[nz[0]], A[r]
            if A[r][c] < 0:
                A[r] = [-u for u in A[r]]
            r += 1
            if r == m:
                break
    return A[:r]


COORD, off = {}, 0
for name in ("ell", "old", "q4", "q6"):
    dm = dims[name]
    sub = ALPHA[:, off:off + dm]; off += dm
    den = sy.ilcm(*[sy.Rational(v).q for v in sub])
    basis = sy.Matrix(lattice_basis([[int(v * den) for v in sub.row(i)] for i in range(B1)])) / den      # dm x dm
    C = sub * basis.inv()                                                                                 # 29 x dm, integers
    assert basis.shape == (dm, dm) and all(sy.Rational(v).q == 1 for v in C)
    COORD[name] = [[int(v) for v in C.row(i)] for i in range(B1)]


def block_coords(c_):
    return {name: [sum(c_[e] * COORD[name][e][i] for e in range(B1)) for i in range(dims[name])] for name in COORD}


BIG = sy.Matrix([[int(v) for v in PHI_GEN[e]] + [v for name in ("ell", "old", "q4", "q6") for v in COORD[name][e]] for e in range(B1)])
index_all = abs(BIG.det())
CUSPM = sy.Matrix([[v for name in ("ell", "old", "q4", "q6") for v in COORD[name][e]] for e in range(B1)])
for l in LINES:
    l["blocks"] = block_coords(l["cls"])
tr2_ok = all(all(sum(l["blocks"][name][i] for l in LINES if l["motion"] == m_) == 0 for name in COORD for i in range(dims[name])) and all(sum(l["phi"][i] for l in LINES if l["motion"] == m_) == 0 for i in range(3)) for m_ in MOTIONS)
check("P5", "so do their windings and their periods, block by block", tr2_ok)
NUM["P8"] = {"charpoly_T2": str(sy.factor(T2c.charpoly(x_).as_expr())), "a_p_block": {str(p): int(v) for p, v in ap_blk.items()}, "index_product_lattice": int(index_all), "index_factors": {str(k): int(v) for k, v in sy.factorint(int(index_all)).items()}}

# ---------------------------------------------------------------- P7: the Fricke involution and the Pell line
pell1 = next(l for l in LINES if l["motion"] == "pell" and l["ell"] == 1)
gP = MOTIONS["pell"]
fricke = ((gP[1][1], -gP[1][0] // N), (-N * gP[0][1], gP[0][0]))       # W gamma W^-1 for gamma in Gamma_0(N)
check("P7", "W_143 reverses the Pell unit; its 143a1 coordinates and its winding phi_1 are zero", fricke == inv(gP) and pell1["gamma"] in (gP, inv(gP)) and pell1["blocks"]["ell"] == [0, 0] and pell1["phi"][0] == 0, f"ell {pell1['blocks']['ell']}, phi {tuple(int(v) for v in pell1['phi'])}")

# ---------------------------------------------------------------- F: the periods of 143a1, numerically
PMAX = 6000
sieve = bytearray([1]) * (PMAX + 1); sieve[0:2] = b"\x00\x00"
for i in range(2, int(PMAX ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
an = [0, 1] + [0] * (PMAX - 1)
for p in range(2, PMAX + 1):
    if sieve[p]:
        ap = ap_curve(p); pk, prev, cur = p, 1, ap
        while pk <= PMAX:
            an[pk] = cur
            prev, cur = cur, (ap * cur - (0 if p in (11, 13) else p) * prev)
            pk *= p
for n in range(2, PMAX + 1):
    if an[n] == 0 and not sieve[n]:
        m, p = n, next(p for p in range(2, n) if n % p == 0)
        pk = 1
        while m % p == 0:
            m //= p; pk *= p
        if m > 1:
            an[n] = an[pk] * an[m]
anf = np.array(an[1:], dtype=float) / np.arange(1, PMAX + 1)
narr = np.arange(1, PMAX + 1)


def period(g):
    (a, b), (c, d) = g
    z0, z1 = complex(-d, 1) / c, complex(a, 1) / c
    return complex(np.sum(anf * (np.exp(2j * np.pi * narr * z1) - np.exp(2j * np.pi * narr * z0))))


# First run (12 November 2026, about 09:19): with the plan's six small elements the 143a1 coordinates had rank 1
# ((0,0) for the two cusp loops and the Pell unit, multiples of one vector for the other three), so the two periods
# could not be solved for and the script stopped here.  A fault of the sample, recorded in MP12; six more elements
# with c = 143 are added (d = 4, 5, 6, 7, 9, 10) and the best pair among the twelve is used.
SMALL9 = list(SMALL)
for d in (4, 5, 6, 7, 9, 10):
    a = pow(d, -1, N); SMALL9.append(((a, (a * d - 1) // N), (N, d)))
PER = [period(g) for g in SMALL9]
ELL_SMALL = [block_coords(cls(g)[0])["ell"] for g in SMALL9]
n9 = len(SMALL9)
best = max(((i, j) for i in range(n9) for j in range(i + 1, n9)), key=lambda ij: abs(ELL_SMALL[ij[0]][0] * ELL_SMALL[ij[1]][1] - ELL_SMALL[ij[0]][1] * ELL_SMALL[ij[1]][0]))
(i_, j_) = best
detb = ELL_SMALL[i_][0] * ELL_SMALL[j_][1] - ELL_SMALL[i_][1] * ELL_SMALL[j_][0]
Om1 = (PER[i_] * ELL_SMALL[j_][1] - PER[j_] * ELL_SMALL[i_][1]) / detb
Om2 = (-PER[i_] * ELL_SMALL[j_][0] + PER[j_] * ELL_SMALL[i_][0]) / detb
perr = max(abs(PER[k] - (ELL_SMALL[k][0] * Om1 + ELL_SMALL[k][1] * Om2)) for k in range(n9))
tauL = Om2 / Om1 if (Om2 / Om1).imag > 0 else Om1 / Om2
mp.mp.dps = 30
jL = 1728 * mp.kleinj(mp.mpc(tauL))
c4, Delta = 64, -1859                                                   # of y^2 + y = x^3 - x^2 - x - 2
jE = mp.mpf(c4) ** 3 / Delta
check("P7", "the integral of 143a1's form along the Pell unit vanishes", abs(PER[1]) < 1e-10, f"|period| = {abs(PER[1]):.1e}")
check("P9", "periods of twelve small elements = n1 Omega1 + n2 Omega2 with the block's integers; the lattice has the curve's j-invariant", perr < 1e-8 and abs(jL - jE) / abs(jE) < 1e-6 and detb != 0, f"worst {perr:.1e}; (n1, n2) = {ELL_SMALL}; j = {complex(jL).real:.6f} against {float(jE):.6f}; Omega1 = {Om1:.6f}, Omega2 = {Om2:.6f}")
NUM["P9"] = {"ell_small": ELL_SMALL, "Omega1": [Om1.real, Om1.imag], "Omega2": [Om2.real, Om2.imag], "j_lattice": float(complex(jL).real), "j_curve": float(jE), "err": perr}

# ---------------------------------------------------------------- P12: why the ghost coordinates of the Pell lines vanish
NQ11 = 4000
pent = np.zeros(NQ11 + 1)
k = 0
while k * (3 * k - 1) // 2 <= NQ11:
    for e in (k * (3 * k - 1) // 2, k * (3 * k + 1) // 2):
        if e <= NQ11:
            pent[e] = (-1) ** k
    k += 1
sq = np.convolve(pent, pent)[:NQ11 + 1]                                 # prod (1 - q^n)^2
sq11 = np.zeros(NQ11 + 1); sq11[::11] = sq[:len(sq11[::11])]            # prod (1 - q^11n)^2
f11 = np.rint(np.convolve(sq, sq11)[:NQ11]).astype(np.int64)            # eta(tau)^2 eta(11 tau)^2 = q * this
a11 = np.concatenate([[0], f11])[:NQ11 + 1]                              # a_n of 11a1, a_1 = 1
assert list(a11[1:12]) == [1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1]
a11f = a11[1:] / np.arange(1, NQ11 + 1)
n11 = np.arange(1, NQ11 + 1)


def period11(g):
    (a, b), (c, d) = g
    z0, z1 = complex(-d, 1) / c, complex(a, 1) / c
    return complex(np.sum(a11f * (np.exp(2j * np.pi * n11 * z1) - np.exp(2j * np.pi * n11 * z0))))


img13 = lambda g: ((g[0][0], 13 * g[0][1]), (g[1][0] // 13, g[1][1]))
Aprime = img13(MOTIONS["pell"])
p11 = [(abs(period11(g)), abs(period11(img13(g)))) for g in SMALL9]
ghost9 = [block_coords(cls(g)[0])["old"] for g in SMALL9]
equiv = all((not any(gh)) == (pa < 1e-9 and pb < 1e-9) for gh, (pa, pb) in zip(ghost9, p11))
nz_min = min([max(pa, pb) for gh, (pa, pb) in zip(ghost9, p11) if any(gh)] or [0.0])
check("P12", "11a1's form integrates to zero along the Pell unit and along its image [[12, 13], [11, 12]]; over the twelve small elements the ghost coordinates vanish exactly when both integrals do", Aprime == ((12, 13), (11, 12)) and p11[1][0] < 1e-10 and p11[1][1] < 1e-10 and equiv and sum(1 for gh in ghost9 if any(gh)) >= 2, f"Pell: {p11[1][0]:.1e}, {p11[1][1]:.1e}; {sum(1 for gh in ghost9 if any(gh))} elements with ghost coordinates, smallest of their larger integrals {nz_min:.4f}")
NUM["P12"] = {"pell": list(p11[1]), "ghost_small": ghost9}

# ---------------------------------------------------------------- P13: the Atkin-Lehner involutions on the classes
WQ = {11: ((66, 5), (143, 11)), 13: ((78, 7), (143, 13)), 143: ((0, -1), (143, 0))}


def al_conj(Wm, g, Qd):
    adj = ((Wm[1][1], -Wm[0][1]), (-Wm[1][0], Wm[0][0]))
    h = mm(mm(Wm, g), adj)
    assert all(v % Qd == 0 for row in h for v in row)
    return tuple(tuple(v // Qd for v in row) for row in h)


assert all(Wm[0][0] * Wm[1][1] - Wm[0][1] * Wm[1][0] == Qd for Qd, Wm in WQ.items())
AL = {Qd: sy.Matrix([list(cls(al_conj(Wm, g, Qd))[0]) for g in gens]) for Qd, Wm in WQ.items()}      # row e = class of W gamma_e W^-1
I29 = sy.eye(B1)
al_ok = all((AL[Qd] * AL[Qd] - I29).is_zero_matrix for Qd in AL) and (AL[11] * AL[13] - AL[13] * AL[11]).is_zero_matrix and (AL[11] * AL[13] - AL[143]).is_zero_matrix
BLKM = {name: sy.Matrix(COORD[name]) for name in COORD}                 # 29 x dim: block coordinates of the generators


def al_sign(Qd, name):
    img = AL[Qd] * BLKM[name]                                           # block coordinates of the images
    for sg in (1, -1):
        if (img - sg * BLKM[name]).is_zero_matrix:
            return sg
    return 0


signs = {name: (al_sign(11, name), al_sign(13, name), al_sign(143, name)) for name in COORD}
pell_neg = (sy.Matrix([list(pell1["cls"])]) * AL[143] + sy.Matrix([list(pell1["cls"])])).is_zero_matrix
check("P13", "W_11, W_13, W_143 are commuting involutions of Z^29 with W_11 W_13 = W_143; signs on the blocks: 143a1 (+,+,+), f2 (-,+,-), f3 (+,-,-), ghost W_11 = -1; W_143 negates the Pell unit's class", al_ok and signs["ell"] == (1, 1, 1) and signs["q4"] == (-1, 1, -1) and signs["q6"] == (1, -1, -1) and signs["old"][0] == -1 and pell_neg, f"signs (W11, W13, W143): {signs}")
for l in LINES:
    row = sy.Matrix([list(l["cls"])])
    cusp = lambda r: [sum(r[e] * COORD[name][e][i] for e in range(B1)) for name in ("ell", "old", "q4", "q6") for i in range(dims[name])]
    me = cusp(list(row))
    l["al"] = {}
    for Qd in (11, 13, 143):
        im = cusp(list(row * AL[Qd]))
        l["al"][Qd] = -1 if im == [-v for v in me] else (1 if im == me else 0)
zero_ell = [l for l in LINES if not any(l["blocks"]["ell"])]
conj_d = all(-1 in l["al"].values() for l in zero_ell)
check("P13", "conjecture: each line with no 143a1 component is negated in H_1 of the closed curve by at least one of the three involutions", conj_d, "; ".join(f"{l['motion']} {l['ell']}: {l['al']}" for l in zero_ell))
NUM["P13"] = {"block_signs": {k: list(v) for k, v in signs.items()}, "lines": [[l["motion"], l["ell"], l["al"][11], l["al"][13], l["al"][143]] for l in LINES]}

# ---------------------------------------------------------------- P10: independence of the 30 integers
samp = []
while len(samp) < 60:
    g = rand_gamma(random.randint(2, 9))
    if abs(g[0][0] + g[1][1]) > 2:
        samp.append(g)
vecs = [[int(Psi(g))] + list(cls(g)[0]) for g in samp]
cong2 = all(v[0] % 6 == chi(g) == chi_lin(v[1:]) for v, g in zip(vecs, samp))
hb = lattice_basis(vecs)
rank30 = len(hb)
idx30 = abs(sy.Matrix(hb).det()) if rank30 == 30 else None
check("P10", "60 random hyperbolic elements: (Psi; cls) has rank 30, obeys Psi = chi (mod 6), and generates a lattice of index 6 in Z^30", rank30 == 30 and cong2 and idx30 == 6, f"rank {rank30}, index {idx30}")
NUM["P10"] = {"rank": rank30, "index": int(idx30) if idx30 is not None else None}

# ---------------------------------------------------------------- found, not predicted
rank28 = sy.Matrix([list(l["cls"]) for l in LINES]).rank()
rank28f = sy.Matrix([[int(l["psi"])] + list(l["cls"]) for l in LINES]).rank()
rk = {m_: sy.Matrix([list(l["cls"]) for l in LINES if l["motion"] == m_]).rank() for m_ in MOTIONS}
zero_blocks = {name: [(l["motion"], l["ell"]) for l in LINES if not any(l["blocks"][name])] for name in COORD}
zero_phi = [(l["motion"], l["ell"]) for l in LINES if not any(l["phi"])]
print("\nfound, not predicted")
print("  null-homologous lines (cls = 0):", zero_lines)
print("  rank of the 28 classes in Z^29:", rank28, "; with Psi in Z^30:", rank28f, "; by motion:", rk)
print("  index of H_1 of the punctured curve in the product of the block lattices and the three windings:", index_all, "=", dict(sy.factorint(int(index_all))))
idx_cusp = abs(sy.Matrix(lattice_basis([list(CUSPM.row(e)) for e in range(B1)])).det())
print("  index of H_1 of the closed curve in the product of the four block lattices:", idx_cusp, "=", dict(sy.factorint(int(idx_cusp))))
NUM["P8"]["index_cuspidal"] = int(idx_cusp); NUM["P8"]["index_cuspidal_factors"] = {str(k_): int(v_) for k_, v_ in sy.factorint(int(idx_cusp)).items()}
print("  lines with all three windings zero:", zero_phi)
for name in COORD:
    print(f"  lines with zero {name} component:", zero_blocks[name])
print("\n  motion ell   Psi   (phi_1, phi_A, phi_B)   143a1    ghost")
for l in LINES:
    print(f"  {l['motion']:5s} {l['ell']:3d} {int(l['psi']):6d}   {str(tuple(int(v) for v in l['phi'])):24s} {str(l['blocks']['ell']):10s} {l['blocks']['old']}")
NUM["found"] = {"zero_lines": zero_lines, "rank28": int(rank28), "rank28_with_psi": int(rank28f), "rank_by_motion": {k: int(v) for k, v in rk.items()}, "zero_blocks": zero_blocks, "zero_phi": zero_phi}
NUM["lines"] = [{"motion": l["motion"], "ell": l["ell"], "dart": list(pts[l["dart"]]), "psi": int(l["psi"]), "phi": [int(v) for v in l["phi"]], "cls": list(l["cls"]), "blocks": l["blocks"],
                 "gamma": [[str(v) for v in row] for row in l["gamma"]]} for l in LINES]
NUM["coord"] = {name: COORD[name] for name in COORD}
NUM["an_143a1"] = an[1:201]
NUM["small"] = [[list(r) for r in g] for g in SMALL9]

# ---------------------------------------------------------------- post hoc: after the fourth run, not predictions
# 12 November 2026, 09:39 New York time.  The fourth run missed P13's conjecture.  What follows was thought of after
# that miss, to learn why it missed and what holds in its place.  It is kept apart from the predictions: tags H1 to
# H5, a list and a count of their own.  Expected before this section first ran: H1, H2 and H4 as stated; in H3 the
# sign + for the cat lines of length 7 and - for Sphereland 10's lines of length 30; in H5 that W_11 keeps the four
# cat lines of length 7 among themselves.
POST = []


def post(tag, name, cond, detail=""):
    POST.append((tag, name, bool(cond), detail))
    print(("ok   " if cond else "FAIL ") + f"{tag:4s} {name}" + (f"   [{detail}]" if detail else ""), flush=True)


def p1q(v, q):
    c, d = v[0] % q, v[1] % q
    return (1, d * pow(c, -1, q) % q) if c else (0, 1)


def orbit_q(v, A, q):
    """The orbit of the dart v, read modulo q, under the motion A on P^1(F_q)."""
    o, w = [], (v[0] % q, v[1] % q)
    while p1q(w, q) not in o:
        o.append(p1q(w, q))
        w = ((w[0] * A[0][0] + w[1] * A[1][0]) % q, (w[0] * A[0][1] + w[1] * A[1][1]) % q)
    return frozenset(o)


def qform(v, A):
    """The motion's quadratic form on a dart: det[v; vA], unchanged by v -> vA and zero on the darts A fixes."""
    w = (v[0] * A[0][0] + v[1] * A[1][0], v[0] * A[0][1] + v[1] * A[1][1])
    return v[0] * w[1] - v[1] * w[0]


legendre = lambda a, q: 0 if a % q == 0 else (1 if pow(a % q, (q - 1) // 2, q) == 1 else -1)
neg = lambda c_: tuple(-v for v in c_)
app = lambda c_, M: tuple(sum(c_[e] * M[e][j] for e in range(B1)) for j in range(B1))      # the class c_ through a matrix whose row e is the image of generator e
print("\npost hoc (after the fourth run; not predictions)")
for l in LINES:
    A = MOTIONS[l["motion"]]; v = pts[l["dart"]]
    l["orb"] = {q: orbit_q(v, A, q) for q in (11, 13)}
    l["m"] = {q: len(l["orb"][q]) for q in (11, 13)}
    l["eps"] = {q: legendre(qform(v, A), q) for q in (11, 13)}
    assert l["m"][11] * l["m"][13] // math.gcd(l["m"][11], l["m"][13]) == l["ell"]

# H1: the transfer through X_0(13) and X_0(11)
names4, new3 = ("ell", "old", "q4", "q6"), ("ell", "q4", "q6")
bsum = lambda grp, b: [sum(l["blocks"][b][i] for l in grp) for i in range(dims[b])]
h1_ok, nrel, ghost_nz, n11 = True, 0, 0, 0
for name in MOTIONS:
    fam = [l for l in LINES if l["motion"] == name]
    for q in (11, 13):
        for o in set(l["orb"][q] for l in fam):
            grp = [l for l in fam if l["orb"][q] == o]
            nrel += 1
            h1_ok &= all(not any(bsum(grp, b)) for b in (new3 if q == 11 else names4))
            if q == 11:
                n11 += 1; ghost_nz += any(bsum(grp, "old"))
post("H1", "the finer transfer: the lines over one orbit modulo 13 add up to zero in all four blocks (X_0(13) has genus 0); over one orbit modulo 11 they add up to zero in the three new blocks (a newform has trace zero)", h1_ok and nrel == 18, f"{nrel} relations; over {ghost_nz} of the {n11} orbits modulo 11 the ghost components do not add up to zero")


# H2: where an Atkin-Lehner involution sends a closed line
def disc0(g):
    (a, b), (c, d) = g
    k = math.gcd(math.gcd(abs(c), abs(d - a)), abs(b))
    return ((a + d) ** 2 - 4) // (k * k)


DA = {k: disc0(A) for k, A in MOTIONS.items()}
h2_ok = DA == {"cat": 5, "s10": 12, "pell": 572} and all(disc0(l["gamma"]) == DA[l["motion"]] for l in LINES)
kept = {11: 0, 13: 0, 143: 0}
for l in LINES:
    l["cond"] = {}
    for Qd, Wm in WQ.items():
        f = 1
        for q in (11, 13):
            if Qd % q == 0 and l["m"][q] > 1:
                f *= q
        l["cond"][Qd] = f
        h2_ok &= disc0(al_conj(Wm, l["gamma"], Qd)) == DA[l["motion"]] * f * f
        kept[Qd] += f == 1
post("H2", "W_Q carries a lift of a motion of discriminant D to a closed line of discriminant D f^2, f the product of the primes q | Q at which the motion moves the dart: it keeps a line among the lifts of the same order only where the motion fixes the Q-part of the dart", h2_ok, f"D = {DA}; lines kept by W_11, W_13, W_143: {kept[11]}, {kept[13]}, {kept[143]} of 28")

# H3: the mirror tau -> -conj(tau), conjugation by diag(1, -1)
cj = lambda g: ((g[0][0], -g[0][1]), (-g[1][0], g[1][1]))
CJ = [list(cls(cj(g))[0]) for g in gens]                                # row e = class of the mirror image of generator e
CJm = sy.Matrix(CJ)
h3_gen = (CJm * CJm - I29).is_zero_matrix and all(phis(cj(g)) == neg(phis(g)) and Psi(cj(g)) == -Psi(g) for g in gens)
h3_ell = (CJm * BLKM["ell"] - BLKM["ell"] * sy.Matrix([[-1, 1], [0, 1]])).is_zero_matrix
for l in LINES:
    im = app(l["cls"], CJ)
    fam = [k for k in LINES if k["motion"] == l["motion"] and k["ell"] == l["ell"]]
    l["mirror"] = 1 if im == l["cls"] else (-1 if im == neg(l["cls"]) else 0)
    l["mirror_in"] = any(im in (k["cls"], neg(k["cls"])) for k in fam)
mir = {}
for l in LINES:
    mir.setdefault((l["motion"], l["ell"]), []).append(l["mirror"])
cat7 = [l for l in LINES if l["motion"] == "cat" and l["ell"] == 7]
s30 = [l for l in LINES if l["motion"] == "s10" and l["ell"] == 30]
h3_ok = (h3_gen and h3_ell and all(l["mirror_in"] for l in LINES)
         and all(l["mirror"] == 1 and not any(l["phi"]) and l["blocks"]["ell"][0] == 0 for l in cat7)
         and all(l["mirror"] == -1 and l["blocks"]["ell"][0] == -2 * l["blocks"]["ell"][1] for l in s30))
post("H3", "the mirror is an involution of Z^29 that reverses Psi and the three windings and acts on the 143a1 block by (n1, n2) -> (-n1, n1 + n2); it carries each line to a line of its own family, up to direction; the cat lines of length 7 are their own mirror images, direction kept, so their windings vanish and their 143a1 period is real; Sphereland 10's lines of length 30 are their own, direction reversed, so theirs is imaginary: a multiple of (2, -1)", h3_ok, "sign of the mirror on each class (0: carried to another line): " + "; ".join(f"{k[0]} {k[1]}: {v}" for k, v in mir.items()))

# H4: two windings and their shadows on X_0(11) and X_0(13)
h4_ok, shadow = True, {}
for l in LINES:
    A = MOTIONS[l["motion"]]; x = l["dart"]
    for q, key, i in ((11, "phiA", 1), (13, "phiB", 2)):
        mq = l["m"][q]
        sh = mm(mm(grep[x], mpow(A, mq)), inv(grep[x]))                 # the line's shadow: an element of Gamma_0(q)
        w = winding(sh, UNITS[key])
        h4_ok &= sh[1][0] % q == 0 and w.denominator == 1 and l["phi"][i] == (l["ell"] // mq) * w
        tab = shadow.setdefault(f"{l['motion']} mod {q}", {})
        okey = f"{mq}:{sorted(l['orb'][q])[0]}"
        h4_ok &= tab.setdefault(okey, int(w)) == int(w)
cat = [l for l in LINES if l["motion"] == "cat"]
cat_self13 = all(orbit_q(pts[pS[l["dart"]]], MOTIONS["cat"], 13) == l["orb"][13] for l in cat)
cat_swap11 = all(orbit_q(pts[pS[l["dart"]]], MOTIONS["cat"], 11) != l["orb"][11] for l in cat)
post("H4", "phi_A of a line is (ell / m_11) times the winding of its shadow on X_0(11), phi_B is (ell / m_13) times that of its shadow on X_0(13), and lines over one orbit have one shadow winding", h4_ok, "shadow windings by orbit (length:first point): " + "; ".join(f"{k}: {v}" for k, v in shadow.items()))
post("H4", "the cat map's partner through S keeps each orbit modulo 13 (-1 is a square) and exchanges the orbits modulo 11 (it is not): its shadows on X_0(13) are their own reverses, so phi_B = 0 on all eight cat lines", cat_self13 and cat_swap11 and all(l["phi"][2] == 0 for l in cat))

# H5: why the cat lines of length 7 have no f2 component
AL11 = [[int(v) for v in AL[11].row(i)] for i in range(B1)]
newpart = lambda l: tuple(l["blocks"]["ell"] + l["blocks"]["q4"] + l["blocks"]["q6"])
n0 = newpart(cat7[0])
h5_family = all(any(app(l["cls"], AL11) in (k["cls"], neg(k["cls"])) for k in cat7) for l in cat7)
h5_pm = all(newpart(l) in (n0, neg(n0)) for l in cat7)
post("H5", "W_11 carries the four cat lines of length 7 among themselves, up to direction; their new parts are +n or -n (H1 with P6); W_11 is +1 on 143a1 and f3 and -1 on f2, so n has no f2 component or else neither of the other two: it has no f2 component", h5_family and h5_pm and all(l["cond"][11] == 1 for l in cat7) and any(cat7[0]["blocks"]["ell"]) and any(cat7[0]["blocks"]["q6"]) and not any(cat7[0]["blocks"]["q4"]))

# where W_Q keeps a line in its family: the image as + or - the class of a line (index within the 28)
print("\n  images under the involutions that keep a line in its family (line -> sign and line)")
al_int = {Qd: [[int(v) for v in AL[Qd].row(i)] for i in range(B1)] for Qd in AL}
perm = []
for i, l in enumerate(LINES):
    for Qd in (11, 13, 143):
        if l["cond"][Qd] == 1:
            im = app(l["cls"], al_int[Qd])
            hits = [(s_, j) for j, k in enumerate(LINES) if k["motion"] == l["motion"] and k["ell"] == l["ell"] for s_ in (1, -1) if im == (k["cls"] if s_ == 1 else neg(k["cls"]))]
            perm.append([i, l["motion"], l["ell"], Qd, hits])
            print(f"  {i:2d} {l['motion']:5s} {l['ell']:3d}  W_{Qd}: {hits}")
print("\n  Sphereland 10's four lines of length 30 by their labels (Legendre symbol of the form modulo 11, modulo 13)")
for l in s30:
    print(f"   eps = ({l['eps'][11]:+d}, {l['eps'][13]:+d})   143a1 {l['blocks']['ell']}   windings {tuple(int(v) for v in l['phi'])}")
NUM["post"] = {"relations_H1": nrel, "ghost_nonzero_over_11": int(ghost_nz), "kept": {str(k): int(v) for k, v in kept.items()}, "mirror": {f"{k[0]} {k[1]}": v for k, v in mir.items()}, "shadow": shadow,
               "s10_30": [{"eps11": l["eps"][11], "eps13": l["eps"][13], "ell": l["blocks"]["ell"]} for l in s30], "al_images": perm}
for l, row in zip(LINES, NUM["lines"]):
    row["m11"], row["m13"], row["eps"], row["mirror"], row["conductor"] = l["m"][11], l["m"][13], [l["eps"][11], l["eps"][13]], l["mirror"], {str(k): v for k, v in l["cond"].items()}

npass = sum(1 for r in RESULTS if r[2])
ppass = sum(1 for r in POST if r[2])
print(f"\n{npass} of {len(RESULTS)} pre-registered checks passed; post hoc {ppass} of {len(POST)}; {time.time() - t_start:.0f} s")
for tag, name, okk, det in RESULTS + POST:
    if not okk:
        print("  MISS", tag, name, det)
if len(sys.argv) > 1:
    with open(sys.argv[1], "w") as fh:
        json.dump({"checks": [{"tag": t, "name": n, "pass": o, "detail": d} for t, n, o, d in RESULTS],
                   "post_hoc": [{"tag": t, "name": n, "pass": o, "detail": d} for t, n, o, d in POST], "numbers": NUM}, fh, indent=1)
sys.exit(0 if npass == len(RESULTS) and ppass == len(POST) else 1)
