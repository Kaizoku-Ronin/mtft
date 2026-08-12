# M8/M8b/M9-R1R2 independent verification -- my own code throughout.
import sys, json, time, math
sys.path.insert(0, "/mnt/agents/output/mtft-v0140-work/src")
import numpy as np
from fractions import Fraction as Fr
from math import gcd
from mtft import hecke as H

t0 = time.time()
out = {}
m = H.model()
E, D2, Kk, free, nq = m["E"], m["D2"], m["K"], m["free"], m["nq"]
tris, tri_of, sS = m["tris"], m["tri_of"], m["sS"]
erep, cusp_of, fans = m["erep"], m["cusp_of"], m["fans"]

# ---------- my own exact RREF over Q (no use of H._rref) ----------
def rref_q(A):
    A = [row[:] for row in A]
    nR, nC = len(A), len(A[0])
    piv, r = [], 0
    for c in range(nC):
        p = next((i for i in range(r, nR) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        f = A[r][c]
        A[r] = [x / f for x in A[r]]
        for i in range(nR):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(nC)]
        piv.append(c); r += 1
        if r == nR: break
    return A, piv

# ---------- harmonic embedding, my route ----------
G56 = [[sum(D2[e][i]*D2[e][j] for e in range(E)) for j in range(56)]
       for i in range(56)]
def harmonic(v26):
    v = [sum(Kk[a][j]*v26[a] for a in range(26)) for j in range(nq)]
    v84 = [Fr(0)]*E
    for j, e in enumerate(free): v84[e] = v[j]
    rhs = [sum(Fr(D2[e][i])*v84[e] for e in range(E)) for i in range(56)]
    Aug = [[Fr(G56[i][j]) for j in range(56)] + [rhs[i]] for i in range(56)]
    R, pv = rref_q(Aug)
    alpha = [Fr(0)]*56
    for r_, c_ in enumerate(pv):
        if c_ < 56: alpha[c_] = R[r_][56]
    return [v84[e] - sum(Fr(D2[e][i])*alpha[i] for i in range(56))
            for e in range(E)]

basis = []
for a in range(26):
    v = [Fr(0)]*26; v[a] = Fr(1)
    basis.append(harmonic(v))
print("harmonic basis built", time.time()-t0, flush=True)
Gm = [[sum(basis[a][e]*basis[b][e] for e in range(E)) for b in range(26)]
      for a in range(26)]
# invert G (26x26 exact)
Aug = [[Gm[i][j] for j in range(26)] + [Fr(1) if i == j else Fr(0)
       for j in range(26)] for i in range(26)]
R, _ = rref_q(Aug)
Ginv = [row[26:] for row in R]

width = {k: len(o) for k, o in enumerate(fans)}
g_w = [sum(width[cusp_of[f]] for f in tris[t]) for t in range(56)]
def coupling(g):
    gavg = [Fr(g[tri_of[erep[k]]] + g[tri_of[sS[erep[k]]]], 2)
            for k in range(E)]
    M = [[sum(basis[a][e]*gavg[e]*basis[b][e] for e in range(E))
          for b in range(26)] for a in range(26)]
    return [[sum(Ginv[i][k]*M[k][j] for k in range(26))
             for j in range(26)] for i in range(26)]
V = coupling(g_w)
print("V built", time.time()-t0, flush=True)

# ---------- anchors vs M7 ledger ----------
blocks = H.blocks()
order = ["ell", "old", "q4", "q6"]
idx = {nm: list(range(sum(len(blocks[o]) for o in order[:k]),
                      sum(len(blocks[o]) for o in order[:k+1])))
       for k, nm in enumerate(order)}
cols = [list(v) for nm in order for v in blocks[nm]]
S = [[cols[b][i] for b in range(26)] for i in range(26)]
RS, _ = rref_q([[Fr(x) for x in row] for row in S])
Sinv = [row[:] for row in RS]  # rref_q reduced [S]? need inverse
# recompute inverse properly
AugS = [[Fr(S[i][j]) for j in range(26)] + [Fr(1) if i == j else Fr(0)
        for j in range(26)] for i in range(26)]
RS, _ = rref_q(AugS)
Sinv = [row[26:] for row in RS]
def mm(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]
Vb = mm(Sinv, mm(V, S))
off2 = Fr(0)
quiet = {}
for nm in order:
    tot = Fr(0)
    for i in idx[nm]:
        for j in idx[nm]:
            tot += Vb[i][j]**2
    quiet[nm] = tot
for a, na in enumerate(order):
    for b, nb in enumerate(order):
        if a != b:
            for i in idx[na]:
                for j in idx[nb]:
                    off2 += Vb[i][j]**2
import math as _m
out["anchor_offblock_norm"] = _m.sqrt(float(off2))
out["anchor_quietness"] = {k: _m.sqrt(float(v)) for k, v in quiet.items()}
print("anchors:", out["anchor_offblock_norm"], out["anchor_quietness"],
      time.time()-t0, flush=True)

# ---------- my own GF(p) rank ----------
def rank_mod_p(A, p):
    A = np.mod(A, p).astype(np.int64)
    nR, nC = A.shape
    r = 0
    for c in range(nC):
        piv = -1
        col = A[r:, c]
        nz = np.nonzero(col)[0]
        if nz.size == 0: continue
        piv = r + nz[0]
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p-2, p)) % p
        rows = np.arange(r+1, nR)
        fac = A[rows, c]
        nzm = fac != 0
        rr = rows[nzm]
        if rr.size:
            A[rr] = (A[rr] - np.outer(fac[nzm], A[r])) % p
        r += 1
        if r == nR: break
    return r

def comm_map(A, p):
    # X -> [X, A] = XA - AX ; column-major vec: (A^T (x) I - I (x) A)
    Ap = np.mod(A, p).astype(np.int64)
    I = np.eye(26, dtype=np.int64)
    return (np.kron(Ap.T, I) - np.kron(I, Ap)) % p

def T_int(pr):
    T = H.cuspidal_hecke(pr)
    return [[Fr(T[i][j]) for j in range(26)] for i in range(26)]
Ts = {pr: T_int(pr) for pr in (2, 3, 5, 11, 13)}

def to_mod(A, p):
    B = np.zeros((26, 26), dtype=np.int64)
    for i in range(26):
        for j in range(26):
            x = A[i][j]
            B[i, j] = (x.numerator % p) * pow(x.denominator % p, p-2, p) % p
    return B

PRIMES = (1000003, 999983, 1000033, 1000037, 999999937, 1000000007)
dimsH, dimsHV, dimsV, dimsFull = {}, {}, {}, {}
for p in PRIMES:
    mats = [to_mod(Ts[pr], p) for pr in (2, 3, 5)]
    stack = np.vstack([comm_map(A, p) for A in mats])
    dimsH[p] = 676 - rank_mod_p(stack, p)
    print("p", p, "hecke commutant", dimsH[p], time.time()-t0, flush=True)
    Vp = to_mod(V, p)
    stackV = np.vstack([stack, comm_map(Vp, p)])
    dimsHV[p] = 676 - rank_mod_p(stackV, p)
    print("p", p, "hecke+V commutant", dimsHV[p], time.time()-t0, flush=True)
    dimsV[p] = 676 - rank_mod_p(comm_map(Vp, p), p)
    print("p", p, "V-alone commutant", dimsV[p], time.time()-t0, flush=True)
    matsF = mats + [to_mod(Ts[11], p), to_mod(Ts[13], p)]
    stackF = np.vstack([comm_map(A, p) for A in matsF])
    dimsFull[p] = 676 - rank_mod_p(stackF, p)
    print("p", p, "full-hecke commutant", dimsFull[p], time.time()-t0, flush=True)
out["commutant_hecke"] = dimsH
out["commutant_hecke_plus_V"] = dimsHV
out["commutant_V_alone"] = dimsV
out["commutant_full_hecke"] = dimsFull

# ---------- M8b: JJ and exact amplitude ----------
U13 = Ts[13]
# ghost block = 'old' (4-dim); U13 restricted: pull back block basis
U13b = mm(Sinv, mm(U13, S))
i0, i1 = idx["old"][0], idx["old"][-1]+1
Ug = [[U13b[i][j] for j in range(i0, i1)] for i in range(i0, i1)]
tr = sum(Ug[i][i] for i in range(4))
disc = tr*tr - 4*13  # per 2-dim factor: tr/2? study uses tr over d//2
tr2 = Fr(sum(Ug[i][i] for i in range(4)), 2)
disc2 = tr2*tr2 - 4*13
JJ = [[(Ug[i][j] - (2 if i == j else 0))/Fr(3) for j in range(4)]
      for i in range(4)]
J2 = mm(JJ, JJ)
sq_ok = all(J2[i][j] == (-1 if i == j else 0) for i in range(4) for j in range(4))
# M9-R1: (2U - a)^2 = (a^2 - 4q) I with a=4, q=13
R1 = mm([[2*Ug[i][j] - (4 if i == j else 0) for j in range(4)] for i in range(4)],
        [[2*Ug[i][j] - (4 if i == j else 0) for j in range(4)] for i in range(4)])
r1_ok = all(R1[i][j] == (-36 if i == j else 0) for i in range(4) for j in range(4))
Vold = [[Vb[i][j] for j in range(i0, i1)] for i in range(i0, i1)]
JVJ = mm(JJ, mm(Vold, JJ))
Vminus = [[(Vold[i][j] + JVJ[i][j])/2 for j in range(4)] for i in range(4)]
Vplus  = [[(Vold[i][j] - JVJ[i][j])/2 for j in range(4)] for i in range(4)]
fminus = sum(Vminus[i][j]**2 for i in range(4) for j in range(4))
fplus  = sum(Vplus[i][j]**2 for i in range(4) for j in range(4))
ftot   = sum(Vold[i][j]**2 for i in range(4) for j in range(4))
comm = mm(Vold, JJ)
commJ = mm(JJ, Vold)
nz = sum(1 for i in range(4) for j in range(4) if comm[i][j] + commJ[i][j] != 0)
# Hecke null control: T_p old-block antilinear part
ctrl = {}
for pr in (2, 3, 5):
    Tb = mm(Sinv, mm(Ts[pr], S))
    Told = [[Tb[i][j] for j in range(i0, i1)] for i in range(i0, i1)]
    TJT = mm(JJ, mm(Told, JJ))
    Tm = [[(Told[i][j] + TJT[i][j])/2 for j in range(4)] for i in range(4)]
    ctrl[pr] = sum(abs(Tm[i][j]) for i in range(4) for j in range(4)) == 0
out["m8b"] = dict(trace_per_factor=str(tr2), discriminant=str(disc2),
                  JJ_squared_minus_I=sq_ok, R1_minus36I=r1_ok,
                  nonzero_anticommutator_entries=nz,
                  fminus=str(fminus), ftot=str(ftot),
                  antilinear_fraction=str(Fr(fminus, ftot)),
                  antilinear_float=float(Fr(fminus, ftot)),
                  hecke_null_control=ctrl)
print("M8b:", out["m8b"], time.time()-t0, flush=True)

out["runtime_s"] = round(time.time()-t0, 1)
json.dump(out, open("m8_verify.json", "w"), indent=1, default=str)
print("DONE", out["runtime_s"], flush=True)
