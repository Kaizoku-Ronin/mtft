import sys, json, time, pickle
sys.path.insert(0, "/mnt/agents/output/mtft-v0140-work/src")
import numpy as np
from fractions import Fraction as Fr
from mtft import hecke as H
t0 = time.time()
V = pickle.load(open("V_exact.pkl","rb"))
def T_int(pr):
    T = H.cuspidal_hecke(pr)
    return [[Fr(T[i][j]) for j in range(26)] for i in range(26)]
Ts = {pr: T_int(pr) for pr in (2,3,5,11,13)}
def to_mod(A, p):
    B = np.zeros((26,26), dtype=np.int64)
    for i in range(26):
        for j in range(26):
            x = A[i][j]
            B[i,j] = (x.numerator % p) * pow(x.denominator % p, p-2, p) % p
    return B
def comm_map(A, p):
    Ap = np.mod(A, p).astype(np.int64)
    I = np.eye(26, dtype=np.int64)
    return (np.kron(Ap.T, I) - np.kron(I, Ap)) % p

def rref_nullity_nullspace(A, p, want_null=False):
    A = np.mod(A, p).astype(np.int64)
    nR, nC = A.shape
    piv_cols = []
    r = 0
    for c in range(nC):
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0: continue
        piv = r + nz[0]
        if piv != r: A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p-2, p)) % p
        rows = np.nonzero(A[:, c])[0]
        rows = rows[rows != r]
        if rows.size:
            fac = A[rows, c].copy()
            A[rows] = (A[rows] - np.outer(fac, A[r])) % p
        piv_cols.append(c); r += 1
        if r == nR: break
    nullity = nC - len(piv_cols)
    ns = []
    if want_null:
        free_cols = [c for c in range(nC) if c not in piv_cols]
        for fc in free_cols:
            v = np.zeros(nC, dtype=np.int64); v[fc] = 1
            for ri, pc in enumerate(piv_cols):
                v[pc] = (-A[ri, fc]) % p
            ns.append(v)
    return nullity, ns

p = 1000037
mats235 = [to_mod(Ts[pr], p) for pr in (2,3,5)]
matsAll = mats235 + [to_mod(Ts[11], p), to_mod(Ts[13], p)]
Vp = to_mod(V, p)
stack235V = np.vstack([comm_map(A, p) for A in mats235] + [comm_map(Vp, p)])
n235V, ns = rref_nullity_nullspace(stack235V, p, want_null=True)
print("{2,3,5}+V nullity:", n235V, f"[{time.time()-t0:.0f}s]", flush=True)
stackAllV = np.vstack([comm_map(A, p) for A in matsAll] + [comm_map(Vp, p)])
nAllV, _ = rref_nullity_nullspace(stackAllV, p)
print("{2,3,5,11,13}+V nullity:", nAllV, f"[{time.time()-t0:.0f}s]", flush=True)

# identify the second nullvector Z (first should be identity)
if len(ns) >= 2:
    for k, v in enumerate(ns):
        Z = v.reshape(26, 26)  # column-major vec
        print(f"nullvector {k}: diag entries mod p:", [int(Z[i,i]) for i in range(0,26,5)], flush=True)
    pickle.dump([v for v in ns], open("nullvecs_p1000037.pkl","wb"))
out = {"nullity_235_plus_V": n235V, "nullity_all_plus_V": nAllV}
json.dump(out, open("m8_deep.json","w"), indent=1)
print("DONE", flush=True)
