# C9B adjudication: the [inf, 2/77] sign via EXACT RATIONAL modular symbols.
# Route (no PARI anywhere): eigenfunctionals of the Hecke algebra on the
# 29-dim relative Manin-symbol quotient of X0(143), split by iota*,
# evaluated on CF decompositions of {inf,1/11} and {inf,2/77}.
# All matrix/data inputs validated in-script; final verdict assembled with
# independent period integrals (mpmath) and checked against agreed anchors.
import sys, json, time
sys.path.insert(0, "/mnt/agents/output/mtft-v0140-work/src")
from fractions import Fraction as Fr
import mpmath as mp
from mtft import hecke as H

t0 = time.time()
mp.mp.dps = 60
out = {"route": "exact rational Manin symbols, own code; module used for raw model data only"}

m = H.model()
N, E, nq = m["N"], m["E"], m["nq"]
P1, idx, canon = m["P1"], m["idx"], m["canon"]
eid, esign, erep = m["eid"], m["esign"], m["erep"]
free, Binv, cols, D2r = m["free"], m["Binv"], m["cols"], m["D2r"]
ncols = len(cols)

# ---------------- own exact linear algebra over Q ----------------
def matmul(A, B):
    n, k, l = len(A), len(B), len(B[0])
    C = [[Fr(0)] * l for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        for kk in range(k):
            a = Ai[kk]
            if a:
                Bk = B[kk]
                C[i] = [c + a * b for c, b in zip(C[i], Bk)]
    return C

def matvec(A, v):
    return [sum((a * b for a, b in zip(row, v)), Fr(0)) for row in A]

def rref_null(A):
    """RREF over Q; return (rref, pivots, nullspace basis as columns)."""
    A = [row[:] for row in A]
    nR, nC = len(A), len(A[0])
    piv, r = [], 0
    for c in range(nC):
        p = next((i for i in range(r, nR) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        f = A[r][c]
        A[r] = [x / f for x in A[r]]
        for i in range(nR):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(nC)]
        piv.append(c)
        r += 1
    freecols = [c for c in range(nC) if c not in set(piv)]
    basis = []
    for fc in freecols:
        v = [Fr(0)] * nC
        v[fc] = Fr(1)
        for i, pc in enumerate(piv):
            v[pc] = -A[i][fc]
        basis.append(v)
    return A, piv, basis

def transpose(A):
    return [list(r) for r in zip(*A)]

def eye(n):
    return [[Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]

# ---------------- validate Binv data: B . Binv == I ----------------
B = [[Fr(0)] * E for _ in range(E)]
for i in range(E):
    for j in range(ncols):
        B[i][j] = Fr(D2r[i][j])
    for j, e in enumerate(free):
        if i == e:
            B[i][ncols + j] = Fr(1)
prod = matmul(B, Binv)
out["Binv_is_inverse"] = all(prod[i][j] == (1 if i == j else 0) for i in range(E) for j in range(E))

# my own qcoords: x = Binv . v ; free coords are x[ncols:]
def qcoords(v):
    return matvec(Binv, v)[ncols:]

# ---------------- validate edge data: S-pairing consistency ----------------
sS = m["sS"]
nP1 = len(P1)
okS = all(eid[sS[i]] == eid[i] and esign[sS[i]] == -esign[i] for i in range(nP1))
out["edge_S_pairing_ok"] = bool(okS)

# ---------------- Hecke matrices on the 29-dim quotient ----------------
Tp = {p: [[Fr(x) for x in row] for row in H.hecke_matrix(p)] for p in (2, 3, 5, 7)}
# commutativity T_p T_q = T_q T_p
comm = all(matmul(Tp[p], Tp[q]) == matmul(Tp[q], Tp[p]) for p in (2, 3, 5, 7) for q in (2, 3, 5, 7))
out["hecke_commute_2357"] = bool(comm)

# ---------------- iota* on the 29-dim quotient ----------------
iota = m["iota"]
IS = [[0] * E for _ in range(E)]
for k in range(E):
    y = iota[erep[k]]
    IS[eid[y]][k] += esign[y]
iota29 = [[Fr(x) for x in row] for row in H._quotient(m, IS)]
I29 = eye(nq)
out["iota29_involution"] = matmul(iota29, iota29) == I29
out["iota29_hecke_commute"] = all(matmul(iota29, Tp[p]) == matmul(Tp[p], iota29) for p in (2, 3, 5, 7))

# ---------------- left common eigenspaces ----------------
def left_common_eig(ap):
    """ap: dict p -> eigenvalue.  Stack (T_p^T - a_p I); return nullspace."""
    rows = []
    for p, a in ap.items():
        Tt = transpose(Tp[p])
        for i in range(nq):
            row = Tt[i][:]
            row[i] -= Fr(a)
            rows.append(row)
    _, _, basis = rref_null(rows)
    return basis  # list of columns (29-vectors)

f1_ap = {2: 0, 3: -1, 5: -1, 7: -2}        # 143a1 eigenform
gh_ap = {2: -2, 3: -1, 5: 1, 7: -2}        # 11a1 ghost (11.a class)
L_f1 = left_common_eig(f1_ap)
L_gh = left_common_eig(gh_ap)
out["dim_L_f1"] = len(L_f1)
out["dim_L_ghost"] = len(L_gh)

# ---------------- split L by iota* (columns: iota action l -> iota29^T l) ----------------
def iota_eigspaces(L):
    """L: list of 29-vectors spanning an iota-invariant subspace.
    Returns (plus_basis, minus_basis) of the +/-1 eigespaces, plus the op matrix."""
    def dot(u, v):
        return sum((x * y for x, y in zip(u, v)), Fr(0))
    nL = len(L)
    il = [matvec(transpose(iota29), l) for l in L]
    A1 = [[dot(L[i], L[j]) for j in range(nL)] for i in range(nL)]
    A2 = [[dot(L[i], il[j]) for j in range(nL)] for i in range(nL)]
    # solve A1 . X = A2 over Q (A1 invertible: Gram of independent vectors)
    A1aug = [A1[i][:] + A2[i][:] for i in range(nL)]
    R, piv, _ = rref_null(A1aug)
    A = [[R[i][nL + j] for j in range(nL)] for i in range(nL)]
    def eig(sig):
        M = [[A[i][j] - (Fr(sig) if i == j else Fr(0)) for j in range(nL)] for i in range(nL)]
        _, _, ns = rref_null(M)
        return [[sum((L[k][i] * x[k] for k in range(nL)), Fr(0)) for i in range(nq)] for x in ns]
    return eig(1), eig(-1), A

plus_f1, minus_f1, A_f1 = iota_eigspaces(L_f1)
plus_gh, minus_gh, A_gh = iota_eigspaces(L_gh)
lp_f1, lm_f1 = plus_f1[0], minus_f1[0]
out["iota_split_f1_dims"] = [len(plus_f1), len(minus_f1)]
out["iota_split_ghost_dims"] = [len(plus_gh), len(minus_gh)]
out["iota_on_Lf1"] = [[str(x) for x in r] for r in A_f1]
out["iota_on_Lf1_sq_is_I"] = matmul(A_f1, A_f1) == eye(2)
# eigen check: lp/lm are iota eigenvectors and Hecke eigenfunctionals
def is_left_eig(lam, ap):
    for p, a in ap.items():
        r = matvec(transpose(Tp[p]), lam)
        if any(r[i] != Fr(a) * lam[i] for i in range(nq)):
            return False
    return True
# preservation check: iota maps span(L) into itself exactly
def preserves(L, A):
    il = [matvec(transpose(iota29), l) for l in L]
    nL = len(L)
    for j in range(nL):
        rec = [sum((L[k][i] * A[k][j] for k in range(nL)), Fr(0)) for i in range(nq)]
        if rec != il[j]:
            return False
    return True
out["iota_preserves_Lf1"] = preserves(L_f1, A_f1)
out["iota_preserves_Lgh"] = preserves(L_gh, A_gh)
out["iota_on_Lgh_sq_is_I"] = matmul(A_gh, A_gh) == eye(len(L_gh))
out["lp_f1_eig_ok"] = is_left_eig(lp_f1, f1_ap) and matvec(transpose(iota29), lp_f1) == lp_f1
out["lm_f1_eig_ok"] = is_left_eig(lm_f1, f1_ap) and all(x == -y for x, y in zip(matvec(transpose(iota29), lm_f1), lm_f1))
out["plus_gh_eig_ok"] = all(is_left_eig(l, gh_ap) for l in plus_gh)
out["minus_gh_eig_ok"] = all(is_left_eig(l, gh_ap) for l in minus_gh)

# primitive integer scaling for reporting
def prim(v):
    from math import gcd
    den = 1
    for x in v:
        den = den * x.denominator // gcd(den, x.denominator)
    w = [int(x * den) for x in v]
    g = 0
    for x in w:
        g = gcd(g, abs(x))
    sgn = 1
    for x in w:
        if x:
            sgn = 1 if x > 0 else -1
            break
    return [sgn * x // g for x in w]
out["lambda_plus_f1_int"] = prim(lp_f1)
out["lambda_minus_f1_int"] = prim(lm_f1)

# ---------------- CF decomposition of {inf, r} ----------------
def convergents(p, q):
    convs = []
    pm2, pm1, qm2, qm1 = 0, 1, 1, 0
    a, b = p, q
    while True:
        ak = a // b
        pk, qk = ak * pm1 + pm2, ak * qm1 + qm2
        convs.append((pk, qk))
        if (pk, qk) == (p, q):
            break
        a, b = b, a - ak * b
        pm2, pm1, qm2, qm1 = pm1, pk, qm1, qk
    return convs

def msym_edge(c, d):
    i = idx[canon(c % N, d % N)]
    return eid[i], esign[i]

def edge_vec_symbol_0_r(p, q):
    """the relative symbol {0, p/q} as an edge-space vector (len E)."""
    v = [Fr(0)] * E
    convs = convergents(p, q)
    for k in range(1, len(convs)):
        qk = convs[k][1]
        qkm1 = convs[k - 1][1]
        d = qkm1 if k % 2 == 1 else -qkm1   # (-1)^{k+1} q_{k-1}
        e, s = msym_edge(qk, d)
        v[e] += Fr(s)
    return v

def edge_vec_symbol_0_inf():
    v = [Fr(0)] * E
    e, s = msym_edge(0, 1)
    v[e] += Fr(s)
    return v

vinf = edge_vec_symbol_0_inf()
def b_symbol(p, q):
    v = edge_vec_symbol_0_r(p, q)
    return qcoords([v[i] - vinf[i] for i in range(E)])

b11 = b_symbol(1, 11)
b277 = b_symbol(2, 77)
binf_q = qcoords(vinf)   # {0, inf} itself, for the L(f,1) control

def pair(lam, b):
    return sum((x * y for x, y in zip(lam, b)), Fr(0))

# ---------------- ghost canonical line: pin by factoring through level 11 ----------------
# pi: oriented symbol (c:d) at 143 -> oriented class at 11 (commutes with canon and iota)
import math as _math
units11 = [u for u in range(1, 11) if _math.gcd(u, 11) == 1]
def canon11(c, d):
    return min(((c * u) % 11, (d * u) % 11) for u in units11)
nP1 = len(P1)
pi_class = [canon11(P1[i][0] % 11, P1[i][1] % 11) for i in range(nP1)]
def F_vec(i):
    v = [Fr(0)] * E
    v[eid[i]] = Fr(esign[i])
    return qcoords(v)
first_of, eq_rows = {}, []
for i in range(nP1):
    cl = pi_class[i]
    if cl in first_of:
        Fi, Fj = F_vec(i), F_vec(first_of[cl])
        eq_rows.append([Fi[k] - Fj[k] for k in range(nq)])
    else:
        first_of[cl] = i
# pinned subspace of L_gh: lambda in span(L_gh) with lambda . (F_i - F_j) = 0
nG = len(L_gh)
EqL = [[pair(row, L_gh[k]) for k in range(nG)] for row in eq_rows]
_, _, pinned_coords = rref_null(EqL)
pinned = [[sum((L_gh[k][i] * x[k] for k in range(nG)), Fr(0)) for i in range(nq)]
          for x in pinned_coords]
out["ghost_pinned_dim"] = len(pinned)
plus_g, minus_g, A_g = iota_eigspaces(pinned)
out["ghost_pinned_iota_split"] = [len(plus_g), len(minus_g)]
lam_ghp, lam_ghm = plus_g[0], minus_g[0]
# canonical (Manin-primitive) normalization of the REAL-direction functional
# lam_ghm: values on the 26 cuspidal K-cycles scaled to coprime integers
Kcyc = m["K"]
raw_vals = [pair(lam_ghm, [Fr(x) for x in Kcyc[i]]) for i in range(len(Kcyc))]
nz = [v for v in raw_vals if v != 0]
gg = nz[0]
for v in nz[1:]:
    from fractions import Fraction as _Fr
    gg = _Fr(abs(_math.gcd(gg.numerator * v.denominator, v.numerator * gg.denominator)),
             gg.denominator * v.denominator)
lam_ghm = [x / gg for x in lam_ghm]
out["ghost_canon_plus_int"] = prim(lam_ghp)
out["ghost_canon_minus_int"] = prim(lam_ghm)

# CONVENTION (verified in-script): the iota-EVEN functional lam+ carries the
# IMAGINARY part of the period map; the iota-ODD functional lam- carries the
# REAL part.  (phi o iota = -conj(phi) for real-coefficient f.)
vals = {
    "lp_11": pair(lp_f1, b11), "lm_11": pair(lm_f1, b11),   # lp=imag dir, lm=real dir
    "lp_277": pair(lp_f1, b277), "lm_277": pair(lm_f1, b277),
    "lp_inf_f1": pair(lp_f1, binf_q), "lm_inf_f1": pair(lm_f1, binf_q),
    "lp_inf_gh_canon": pair(lam_ghp, binf_q), "lm_inf_gh_canon": pair(lam_ghm, binf_q),
    "lp_11_gh_canon": pair(lam_ghp, b11), "lm_11_gh_canon": pair(lam_ghm, b11),
}
out["symbol_values"] = {k: str(v) for k, v in vals.items()}
out["check_per11_real_part_zero"] = vals["lm_11"] == 0     # per11 pure imaginary
out["check_per11_imag_part_nonzero"] = vals["lp_11"] != 0
out["check_f1_rank1_L1_zero"] = vals["lp_inf_f1"] == 0 and vals["lm_inf_f1"] == 0
# ghost BSD anchor: s({0,inf}) = -L(11a1,1) is REAL (rank 0): imag dir vanishes,
# real dir (canonically normalized on cuspidal cycles) is a small rational
out["check_ghost_rank0_imag_vanishes"] = vals["lp_inf_gh_canon"] == 0
out["ghost_BSD_symbol_value_realdir"] = str(vals["lm_inf_gh_canon"])
out["check_ghost_BSD_small_rational"] = vals["lm_inf_gh_canon"] in (Fr(-1, 5), Fr(-2, 5), Fr(1, 5), Fr(2, 5))

# decisive exact rationals (with the corrected convention):
#   per(r) = Omega_re * lm(r) + i * Omega_im * lp(r)
#   lambda1 = per277/per11 ;  per11 = i*Omega_im*lp_11 (lm_11 = 0)
#   => Re(lambda1) = lp_277/lp_11  (EXACT, orientation-free)
#      |Im(lambda1)| = |lm_277/lp_11| * Omega_re/Omega_im  (magnitude check)
Re_lam1 = vals["lp_277"] / vals["lp_11"]
Im_coeff = vals["lm_277"] / vals["lp_11"]
out["Re_lambda1_exact"] = str(Re_lam1)
out["Im_lambda1_coeff_exact"] = str(Im_coeff)   # = lm_277/lp_11

# ---------------- independent periods via contour integrals ----------------
def periods_from_g(g, guess, label):
    e = mp.findroot(g, guess)
    O1 = 2 * mp.quad(lambda x: 1 / mp.sqrt(g(x)), [e, mp.inf])
    Om = 2 * mp.quad(lambda x: 1 / mp.sqrt(-g(x)), [-mp.inf, e])
    out[label] = {"real_root": mp.nstr(e, 50), "Omega_real": mp.nstr(O1, 50), "Omega_imag_abs": mp.nstr(Om, 50)}
    return e, O1, Om

g1 = lambda x: 4 * x**3 - 4 * x**2 - 4 * x - 7          # 143a1: (2y+1)^2
e1, O1p, O1m = periods_from_g(g1, (mp.mpf("1.5"), mp.mpf("2.0")), "periods_143a1")
# ghost 11a1: y^2+y = x^3 - x^2 - 10x - 20 -> (2y+1)^2 = 4x^3-4x^2-40x-79
g2 = lambda x: 4 * x**3 - 4 * x**2 - 40 * x - 79
e2, O2p, O2m = periods_from_g(g2, (mp.mpf("4.0"), mp.mpf("5.0")), "periods_11a1")

# ---------------- ghost full-stack anchor: lambda+_gh(b_inf) . Omega+ == -L(11a1,1) ----------------
# L(11a1,1) from my own sieve a_p (ap_sieve.json): epsilon=+1, L(1) = 2 sum a_n/n e^{-2 pi n/sqrt(11)}
sieve = json.load(open("/mnt/agents/output/scratch_v0150/ap_sieve.json"))
ap11 = {int(k): v for k, v in sieve["11a1"].items()}
primes_sorted = sorted(ap11)
M = 6000
an = [0] * (M + 1)
an[1] = 1
for n in range(2, M + 1):
    f = 0
    for p in primes_sorted:
        if p * p > n:
            break
        if n % p == 0:
            f = p
            break
    if f == 0:
        an[n] = ap11[n]          # n is prime and <= 6000 << 60000 table
    else:
        e_, mm = 0, n
        while mm % f == 0:
            mm //= f
            e_ += 1
        seq = [1, ap11[f]]
        for k in range(2, e_ + 1):
            if f == 11:
                seq.append(ap11[f] * seq[-1])          # bad prime (11 | level): a_{p^k} = a_p^k
            else:
                seq.append(ap11[f] * seq[-1] - f * seq[-2])
        an[n] = seq[e_] * an[mm]
s = mp.mpf(0)
for n in range(1, M + 1):
    s += mp.mpf(an[n]) / n * mp.exp(-2 * mp.pi * n / mp.sqrt(11))
L11 = 2 * s
out["L_11a1_1_own"] = mp.nstr(L11, 40)
out["ghost_L_over_Omega_own"] = mp.nstr(L11 / O2p, 30)   # expect 0.2 (BSD 1/5)
rhs = mp.mpf(vals["lm_inf_gh_canon"].numerator) / vals["lm_inf_gh_canon"].denominator
out["ghost_canon_symbol_vs_L_over_Omega"] = mp.nstr(rhs / (L11 / O2p), 20)  # expect -1 or -2

# ---------------- assemble per11, per277, lambda1 from MY data ----------------
# per(r) = Omega_re*lm(r) + i*Omega_im*lp(r); per11 = i*Omega_im*lp_11.
# Re(lambda1) = lp_277/lp_11 EXACT.  |Im| uses Omega_re/Omega_im from my own
# integrals (note: for a one-component curve the staggered lattice can put a
# factor 1/2 in the effective real period -- flagged if off by ~2x).
lp11, lm11 = vals["lp_11"], vals["lm_11"]
lp277, lm277 = vals["lp_277"], vals["lm_277"]
f_lp11 = mp.mpf(lp11.numerator) / lp11.denominator
f_lp277 = mp.mpf(lp277.numerator) / lp277.denominator
f_lm277 = mp.mpf(lm277.numerator) / lm277.denominator
sgn = -1 if f_lp11 > 0 else 1      # orient i*Omega_im so per11 matches archive sign
Oim = sgn * O1m
per11_mine = 1j * Oim * f_lp11
per11_archive = -mp.mpf("0.31352300915287423503215575582714111095417953831197") * 1j
out["per11_mine"] = mp.nstr(per11_mine, 50)
out["per11_abs_err_vs_archive"] = mp.nstr(abs(per11_mine - per11_archive), 8)
out["Omega_im_over_abs_per11"] = mp.nstr(O1m / mp.mpf("0.31352300915287423503215575582714111095417953831197"), 30)
out["lp_11_exact"] = str(lp11)
lam1_Re = mp.mpf(lp277.numerator) / lp277.denominator
lam1_Im_mag = abs(mp.mpf(lm277.numerator) / lm277.denominator) * (O1p / O1m)
out["lambda1_Re_exact"] = str(Fr(lp277) / Fr(lp11))
out["abs_Im_lambda1_mine"] = mp.nstr(lam1_Im_mag, 50)
out["abs_Im_lambda1_C6_anchor"] = "1.0232745926964612"
out["abs_Im_err"] = mp.nstr(abs(lam1_Im_mag - mp.mpf("1.0232745926964612")), 8)
out["abs_Im_err_with_halved_real_period"] = mp.nstr(abs(lam1_Im_mag / 2 - mp.mpf("1.0232745926964612")), 8)
Re_exact = Fr(lp277) / Fr(lp11)
if Re_exact == Fr(1, 2):
    verdict = "ARCHIVE (PARI 2.17) sign CORRECT: Re lambda1 = +1/2 exactly; v5/v6 (PARI 2.15.4) sign is the artifact"
elif Re_exact == Fr(-1, 2):
    verdict = "V6 (PARI 2.15.4) sign CORRECT: Re lambda1 = -1/2 exactly; Paper 33 v2 archive (PARI 2.17) needs the sign correction"
else:
    verdict = f"BOTH WRONG: Re lambda1 = {Re_exact} (not +/-1/2) -- needs deeper investigation"
out["VERDICT"] = verdict
out["elapsed_s"] = round(time.time() - t0, 1)

with open("/mnt/agents/output/scratch_v0150/c9b_exact_symbol.json", "w") as fh:
    json.dump(out, fh, indent=2, default=str)
for k, v in out.items():
    print(f"{k}: {v}")
print("DONE", out["elapsed_s"], "s")
