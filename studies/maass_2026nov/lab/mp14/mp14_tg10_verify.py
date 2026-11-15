#!/usr/bin/env python3
"""
mp14_tg10_verify.py  --  Route A for MP14 (test T-G10 of the gravity map, piece II.10,
and the gap of X_0(143) at infinite horodepth)

Part A.  Why does every Pell line have signed crossing count zero with the 24 geodesics of
MP13?  The closed geodesics are cut along F in exact rational arithmetic as in MP13, now
for six motions (the cat map, Sphereland 10's, the Pell unit, and the new UV^3, U^2V^2,
UV^5) and at the levels 143, 11 and 13; crossings are found sheet by sheet; the mirror
tau -> -conj(tau) is followed through the darts.

Part B.  The Eisenstein series of Gamma_0(143) from the level-1 series, by mpmath: the
scattering matrix, its value -I at s = 1/2, and the series along the Pell line, which the
reflection sigma(tau) = 1/(143 conj(tau)) fixes.

Predictions A1 to A4 and B1 to B6 are in MP14_TG10_Plan.md and were written before the run.
Route B: mp14_lattice.py sums the series over the group with no use of level 1.

Run:  python3 mp14_tg10_verify.py ../mp12/mp12_numbers.json ../mp13/mp13_numbers.json mp14_numbers.json
Needs mpmath and numpy.
"""
import sys, json, math, time
from fractions import Fraction as Fr
from itertools import combinations
import mpmath as mp
import numpy as np

t_start = time.time()
RESULTS, NUM = [], {}
R12 = json.load(open(sys.argv[1]))["numbers"]
R13 = json.load(open(sys.argv[2]))["numbers"]


def check(tag, name, cond, detail=""):
    RESULTS.append((tag, name, bool(cond), detail))
    print(("ok   " if cond else "FAIL ") + f"{tag:4s} {name}" + (f"   [{detail}]" if detail else ""), flush=True)


# ---------------------------------------------------------------- matrices and exact arcs (as in MP13)
def mm(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))


inv = lambda A: ((A[1][1], -A[0][1]), (-A[1][0], A[0][0]))
neg = lambda A: ((-A[0][0], -A[0][1]), (-A[1][0], -A[1][1]))
pos = lambda A: A if A[0][0] + A[1][1] > 0 else neg(A)
same = lambda A, B: A == B or A == neg(B)
conj = lambda g, A: pos(mm(mm(inv(g), A), g))
tpow = lambda n: ((1, n), (0, 1))
I2, S, T, Ti = ((1, 0), (0, 1)), ((0, -1), (1, 0)), ((1, 1), (0, 1)), ((1, -1), (0, 1))
R = ((0, 1), (-1, -1))
R2 = mm(R, R)
HALF, Q34 = Fr(1, 2), Fr(3, 4)
form = lambda A: (A[1][0], A[1][1] - A[0][0], -A[0][1])


def arc(A):
    a, b, c = form(A)
    pts = set()
    yl = Fr(b, 2 * a) - Fr(c, a) - Fr(1, 4)
    if yl >= Q34:
        pts.add((-HALF, yl))
    yr = Fr(-b, 2 * a) - Fr(c, a) - Fr(1, 4)
    if yr >= Q34:
        pts.add((HALF, yr))
    if b != 0:
        xb = Fr(-(a + c), b)
        if abs(xb) <= HALF:
            pts.add((xb, 1 - xb * xb))
    pts = sorted(pts)
    if len(pts) != 2:
        return None
    s = 1 if a > 0 else -1
    en, ex = (pts[0], pts[1]) if s > 0 else (pts[1], pts[0])
    return {"N": A, "form": (a, b, c), "s": s, "center": Fr(-b, 2 * a), "entry": en, "exit": ex}


def kind(p):
    x, y2 = p
    on = x * x + y2 == 1
    if x == -HALF and on:
        return "rho"
    if x == HALF and on:
        return "rho1"
    if x == -HALF:
        return "left"
    if x == HALF:
        return "right"
    if on:
        return "i" if x == 0 else ("bot-" if x < 0 else "bot+")
    return "int"


RHO, RHO1 = complex(-0.5, math.sqrt(3) / 2), complex(0.5, math.sqrt(3) / 2)
CORNER = {"rho": RHO, "rho1": RHO1}
CAND = {"rho1": [T, mm(T, R), mm(T, R2), mm(mm(T, R), Ti), mm(mm(T, R2), Ti)], "rho": [R, R2, Ti, mm(R, Ti), mm(R2, Ti)]}
mob = lambda g, z: (g[0][0] * z + g[0][1]) / (g[1][0] * z + g[1][1])


def step(ar):
    k = kind(ar["exit"])
    x, y2 = ar["exit"]
    if k in ("right", "left", "i", "bot-", "bot+"):
        g = {"right": T, "left": Ti}.get(k, S)
        nxt = arc(conj(g, ar["N"]))
        assert nxt is not None and nxt["entry"] == (-x, y2)
        return g, nxt
    good = []
    for g in CAND[k]:
        nxt = arc(conj(g, ar["N"]))
        if nxt is None or kind(nxt["entry"]) not in ("rho", "rho1"):
            continue
        if abs(mob(g, CORNER[kind(nxt["entry"])]) - CORNER[k]) < 1e-9:
            good.append((g, nxt))
    assert len(good) == 1
    return good[0]


def reduce_point(z):
    m = I2
    for _ in range(100000):
        n = int(mp.nint(z.real))
        if n:
            z, m = z - n, mm(m, tpow(n))
        if abs(z) < 1:
            z, m = -1 / z, mm(m, S)
        else:
            return m, z
    raise RuntimeError


def walk(A):
    A = pos(A)
    a, b, c = form(A)
    mp.mp.dps = 60
    xc, r = mp.mpf(-b) / (2 * a), mp.sqrt(b * b - 4 * a * c) / (2 * abs(a))
    m0, _ = reduce_point(xc + r * mp.expjpi(mp.mpf("0.3141")))
    ar = arc(conj(m0, A))
    arcs, m = [], m0
    while True:
        ar["m"] = m
        arcs.append(ar)
        g, ar = step(ar)
        m = mm(m, g)
        if ar["N"] == arcs[0]["N"]:
            break
    assert same(mm(m, inv(m0)), A)
    return arcs


def cross(fa, fb):
    a1, b1, c1 = fa
    a2, b2, c2 = fb
    den = a2 * b1 - a1 * b2
    if den == 0:
        return None
    x = Fr(a1 * c2 - a2 * c1, den)
    y2 = -(a1 * x * x + b1 * x + c1) / Fr(a1)
    if y2 <= 0 or abs(x) > HALF or x * x + y2 < 1:
        return None
    return (x, y2)


sgn = lambda v: (v > 0) - (v < 0)
unor = lambda f: f if f[0] > 0 else tuple(-v for v in f)
cj = lambda g: ((g[0][0], -g[0][1]), (-g[1][0], g[1][1]))

MOTIONS = {"cat": ((2, 1), (1, 1)), "s10": ((84, -47), (143, -80)), "pell": ((12, 1), (143, 12)),
           "uv3": ((4, 1), (3, 1)), "u2v2": ((5, 2), (2, 1)), "uv5": ((6, 1), (5, 1))}
ARCS = {k: walk(A) for k, A in MOTIONS.items()}
XT = {}
ALLARC = [(k, j, ar) for k in ARCS for j, ar in enumerate(ARCS[k])]
for (k1, j1, a1_), (k2, j2, a2_) in combinations(ALLARC, 2):
    p = cross(a1_["form"], a2_["form"])
    if p is not None:
        sg = a1_["s"] * a2_["s"] * sgn(a2_["center"] - a1_["center"])
        XT[(k1, j1, k2, j2)] = (kind(p), sg)
        XT[(k2, j2, k1, j1)] = (kind(p), -sg)


def mirror_data(k):
    A = MOTIONS[k]
    for ar_m in walk(cj(A)):
        for ar in ARCS[k]:
            for e, Nm in ((1, ar["N"]), (-1, pos(inv(ar["N"])))):
                if ar_m["N"] == Nm:
                    h = mm(ar_m["m"], inv(ar["m"]))
                    assert same(mm(mm(h, A if e == 1 else inv(A)), inv(h)), cj(A))
                    return e, h
    return None                                                          # the mirror image is another orbit of the modular surface


MIR = {k: mirror_data(k) for k in MOTIONS}


# ---------------------------------------------------------------- the curve at level N: darts, lifts, geodesics, crossings
class Level:
    def __init__(self, N, motions):
        self.N = N
        units = [u for u in range(1, N) if math.gcd(u, N) == 1]
        CAN = {}
        for c in range(N):
            for d in range(N):
                if math.gcd(math.gcd(c, d), N) == 1:
                    CAN[(c, d)] = min(((u * c) % N, (u * d) % N) for u in units)
        self.CAN, self.pts = CAN, sorted(set(CAN.values()))
        self.idx = {p: i for i, p in enumerate(self.pts)}
        self.lifts, self.geo = [], []
        gid = {}
        for k in motions:
            A, seen = MOTIONS[k], set()
            for x in range(len(self.pts)):
                if x in seen:
                    continue
                cyc = [x]
                while True:
                    y = self.act(cyc[-1], A)
                    if y == x:
                        break
                    cyc.append(y)
                seen.update(cyc)
                key = frozenset((self.act(xp, ar["m"]), unor(ar["form"])) for xp in cyc for ar in ARCS[k])
                L = {"motion": k, "cyc": cyc, "dart": self.pts[x], "ell": len(cyc), "key": key}
                if key not in gid:
                    gid[key] = len(self.geo)
                    self.geo.append({"motion": k, "ell": len(cyc), "rep": L, "lifts": [L]})
                else:
                    self.geo[gid[key]]["lifts"].append(L)
                L["geo"] = gid[key]
                self.lifts.append(L)
        self.lift_of = {(L["motion"], x): L for L in self.lifts for x in L["cyc"]}
        self.NG = len(self.geo)

    def act(self, i, g):
        c, d = self.pts[i]
        N = self.N
        return self.idx[self.CAN[((c * g[0][0] + d * g[1][0]) % N, (c * g[0][1] + d * g[1][1]) % N)]]

    def mirror(self):
        """For each geodesic: its mirror image and the orientation factor (None if the image is not among them)."""
        out = []
        for g, G in enumerate(self.geo):
            L = G["rep"]
            if MIR[L["motion"]] is None:
                out.append((None, 0))
                continue
            e, h = MIR[L["motion"]]
            c, d = self.pts[L["cyc"][0]]
            xbar = self.idx[self.CAN[(c % self.N, (-d) % self.N)]]
            tgt = self.lift_of[(L["motion"], self.act(xbar, h))]
            out.append((tgt["geo"], e * (1 if self.geo[tgt["geo"]]["rep"] is tgt else -1)))
        return out

    def crossings(self, conv="default"):
        inc = {"default": {"int", "left", "bot-"}, "mirrored": {"int", "right", "bot+"}}[conv]
        corner = {"default": "rho", "mirrored": "rho1"}[conv]
        NG = self.NG
        Im = [[0] * NG for _ in range(NG)]
        Io = [[0] * NG for _ in range(NG)]
        by_sheet, inst = {}, []
        for g, G in enumerate(self.geo):
            L = G["rep"]
            for xp in L["cyc"]:
                for j, ar in enumerate(ARCS[L["motion"]]):
                    it = (g, L["motion"], j, self.act(xp, ar["m"]))
                    inst.append(it)
                    by_sheet.setdefault(it[3], []).append(it)

        def add(g, h, sg):
            if g == h:
                Im[g][g] += 1
            else:
                Im[g][h] += 1; Im[h][g] += 1
                Io[g][h] += sg; Io[h][g] -= sg
        for s, L in by_sheet.items():
            sS = self.act(s, S)
            for p_, q_ in combinations(L, 2):
                hit = XT.get((p_[1], p_[2], q_[1], q_[2]))
                if hit is None:
                    continue
                loc, sg = hit
                if loc == "i":
                    if (s < sS) if conv == "default" else (s > sS):
                        add(p_[0], q_[0], sg)
                elif loc in inc:
                    add(p_[0], q_[0], sg)
        groups = {}
        for g, k, j, s in inst:
            ar = ARCS[k][j]
            if corner not in (kind(ar["entry"]), kind(ar["exit"])):
                continue
            Nm = ar["N"]
            if corner == "rho1":
                s, Nm = self.act(s, T), conj(T, Nm)
            orb = [s, self.act(s, R), self.act(s, R2)]
            u = min(orb)
            r = orb.index(u)
            Nm = conj(R if r == 1 else R2, Nm) if r else Nm
            a, b, c = form(Nm)
            assert a - Fr(b, 2) + c == 0
            groups.setdefault(u, []).append((g, (1 if a > 0 else -1), Fr(-b, 2 * a)))
        for u, L in groups.items():
            for (g, sa, ca), (h, sb, cb) in combinations(L, 2):
                add(g, h, sa * sb * sgn(cb - ca))
        return Im, Io


# ---------------------------------------------------------------- Part A at level 143
L143 = Level(143, ["cat", "s10", "pell", "uv3", "u2v2", "uv5"])
IM, IO = L143.crossings("default")
IM2, IO2 = L143.crossings("mirrored")
MIRROR = L143.mirror()
geo = L143.geo
# name the old geodesics by MP12's line numbers
line_geo = {}
for n_, l in enumerate(R12["lines"]):
    x = L143.idx[tuple(l["dart"])]
    line_geo[n_ + 1] = L143.lift_of[(l["motion"], x)]["geo"]
old = sorted(set(line_geo.values()))
new = [g for g in range(L143.NG) if g not in old]
P25, P26, P27, P28 = (line_geo[n] for n in (25, 26, 27, 28))
# A0: the old 24 reproduce MP13 (MP13 oriented each geodesic by its lowest-numbered line; here by the first lift found)
map13, eps13 = {}, {}
for G in R13["geodesics"]:
    l = R12["lines"][min(G["lines"]) - 1]
    Lft = L143.lift_of[(l["motion"], L143.idx[tuple(l["dart"])])]
    map13[G["id"]] = Lft["geo"]
    eps13[G["id"]] = 1 if geo[Lft["geo"]]["rep"] is Lft else -1
a0 = all(IM[map13[a]][map13[b]] == R13["I"][a][b] for a in range(24) for b in range(24)) and \
     all(R13["iota"][a][b] == eps13[a] * eps13[b] * IO[map13[a]][map13[b]] for a in range(24) for b in range(24))
check("A0", "the old 24 geodesics give MP13's crossing and signed matrices again, and the mirrored conventions give the same counts on all the lines, old and new", a0 and IM == IM2 and IO == IO2,
      f"{L143.NG} geodesics in all: 24 old, {len(new)} new")
kindm = {}
for g in range(L143.NG):
    tg, e = MIRROR[g]
    kindm[g] = "reversed" if (tg == g and e == -1) else ("kept" if (tg == g and e == 1) else ("swapped" if tg is not None else "none"))
old_kinds = {k_: sum(1 for g in old if kindm[g] == k_) for k_ in ("reversed", "kept", "swapped")}
lab = lambda g: f"{geo[g]['motion']}{geo[g]['ell']}" + (f"(lines {'-'.join(str(n) for n in sorted(n for n, gg in line_geo.items() if gg == g))})" if g in old else f"({geo[g]['rep']['dart'][0]}:{geo[g]['rep']['dart'][1]})")
newinfo = [(lab(g), kindm[g]) for g in new]
rev_old = [g for g in old if kindm[g] == "reversed"]
rev_new = [g for g in new if kindm[g] == "reversed"]
pell = [P25, P26, P27, P28]
a1 = len(rev_new) > 0 and all(IO[g][h] == 0 for g in rev_new for h in pell + rev_old + rev_new if g != h)
check("A1", "every new geodesic the mirror reverses has signed count zero with the four Pell lines, with the 10 old reversed geodesics and with the other reversed new ones", a1,
      f"old: {old_kinds}; new: {newinfo}")
a2 = all(IO[P25][d] == IO[P28][d] == -IO[P26][d] == -IO[P27][d] for d in new)
check("A2", "for every new geodesic d: iota(25, d) = iota(28, d) = -iota(26, d) = -iota(27, d), crossing by crossing", a2,
      "; ".join(f"{lab(d)}: {IO[P25][d]}, {IO[P26][d]}, {IO[P27][d]}, {IO[P28][d]}" for d in new))
nz = [(lab(d), kindm[d], IO[P25][d], IM[P25][d]) for d in new if IO[P25][d] != 0]
a3 = any(k_ != "reversed" for _, k_, _, _ in nz)
check("A3", "some new geodesic the mirror does not reverse has iota(25, d) != 0: the Pell class pairs with the rest of homology", a3, f"non-zero: {nz}")

# A4: the shadows on X_0(11) and X_0(13)
shadow = {}
for N_ in (11, 13):
    LN = Level(N_, ["cat", "s10", "pell"])
    Imn, Ion = LN.crossings("default")
    pg = LN.lift_of[("pell", LN.idx[LN.CAN[(0, 1)]])]["geo"]
    rows = {}
    for k in ("cat", "s10"):
        for g, G in enumerate(LN.geo):
            if G["motion"] == k and G["ell"] == 1:
                rows[f"{k} through ({G['rep']['dart'][0]}:{G['rep']['dart'][1]})"] = Imn[pg][g]
    shadow[N_] = {"pell_shadow_length": LN.geo[pg]["ell"], "crossings": rows,
                  "geodesics": [f"{G['motion']}{G['ell']}" for G in LN.geo]}
a4 = all(v == 0 for v in shadow[11]["crossings"].values() if True) and all(v == 0 for k_, v in shadow[13]["crossings"].items() if k_.startswith("s10"))
check("A4", "conjecture (low confidence): the Pell shadow through (0:1) misses the cat's and Sphereland 10's shadows of length 1 on X_0(11), and Sphereland 10's on X_0(13)", a4,
      f"X_0(11): {shadow[11]['crossings']}; X_0(13): {dict((k_, v) for k_, v in shadow[13]['crossings'].items() if k_.startswith('s10'))}")
NUM["A"] = {"geodesics": [{"id": g, "label": lab(g), "motion": geo[g]["motion"], "ell": geo[g]["ell"], "mirror": kindm[g], "dart": list(geo[g]["rep"]["dart"]),
                           "lifts": len(geo[g]["lifts"])} for g in range(L143.NG)],
            "I": IM, "iota": IO, "pell": pell, "new": new, "shadow": shadow,
            "new_with_pell": {lab(d): {"I": [IM[p][d] for p in pell], "iota": [IO[p][d] for p in pell]} for d in new}}

# ---------------------------------------------------------------- Part B: the Eisenstein series of Gamma_0(143)
mp.mp.dps = 40
xi = lambda s: mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def phi1(s):
    return mp.sqrt(mp.pi) * mp.gamma(s - mp.mpf(1) / 2) * mp.zeta(2 * s - 1) / (mp.gamma(s) * mp.zeta(2 * s))


def Mp(p, s):
    d = p ** (2 * s) - 1
    a, b = (p - 1) / d, (p ** s - p ** (1 - s)) / d
    return mp.matrix([[a, b], [b, a]])


def kron(A, B):
    out = mp.matrix(4, 4)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    out[2 * i + k, 2 * j + l] = A[i, j] * B[k, l]
    return out


Phi = lambda s: phi1(s) * kron(Mp(11, s), Mp(13, s))                  # cusps in the order (inf, inf), (inf, 0), (0, inf), (0, 0) at (11, 13)
h_ = mp.mpf(10) ** -25
mx = lambda M: max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))
# First run (its output written at 12:08:59 by the file times): the residue was taken at 40 digits with h = 1e-25, where 1 + 2h keeps only 15 digits of h,
# so zeta(1 + 2h) and the residue came out good to 7e-17 only, against a tolerance of 1e-20.  A fault of the check's
# precision, recorded in MP14; the residue is now taken at 80 digits.
with mp.workdps(80):
    hr = mp.mpf(10) ** -30
    res1 = Phi(1 + hr) * hr
    res_dev = max(abs(res1[i, j] - 3 / (168 * mp.pi)) for i in range(4) for j in range(4))
res_ok = res_dev < 1e-20
fe = max(mx(Phi(s) * Phi(1 - s) - mp.eye(4)) for s in (mp.mpf("0.31"), mp.mpc("0.7", "1.3"), mp.mpc("2.2", "-0.4")))
un = max(mx(Phi(mp.mpc(0.5, t)) * Phi(mp.mpc(0.5, t)).H - mp.eye(4)) for t in (0.3, 2.0, 7.5))
half = Phi(mp.mpf(1) / 2 + h_)
half_ok = max(abs(half[i, j] - (-1 if i == j else 0)) for i in range(4) for j in range(4)) < 1e-20
# the Dirichlet series over double cosets, summed directly (numpy), at s = 2.25 and 3
Cmax = 2000000
phi_e = np.arange(Cmax + 1, dtype=np.float64)
for p in range(2, Cmax + 1):
    if phi_e[p] == p:
        phi_e[p::p] *= (1 - 1 / p)
cs = np.arange(1, Cmax + 1)
ds_ok, ds_err = True, 0.0
for s in (2.25, 3.0):
    pref = float(mp.sqrt(mp.pi) * mp.gamma(s - 0.5) / mp.gamma(s))
    m143 = (cs % 143 == 0)
    s_inf = pref * np.sum(phi_e[1:][m143] * cs[m143] ** (-2 * s))
    cop = (cs % 11 != 0) & (cs % 13 != 0)
    s_0 = pref * 143 ** (-s) * np.sum(phi_e[1:][cop] * cs[cop] ** (-2 * s))
    P = Phi(mp.mpf(s))
    ds_err = max(ds_err, abs(s_inf - float(P[0, 0])) / abs(float(P[0, 0])), abs(s_0 - float(P[0, 3])) / abs(float(P[0, 3])))
check("B1", "the scattering matrix of Gamma_0(143): every entry has residue 3/(168 pi) at s = 1; Phi(s) Phi(1 - s) = I; unitary on Re s = 1/2; Phi(1/2) = -I; its (inf, inf) and (inf, 0) entries are the Dirichlet series over double cosets at s = 2.25 and 3",
      res_ok and fe < 1e-30 and un < 1e-30 and half_ok and ds_err < 1e-10, f"residue {float(res_dev):.1e}, functional equation {float(fe):.1e}, unitarity {float(un):.1e}, Dirichlet series {ds_err:.1e} (amended before the run from s = 1.5, where direct summation is too slow)")


def E1(z, s):
    """The level-1 series sum over Gamma_inf\\PSL(2, Z) of Im(gamma z)^s, by its Fourier expansion after reducing z to F."""
    _, w = reduce_point(mp.mpc(z))
    x, y = w.real, w.imag
    nu = s - mp.mpf(1) / 2
    out = y ** s + phi1(s) * y ** (1 - s)
    c = 4 / xi(2 * s)
    n, tot = 1, mp.mpc(0)
    while True:
        k = mp.besselk(nu, 2 * mp.pi * n * y)
        term = n ** nu * mp.fsum(mp.mpf(d) ** (1 - 2 * s) for d in range(1, n + 1) if n % d == 0) * mp.sqrt(y) * k * mp.cos(2 * mp.pi * n * x)
        tot += term
        if n > 3 and abs(k) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
        n += 1
    return out + c * tot


def Eis(z, s):
    """The four Eisenstein series of Gamma_0(143), cusps inf, 0, 1/11, 1/13."""
    e1, e11, e13, e143 = (E1(d * z, s) for d in (1, 11, 13, 143))
    a, b = mp.mpf(11) ** s, mp.mpf(13) ** s
    den = (a * a - 1) * (b * b - 1)
    return {"inf": (a * b * e143 - a * e11 - b * e13 + e1) / den, "0": (a * b * e1 - a * e13 - b * e11 + e143) / den,
            "1/11": (a * b * e11 - a * e143 - b * e1 + e13) / den, "1/13": (a * b * e13 - a * e1 - b * e143 + e11) / den}


Z5 = [mp.mpc("0.123", "0.6"), mp.mpc("-0.31", "0.05"), mp.mpc("0.4", "1.7"), mp.mpc("0.0613", "0.0291"), mp.mpc("-0.07", "0.2")]
ratios, b2_ok = {}, True
for z in Z5:
    rs = []
    for hh in (mp.mpf(10) ** -6, mp.mpf(10) ** -9, mp.mpf(10) ** -12):
        E = Eis(z, mp.mpf(1) / 2 + hh)
        rs.append(max(abs(E[a]) / hh for a in E))
    ratios[str(complex(z))] = [float(v) for v in rs]
    b2_ok &= rs[-1] < 1e6 and abs(rs[-1] - rs[-2]) < 1e-3 * (1 + rs[-1])
check("B2", "at s = 1/2 + h every one of the four series is h times a bounded function (h = 1e-6, 1e-9, 1e-12, five points): all of them vanish identically at the threshold", b2_ok,
      "; ".join(f"{k}: {v[-1]:.4g}" for k, v in ratios.items()))
r0 = 1 / mp.sqrt(143)
diff_on, diff_off = mp.mpf(0), mp.mpf(10)
for s in (mp.mpf(2), mp.mpf("0.75"), mp.mpc("0.5", "2")):
    for th in (mp.mpf("0.4"), mp.mpf("1.1"), mp.mpf("2.3")):
        z = r0 * mp.expj(th)
        E = Eis(z, s)
        scale = max(abs(E[a]) for a in E)
        diff_on = max(diff_on, abs(E["inf"] - E["0"]) / scale, abs(E["1/11"] - E["1/13"]) / scale)
        z2 = 1.3 * z
        E2 = Eis(z2, s)
        diff_off = min(diff_off, abs(E2["inf"] - E2["0"]) / max(abs(E2[a]) for a in E2))
check("B3", "on the Pell line |z| = 1/sqrt(143), E_inf = E_0 and E_1/11 = E_1/13 at s = 2, 0.75 and 1/2 + 2i; off the line they differ", diff_on < 1e-25 and diff_off > 1e-3,
      f"on the line {float(diff_on):.1e}; off it at least {float(diff_off):.2e}")
check("B4", "Selberg's eigenvalue conjecture is a theorem for Gamma_1(N), N squarefree and below 857 (Booker and Strombergsson 2007), and Gamma_0(143) contains Gamma_1(143): no eigenvalue of X_0(143) lies in (0, 1/4) [literature, not computed]", True)
NUM["B"] = {"functional_eq": float(fe), "unitarity": float(un), "dirichlet": ds_err, "ratios_at_half": ratios, "pell_line": float(diff_on), "off_line": float(diff_off)}
# values for Route B: the four series at s = 2 at two points
PTS = [mp.mpc("0.3", "0.9"), mp.mpc("-0.1", "0.05")]
NUM["B"]["s2"] = [{"z": [float(z.real), float(z.imag)], "E": {a: float(v.real) for a, v in Eis(z, mp.mpf(2)).items()}, "E1": float(E1(z, mp.mpf(2)).real)} for z in PTS]

print("\nfound, not predicted")
for k in ("uv3", "u2v2", "uv5"):
    gs = [g for g in new if geo[g]["motion"] == k]
    print(f"  {k}: {len(ARCS[k])} arcs a period; lifts {sorted(L['ell'] for L in L143.lifts if L['motion'] == k)}; geodesics {len(gs)}: " +
          ", ".join(f"{lab(g)} {kindm[g]}" for g in gs))
print("  new geodesics against the Pell lines 25, 26, 27, 28 (crossings; signed):")
for d in new:
    print(f"   {lab(d):22s} {kindm[d]:9s} I = {[IM[p][d] for p in pell]}  iota = {[IO[p][d] for p in pell]}")
print("  self-crossings of the new geodesics:", {lab(d): IM[d][d] for d in new})
print("  shadows:", json.dumps(shadow)[:600])
npass = sum(1 for r in RESULTS if r[2])
print(f"\n{npass} of {len(RESULTS)} checks passed in {time.time() - t_start:.0f} s")
for tag, name, okk, det in RESULTS:
    if not okk:
        print("  MISS", tag, name, det)
if len(sys.argv) > 3:
    with open(sys.argv[3], "w") as fh:
        json.dump({"checks": [{"tag": t, "name": n, "pass": o, "detail": d} for t, n, o, d in RESULTS], "numbers": NUM}, fh, indent=1)
sys.exit(0 if npass == len(RESULTS) else 1)
