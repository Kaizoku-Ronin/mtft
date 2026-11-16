#!/usr/bin/env python3
"""
mp13_tg9_verify.py  --  Route A for MP13 (test T-G9 of the gravity map, piece II.8)

Which unions of 12.7's closed lines bound in M = Gamma_0(143)\\PSL(2,R), and what
are their linking numbers?  For a closed geodesic g of the curve, the orbit and
the reversed orbit together, K_g, bound, and (MP13, derived before the run)

        Lk(K_g, K_h) = -I(g, h),

I being the number of crossings of g and h on the curve.  This script finds every
crossing.  Each closed geodesic of the modular surface is cut along the standard
fundamental domain F in exact rational arithmetic (a geodesic in F-coordinates is
a(x^2 + y^2) + bx + c = 0 with integer a, b, c; two of them meet at rational x and
rational y^2).  The 168 sheets of the curve over F are the points of P^1(Z/143), and
the arc j of a lift through the dart x lies, in its k-th period, in the sheet
x A^k m_j.  Crossings are found sheet by sheet, with a half-open fundamental domain;
the points over i and over rho, where several sheets meet, have rules of their own.

Predictions P1 to P7 are in MP13_TG9_Plan.md and were written before the first run.
Route B: test_core_tg9.js (a formula on the darts, on the page's core).  Route C:
mp13_lorenz.py (Ghys's Lorenz template and the transfer to level 1).

Run:  python3 mp13_tg9_verify.py ../mp12/mp12_numbers.json mp13_numbers.json      (about two seconds)
Needs sympy and mpmath.  Nothing is imported from the lab or from mtft.

Outcome on 12 November 2026: 8 of 8 on the first run (about 11:24 New York time).
Before that run, on reading, the frame rotation at the lifts of rho was found to be
the wrong way round (it would have touched only the signs there) and was corrected.
"""
import sys, json, math, time
from fractions import Fraction as Fr
from itertools import combinations
import mpmath as mp
import sympy as sy

t_start = time.time()
N = 143
RESULTS, NUM = [], {}
R12 = json.load(open(sys.argv[1]))["numbers"]


def check(tag, name, cond, detail=""):
    RESULTS.append((tag, name, bool(cond), detail))
    print(("ok   " if cond else "FAIL ") + f"{tag:4s} {name}" + (f"   [{detail}]" if detail else ""), flush=True)


# ---------------------------------------------------------------- 2 x 2 integer matrices, PSL(2, Z)
def mm(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))


def inv(A):
    return ((A[1][1], -A[0][1]), (-A[1][0], A[0][0]))


def neg(A):
    return ((-A[0][0], -A[0][1]), (-A[1][0], -A[1][1]))


def pos(A):                                     # the sign with positive trace (hyperbolic elements only)
    return A if A[0][0] + A[1][1] > 0 else neg(A)


def same(A, B):
    return A == B or A == neg(B)


def conj(g, A):                                 # g^-1 A g: the geodesic of A seen in the frame z' = g^-1 z
    return pos(mm(mm(inv(g), A), g))


tpow = lambda n: ((1, n), (0, 1))
I2, S, T, Ti = ((1, 0), (0, 1)), ((0, -1), (1, 0)), ((1, 1), (0, 1)), ((1, -1), (0, 1))
R = ((0, 1), (-1, -1))                          # R(z) = -1/(z + 1) fixes rho; T = S R
R2 = mm(R, R)
assert same(mm(S, R), T) and same(mm(R, R2), I2)
HALF, Q34 = Fr(1, 2), Fr(3, 4)


# ---------------------------------------------------------------- A: one arc in F, exactly
def form(A):                                    # A = pos(A): the geodesic a(x^2 + y^2) + bx + c = 0 of its axis
    return (A[1][0], A[1][1] - A[0][0], -A[0][1])


def arc(A):
    """The arc of the oriented axis of A inside F: exact endpoints (x, y^2), or None if it has no length."""
    a, b, c = form(A)
    pts = set()
    yl = Fr(b, 2 * a) - Fr(c, a) - Fr(1, 4)                             # on x = -1/2
    if yl >= Q34:
        pts.add((-HALF, yl))
    yr = Fr(-b, 2 * a) - Fr(c, a) - Fr(1, 4)                            # on x = 1/2
    if yr >= Q34:
        pts.add((HALF, yr))
    if b != 0:                                                          # on |z| = 1
        xb = Fr(-(a + c), b)
        if abs(xb) <= HALF:
            pts.add((xb, 1 - xb * xb))
    pts = sorted(pts)
    if len(pts) != 2:
        return None
    s = 1 if a > 0 else -1                                              # the attracting end is to the right when a > 0
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
CAND = {"rho1": [T, mm(T, R), mm(T, R2), mm(mm(T, R), Ti), mm(mm(T, R2), Ti)],      # the other sectors round rho + 1
        "rho": [R, R2, Ti, mm(R, Ti), mm(R2, Ti)]}                                  # and round rho
mob = lambda g, z: (g[0][0] * z + g[0][1]) / (g[1][0] * z + g[1][1])


def step(ar):
    k = kind(ar["exit"])
    x, y2 = ar["exit"]
    if k in ("right", "left", "i", "bot-", "bot+"):
        g = {"right": T, "left": Ti}.get(k, S)
        nxt = arc(conj(g, ar["N"]))
        want = (-x, y2)                                                 # T, T^-1 and S all send the exit to the next entry this way
        assert nxt is not None and nxt["entry"] == want, (ar, nxt)
        return g, nxt
    good = []                                                           # through a corner: into the opposite sector
    for g in CAND[k]:
        nxt = arc(conj(g, ar["N"]))
        if nxt is None or kind(nxt["entry"]) not in ("rho", "rho1"):
            continue
        if abs(mob(g, CORNER[kind(nxt["entry"])]) - CORNER[k]) < 1e-9:
            good.append((g, nxt))
    assert len(good) == 1, (ar, len(good))
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
    D = b * b - 4 * a * c
    mp.mp.dps = 60
    xc, r = mp.mpf(-b) / (2 * a), mp.sqrt(D) / (2 * abs(a))
    m0, w = reduce_point(xc + r * mp.expjpi(mp.mpf("0.3141")))        # a generic point of the axis
    ar = arc(conj(m0, A))
    assert ar is not None
    arcs, m = [], m0
    while True:
        ar["m"] = m
        arcs.append(ar)
        g, ar = step(ar)
        m = mm(m, g)
        if ar["N"] == arcs[0]["N"]:
            break
    assert same(mm(m, inv(m0)), A), "the walk did not close after one period"
    return arcs, m0, m


def arc_length(ar):
    a, b, c = ar["form"]
    out = []
    for p in (ar["entry"], ar["exit"]):
        x, y2 = p
        out.append(mp.atan2(mp.sqrt(mp.mpf(y2.numerator) / y2.denominator), mp.mpf(x.numerator) / x.denominator - mp.mpf(-b) / (2 * a)))
    return abs(mp.log(mp.tan(out[1] / 2) / mp.tan(out[0] / 2)))


MOTIONS = {"cat": ((2, 1), (1, 1)), "s10": ((84, -47), (143, -80)), "pell": ((12, 1), (143, 12))}
WALK = {k: walk(A) for k, A in MOTIONS.items()}
ARCS = {k: v[0] for k, v in WALK.items()}
info, len_err = {}, 0.0
for k, arcs in ARCS.items():
    t = MOTIONS[k][0][0] + MOTIONS[k][1][1]
    L = 2 * mp.acosh(mp.mpf(abs(t)) / 2)
    len_err = max(len_err, float(abs(sum(arc_length(ar) for ar in arcs) - L)))
    exits = [kind(ar["exit"]) for ar in arcs]
    info[k] = {"arcs": len(arcs), "through_i": exits.count("i"), "through_rho": exits.count("rho") + exits.count("rho1"),
               "exits": "".join({"right": "R", "left": "L", "i": "i", "bot-": "B", "bot+": "B", "rho": "p", "rho1": "q"}[e] for e in exits)}
check("P1", "exact cutting sequences: the cat orbit has four arcs a period and passes through i twice and never through rho; Sphereland 10's passes through rho and not through i; the Pell orbit not through i; each walk closes after one period with +-A m0, and the arcs' lengths add up to 2 arccosh(t/2)",
      info["cat"]["arcs"] == 4 and info["cat"]["through_i"] == 2 and info["cat"]["through_rho"] == 0 and info["s10"]["through_rho"] >= 1 and info["s10"]["through_i"] == 0
      and info["pell"]["through_i"] == 0 and len_err < 1e-12, f"{info}; worst length error {len_err:.1e}")
NUM["P1"] = {"walks": info, "length_err": len_err,
             "arcs": {k: [{"form": list(ar["form"]), "m": [list(r) for r in ar["m"]], "entry": [str(v) for v in ar["entry"]], "exit": [str(v) for v in ar["exit"]]} for ar in arcs] for k, arcs in ARCS.items()}}

# ---------------------------------------------------------------- B: the darts, the 28 lines, the 24 geodesics
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


def unor(f):                                                            # the unoriented semicircle
    return f if f[0] > 0 else tuple(-v for v in f)


LINES = R12["lines"]
for n_, l in enumerate(LINES):
    A = MOTIONS[l["motion"]]
    x = idx[tuple(l["dart"])]
    cyc = [x]
    while True:
        y = act(cyc[-1], A)
        if y == x:
            break
        cyc.append(y)
    assert len(cyc) == l["ell"]
    l["num"], l["cyc"] = n_ + 1, cyc
    l["key"] = frozenset((act(xp, ar["m"]), unor(ar["form"])) for xp in cyc for ar in ARCS[l["motion"]])
GEO, gid = [], {}
for l in LINES:
    if l["key"] not in gid:
        gid[l["key"]] = len(GEO)
        GEO.append({"rep": l, "lines": [l["num"]], "motion": l["motion"], "ell": l["ell"]})
    else:
        GEO[gid[l["key"]]]["lines"].append(l["num"])
    l["geo"] = gid[l["key"]]
NG = len(GEO)
pairs_cat = sorted(tuple(g["lines"]) for g in GEO if g["motion"] == "cat")
cls_opp = all(LINES[a - 1]["cls"] == [-v for v in LINES[b - 1]["cls"]] for a, b in pairs_cat)
check("P1", "the 28 lines are 24 geodesics: each cat lift and its partner through S trace one geodesic in opposite directions (pairs 1-4, 2-3, 5-6, 7-8, with opposite classes); the 16 + 4 others are distinct",
      NG == 24 and pairs_cat == [(1, 4), (2, 3), (5, 6), (7, 8)] and cls_opp, f"cat pairs {pairs_cat}")

# the arcs of each geodesic, once round: its representative's ell periods
INST, by_sheet = [], {}
for g, G in enumerate(GEO):
    l = G["rep"]
    for xp in l["cyc"]:
        for j, ar in enumerate(ARCS[l["motion"]]):
            s = act(xp, ar["m"])
            inst = {"geo": g, "motion": l["motion"], "j": j, "sheet": s}
            INST.append(inst)
            by_sheet.setdefault(s, []).append(inst)
per_sheet = {s: sorted((i["motion"]) for i in L) for s, L in by_sheet.items()}
count_ok = len(by_sheet) == ND and all(L.count("cat") == 2 and L.count("s10") == info["s10"]["arcs"] and L.count("pell") == info["pell"]["arcs"] for L in per_sheet.values())
check("P1", "every sheet carries two cat arcs (one of each half of the orbit), and every arc of Sphereland 10's and of the Pell orbit once", count_ok)

# ---------------------------------------------------------------- C: crossings, exactly
ALLARC = [(k, j, ar) for k in ARCS for j, ar in enumerate(ARCS[k])]
sgn = lambda v: (v > 0) - (v < 0)


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


XT = {}
for (k1, j1, a1_), (k2, j2, a2_) in combinations(ALLARC, 2):
    p = cross(a1_["form"], a2_["form"])
    if p is not None:
        sg = a1_["s"] * a2_["s"] * sgn(a2_["center"] - a1_["center"])
        XT[(k1, j1, k2, j2)] = (kind(p), sg, p)
        XT[(k2, j2, k1, j1)] = (kind(p), -sg, p)
locs_seen = {}
for v in XT.values():
    locs_seen[v[0]] = locs_seen.get(v[0], 0) + 1

CONV = {"default": {"in": {"int", "left", "bot-"}, "i": lambda s: s < act(s, S), "corner": "rho"},
        "mirrored": {"in": {"int", "right", "bot+"}, "i": lambda s: s > act(s, S), "corner": "rho1"}}


def crossings(cv):
    Im = [[0] * NG for _ in range(NG)]
    Io = [[0] * NG for _ in range(NG)]
    where = {}

    def add(g, h, sg, loc):
        where[(min(g, h), max(g, h), loc)] = where.get((min(g, h), max(g, h), loc), 0) + 1
        if g == h:
            Im[g][g] += 1
        else:
            Im[g][h] += 1; Im[h][g] += 1
            Io[g][h] += sg; Io[h][g] -= sg
    for s, L in by_sheet.items():
        for p_, q_ in combinations(L, 2):
            hit = XT.get((p_["motion"], p_["j"], q_["motion"], q_["j"]))
            if hit is None:
                continue
            loc, sg, _ = hit
            if loc == "i":
                if cv["i"](s):
                    add(p_["geo"], q_["geo"], sg, loc)
            elif loc in cv["in"]:
                add(p_["geo"], q_["geo"], sg, loc)
    groups = {}                                                         # the corner rule: all branches through a lift of rho
    for it in INST:
        ar = ARCS[it["motion"]][it["j"]]
        ends = {kind(ar["entry"]), kind(ar["exit"])}
        if cv["corner"] not in ends:
            continue
        s, Nm = it["sheet"], ar["N"]
        if cv["corner"] == "rho1":                                      # rho + 1 in the sheet s is rho in the sheet s T
            s, Nm = act(s, T), conj(T, Nm)
        orb = [s, act(s, R), act(s, R2)]
        u = min(orb)
        r = orb.index(u)                                                # u = s R^r, so s = u R^-r and the arc is R^-r of itself in u's frame
        Nm = conj(R if r == 1 else R2, Nm) if r else Nm                 # R^-r N R^r = conj(R^r, N)
        a, b, c = form(Nm)
        assert a - Fr(b, 2) + c == 0                                    # it passes through rho
        groups.setdefault(u, []).append((it["geo"], Nm))
    for u, L in groups.items():
        for (g, Na), (h, Nb) in combinations(L, 2):
            fa, fb = arc_like(Na), arc_like(Nb)
            sg = fa[0] * fb[0] * sgn(fb[1] - fa[1])
            add(g, h, sg, cv["corner"])
    return Im, Io, where, groups


def arc_like(Nm):                                                       # direction and centre of a geodesic through rho
    a, b, c = form(Nm)
    return (1 if a > 0 else -1, Fr(-b, 2 * a))


IM, IO, WHERE, GROUPS = crossings(CONV["default"])
IM2, IO2, WHERE2, _ = crossings(CONV["mirrored"])
check("P6", "the mirrored half-open conventions (right edge, bottom arc x > 0, the higher sheet at i, the corner rho + 1) give the same crossing counts and the same signs", IM == IM2 and IO == IO2,
      f"crossing points by kind of place: {dict(sorted(locs_seen.items()))} pairs of arcs on the modular surface")

# ---------------------------------------------------------------- D: the predictions on the counts
cat = [g for g in range(NG) if GEO[g]["motion"] == "cat"]
s10 = [g for g in range(NG) if GEO[g]["motion"] == "s10"]
pell = [g for g in range(NG) if GEO[g]["motion"] == "pell"]
cat_where = {k: v for k, v in WHERE.items() if k[0] in cat and k[1] in cat}
cat_total = sum(cat_where.values())
row_ok = all(sum(IM[P][Q] for Q in cat if Q != P) + 2 * IM[P][P] == 2 * GEO[P]["ell"] for P in cat)
check("P3", "the four cat geodesics cross only at lifts of i, 84 crossings in all, and each geodesic P has sum over Q != P of I(P, Q) + 2 self(P) = 2 ell_P (14 and 70)",
      set(k[2] for k in cat_where) == {"i"} and cat_total == 84 and row_ok,
      "cat block (rows: lines " + ", ".join("-".join(map(str, GEO[g]["lines"])) for g in cat) + "): " + str([[IM[P][Q] for Q in cat] for P in cat]))
tot = lambda A_, B_: sum(IM[a][b] for a in A_ for b in B_)
T_cs, T_cp, T_sp = tot(cat, s10), tot(cat, pell), tot(s10, pell)
check("P4", "totals of crossings on the curve: cat x Sphereland 10 = 336, cat x Pell = 336, Sphereland 10 x Pell = 1008 (168 times the Lorenz sums 2, 2, 6)", (T_cs, T_cp, T_sp) == (336, 336, 1008), f"found {T_cs}, {T_cp}, {T_sp}")

# P5: the signed counts are a bilinear form in the 26 periods
names = ("ell", "old", "q4", "q6")
Cmat = sy.Matrix([[v for nm in names for v in GEO[g]["rep"]["blocks"][nm]] for g in range(NG)])
rel = (Cmat.T).nullspace()
IOm = sy.Matrix(IO)
rel_ok = all((v.T * IOm).is_zero_matrix for v in rel)
check("P5", "the period vectors of the 24 geodesics have rank 14, and each of their 10 linear relations annihilates the rows of the signed crossing matrix", Cmat.rank() == 14 and len(rel) == 10 and rel_ok, f"rank {Cmat.rank()}, {len(rel)} relations; rank of the signed matrix {IOm.rank()}")

# P7: the mirror and the Atkin-Lehner involutions
cj = lambda g: ((g[0][0], -g[0][1]), (-g[1][0], g[1][1]))


def mirror_data(k):
    A = MOTIONS[k]
    arcs_m, m0m, _ = walk(cj(A))
    for ar_m in arcs_m:
        for ar in ARCS[k]:
            for e, Nm in ((1, ar["N"]), (-1, pos(inv(ar["N"])))):
                if ar_m["N"] == Nm:
                    h = mm(ar_m["m"], inv(ar["m"]))
                    Ae = A if e == 1 else inv(A)
                    assert same(mm(mm(h, Ae), inv(h)), cj(A))
                    return e, h
    raise RuntimeError


MIR = {k: mirror_data(k) for k in MOTIONS}
line_of = {}
for l in LINES:
    for x in l["cyc"]:
        line_of[(l["motion"], x)] = l
img_m, eps_m = [None] * NG, [0] * NG
for g, G in enumerate(GEO):
    l = G["rep"]
    c, d = pts[l["cyc"][0]]
    xbar = idx[CAN[(c % N, (-d) % N)]]
    e, h = MIR[l["motion"]]
    tgt = line_of[(l["motion"], act(xbar, h))]
    img_m[g] = tgt["geo"]
    eps_m[g] = e * (1 if GEO[tgt["geo"]]["rep"] is tgt else -1)
mir_ok = sorted(img_m) == list(range(NG)) and all(IM[img_m[g]][img_m[h]] == IM[g][h] and eps_m[g] * eps_m[h] * IO[img_m[g]][img_m[h]] == -IO[g][h] for g in range(NG) for h in range(NG) if g != h) \
    and all(IM[img_m[g]][img_m[g]] == IM[g][g] for g in range(NG))
al_maps = {}
for i, mo, ell, Qd, hits in R12["post"]["al_images"]:
    l = LINES[i]
    if GEO[l["geo"]]["rep"] is not l or not hits:
        continue
    sg, j = hits[0]
    tgt = LINES[j]
    al_maps.setdefault(Qd, {})[l["geo"]] = (tgt["geo"], sg * (1 if GEO[tgt["geo"]]["rep"] is tgt else -1))
al_ok, al_pairs = True, 0
for Qd, mp_ in al_maps.items():
    for g, (wg, eg) in mp_.items():
        al_ok &= IM[wg][wg] == IM[g][g]
        for h, (wh, eh) in mp_.items():
            if g < h:
                al_pairs += 1
                al_ok &= IM[wg][wh] == IM[g][h] and eg * eh * IO[wg][wh] == IO[g][h]
gA, gB = gid[LINES[0]["key"]], gid[LINES[1]["key"]]
w11_swap = al_maps.get(11, {}).get(gA, (None,))[0] == gB and IM[gA][gA] == IM[gB][gB]
check("P7", "the mirror permutes the 24 geodesics keeping every crossing count and reversing every sign; W_11, W_13, W_143 keep the counts and the signs on every pair whose images are among them; W_11 exchanges the two cat geodesics of length 7, which have the same number of self-crossings",
      mir_ok and al_ok and w11_swap and al_pairs >= 10, f"mirror {[img_m[g] for g in range(NG)]}; {al_pairs} pairs under the Atkin-Lehner involutions; self-crossings of the two: {IM[gA][gA]}, {IM[gB][gB]}")

# ---------------------------------------------------------------- found, not predicted
LK = [[(-IM[g][h] if g != h else None) for h in range(NG)] for g in range(NG)]
lab = lambda g: f"{GEO[g]['motion']}{GEO[g]['ell']}(" + "-".join(map(str, GEO[g]["lines"])) + ")"
print("\nfound, not predicted")
print("  the geodesics:", "; ".join(f"{g}: {lab(g)}" for g in range(NG)))
print("  self-crossings:", [IM[g][g] for g in range(NG)])
print("  linking numbers Lk(K_g, K_h) = -I(g, h):")
for g in range(NG):
    print(f"   {lab(g):18s}", " ".join(f"{'.' if g == h else -IM[g][h]:>4}" for h in range(NG)))
unl = [(lab(g), lab(h)) for g in range(NG) for h in range(g + 1, NG) if IM[g][h] == 0]
print(f"  unlinked pairs: {len(unl)} of {NG * (NG - 1) // 2}")
alg = [[IO[g][h] for h in range(NG)] for g in range(NG)]
gap = sum(IM[g][h] - abs(IO[g][h]) for g in range(NG) for h in range(g + 1, NG))
print("  signed (algebraic) intersection numbers, rows as above:")
for g in range(NG):
    print(f"   {lab(g):18s}", " ".join(f"{'.' if g == h else IO[g][h]:>4}" for h in range(NG)))
print(f"  crossings beyond the algebraic count (sum of I - |iota| over pairs): {gap}")
print("  crossing places on the curve (pairs of geodesics), by kind:", {k: sum(v for kk, v in WHERE.items() if kk[2] == k) for k in sorted(set(kk[2] for kk in WHERE))})
lk_cat_B = {k: -sum(IM[P][d] for P in cat for d in grp) // 2 for k, grp in (("s10", s10), ("pell", pell))}
print("  Lk(all cat lifts, all lifts of B) = -(1/2) N:", lk_cat_B)
NUM["geodesics"] = [{"id": g, "motion": GEO[g]["motion"], "ell": GEO[g]["ell"], "lines": GEO[g]["lines"]} for g in range(NG)]
NUM["I"], NUM["iota"], NUM["self"] = IM, IO, [IM[g][g] for g in range(NG)]
NUM["Lk"] = LK
NUM["totals"] = {"cat_s10": T_cs, "cat_pell": T_cp, "s10_pell": T_sp, "Lk_cat_B": lk_cat_B}
NUM["places"] = {f"{k[0]}-{k[1]}-{k[2]}": v for k, v in WHERE.items()}
NUM["mirror"] = {"image": img_m, "eps": eps_m}
NUM["darts"] = {"cat_cycles": [{"lines": GEO[g]["lines"], "cycles": [[list(pts[x]) for x in LINES[n - 1]["cyc"]] for n in GEO[g]["lines"]]} for g in cat]}

npass = sum(1 for r in RESULTS if r[2])
print(f"\n{npass} of {len(RESULTS)} checks passed in {time.time() - t_start:.0f} s")
for tag, name, okk, det in RESULTS:
    if not okk:
        print("  MISS", tag, name, det)
if len(sys.argv) > 2:
    with open(sys.argv[2], "w") as fh:
        json.dump({"checks": [{"tag": t, "name": n, "pass": o, "detail": d} for t, n, o, d in RESULTS], "numbers": NUM}, fh, indent=1)
sys.exit(0 if npass == len(RESULTS) else 1)
