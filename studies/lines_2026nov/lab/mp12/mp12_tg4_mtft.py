#!/usr/bin/env python3
"""
mp12_tg4_mtft.py  --  Route B for MP12 (test T-G4) on the programme's own package.

mtft.hecke builds its own Manin-symbol model of X_0(143) and its own four Hecke blocks.  This script sends the 28
closed lines of Route A (their matrices, from mp12_numbers.json) through that model and compares what does not
depend on a choice of basis: where a block component vanishes, the rank of the 28 components in each block, and
the linear relations among the 28 lines, block by block and for the whole cuspidal class.

Run:  pip install mtft sympy ;  python3 mp12_tg4_mtft.py mp12_numbers.json
 or:  python3 mp12_tg4_mtft.py mp12_numbers.json /path/to/site-packages        (mtft 1.0.0 used)
"""
import sys, json
from fractions import Fraction as Fr
if len(sys.argv) > 2:
    sys.path.insert(0, sys.argv[2])
import sympy as sy
import mtft
import mtft.hecke as H

R = json.load(open(sys.argv[1]))["numbers"]
res = []


def check(name, cond, detail=""):
    res.append(bool(cond))
    print(("ok   " if cond else "FAIL ") + name + (f"   [{detail}]" if detail else ""), flush=True)


m = H.model()
N, E = m["N"], m["E"]
idx, canon, eid, esign, DEL = m["idx"], m["canon"], m["eid"], m["esign"], m["DEL"]


def symbol(g):
    """{infinity, g(infinity)} as a vector on mtft's 84 edges: Manin symbols ((-1)^(j-1) q_j : q_(j-1)) of the convergents."""
    (a, b), (c, d) = g
    if c < 0:
        a, c = -a, -c
    v = [0] * E

    def put(u, w):
        y = idx[canon(u % N, w % N)]
        v[eid[y]] += esign[y]
    if c:
        pm, qm, p_, q_ = 1, 0, a // c, 1
        num, den = a - (a // c) * c, c
        put(-1, 0)
        j = 0
        while num:
            j += 1
            aj = den // num
            num, den = den - aj * num, num
            pm, qm, p_, q_ = p_, q_, aj * p_ + pm, aj * q_ + qm
            put((-1) ** (j - 1) * q_, qm)
        assert (p_, q_) == (a, c)
    return v


blocks = H.blocks()
names = ("ell", "old", "q4", "q6")
W = sy.Matrix([[sy.Rational(x.numerator, x.denominator) for x in vec] for nm in names for vec in blocks[nm]]).T      # 26 x 26
Winv = W.inv()
dims = [len(blocks[nm]) for nm in names]
lines = R["lines"]
comp, cusp_ok = [], True
for l in lines:
    g = tuple(tuple(int(x) for x in row) for row in l["gamma"])
    v = symbol(g)
    cusp_ok &= all(sum(v[k] * DEL[k][c] for k in range(E)) == 0 for c in range(4))
    q = H._qcoords(m, [Fr(x) for x in v])
    w26 = sy.Matrix([sy.Rational(q[fj].numerator, q[fj].denominator) for fj in m["freeK"]])
    # the class must lie in the span of mtft's cuspidal basis K
    back = [sum(m["K"][j][i] * q[m["freeK"][j]] for j in range(26)) for i in range(m["nq"])]
    cusp_ok &= all(back[i] == q[i] for i in range(m["nq"]))
    comp.append(Winv * w26)
check("mtft 1.0.0: dimensions 84, 29, 26 and blocks 2, 4, 8, 12; every closed line is a cuspidal class in its model", (E, m["nq"], len(m["K"])) == (84, 29, 26) and dims == [2, 4, 8, 12] and cusp_ok, f"mtft {mtft.__version__}")
g4 = [1, 5, -1, -3, 1]; h6 = [-12, -7, 24, 2, -10, 0, 1]
check("mtft's quartic and sextic are the factors Route A found", list(H.G4) == g4 and list(H.H6) == h6 and "x**4 - 3*x**3 - x**2 + 5*x + 1" in R["P8"]["charpoly_T2"] and "x**6 - 10*x**4 + 2*x**3 + 24*x**2 - 7*x - 12" in R["P8"]["charpoly_T2"])


def relations(M):
    """Reduced basis of the linear relations among the rows of M."""
    ns = M.T.nullspace()
    return sy.Matrix.hstack(*ns).T.rref()[0] if ns else sy.zeros(0, M.rows)


off, same_zero, same_rank, same_rel, detail = 0, True, True, True, []
for nm, dm in zip(names, dims):
    Mb = sy.Matrix([[c[off + i] for i in range(dm)] for c in comp]); off += dm
    Ma = sy.Matrix([l["blocks"][nm] for l in lines])
    zb = [i for i in range(len(lines)) if Mb.row(i).is_zero_matrix]
    za = [i for i in range(len(lines)) if Ma.row(i).is_zero_matrix]
    same_zero &= za == zb
    same_rank &= Ma.rank() == Mb.rank()
    same_rel &= relations(Ma) == relations(Mb)
    detail.append(f"{nm}: rank {Mb.rank()}, {len(zb)} zero")
check("block by block: the same lines have zero component, the 28 components have the same rank, and the same linear relations hold", same_zero and same_rank and same_rel, "; ".join(detail))
Mb_all = sy.Matrix([list(c) for c in comp])
Ma_all = sy.Matrix([[v for nm in names for v in l["blocks"][nm]] for l in lines])
check("the whole cuspidal class: same rank and same relations among the 28 lines", Ma_all.rank() == Mb_all.rank() and relations(Ma_all) == relations(Mb_all), f"rank {Mb_all.rank()}")
# the transfer and the Pell lines, read in mtft's model alone
tr = all(sum((comp[i] for i, l in enumerate(lines) if l["motion"] == mo), sy.zeros(26, 1)).is_zero_matrix for mo in ("cat", "s10", "pell"))
pell = [i for i, l in enumerate(lines) if l["motion"] == "pell"]
check("in mtft's model alone: the lifts of each motion add up to zero; the four Pell lines have no 143a1 and no ghost component", tr and all(all(comp[i][j] == 0 for j in range(6)) for i in pell))

# ---------------------------------------------------------------- post hoc (after Route A's fourth run; not predictions)
# The finer transfer, the mirror and the Atkin-Lehner images, read in mtft's model alone.  The orbits modulo 11 and 13
# are recomputed here from each line's dart; the images are recomputed from the large matrices.
MOT = {"cat": ((2, 1), (1, 1)), "s10": ((84, -47), (143, -80)), "pell": ((12, 1), (143, 12))}
WQ = {11: ((66, 5), (143, 11)), 13: ((78, 7), (143, 13)), 143: ((0, -1), (143, 0))}


def orbit(v, A, q):
    seen, w = set(), (v[0] % q, v[1] % q)
    while True:
        p = (1, w[1] * pow(w[0], -1, q) % q) if w[0] else (0, 1)
        if p in seen:
            return frozenset(seen)
        seen.add(p)
        w = ((w[0] * A[0][0] + w[1] * A[1][0]) % q, (w[0] * A[0][1] + w[1] * A[1][1]) % q)


def coords(g):
    q = H._qcoords(m, [Fr(x) for x in symbol(g)])
    return Winv * sy.Matrix([sy.Rational(q[fj].numerator, q[fj].denominator) for fj in m["freeK"]])


def mat(l):
    return tuple(tuple(int(x) for x in row) for row in l["gamma"])


def mul(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]), (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))


fine, nrel, same_m = True, 0, True
for mo, A in MOT.items():
    fam = [i for i, l in enumerate(lines) if l["motion"] == mo]
    orb = {q: {i: orbit(lines[i]["dart"], A, q) for i in fam} for q in (11, 13)}
    same_m &= all(len(orb[11][i]) == lines[i]["m11"] and len(orb[13][i]) == lines[i]["m13"] for i in fam)
    for q, cols in ((13, range(26)), (11, [0, 1] + list(range(6, 26)))):
        for o in set(orb[q].values()):
            grp = [i for i in fam if orb[q][i] == o]
            nrel += 1
            fine &= all(sum(comp[i][j] for i in grp) == 0 for j in cols)
check("post hoc, in mtft's model alone: over each orbit modulo 13 the whole cuspidal class adds up to zero; over each orbit modulo 11 the 143a1, f2 and f3 components do", fine and same_m and nrel == 18, f"{nrel} relations")
mirror_ok = True
for i, l in enumerate(lines):
    (a, b), (c, d) = mat(l)
    im = coords(((a, -b), (-c, d)))
    sg = 1 if im == comp[i] else (-1 if im == -comp[i] else 0)
    # Route A's sign is on the class in Z^29; here only the cuspidal part is seen, where a sign found there must show
    mirror_ok &= (sg == l["mirror"]) if l["mirror"] else True
check("post hoc: the mirror image of each cat line of length 7 has the same cuspidal class, of each cat line of length 35, each Sphereland 10 line of length 30 and each Pell line the opposite one", mirror_ok and all(l["mirror"] == (1 if (l["motion"], l["ell"]) == ("cat", 7) else -1) for l in lines if (l["motion"], l["ell"]) in (("cat", 7), ("cat", 35), ("s10", 30)) or l["motion"] == "pell"))
al_ok, nimg = True, 0
for i, mo, ell, Qd, hits in R["post"]["al_images"]:
    Wm = WQ[Qd]
    h = mul(mul(Wm, mat(lines[i])), ((Wm[1][1], -Wm[0][1]), (-Wm[1][0], Wm[0][0])))
    im = coords(tuple(tuple(x // Qd for x in row) for row in h))
    for s_, j in hits:
        nimg += 1
        al_ok &= im == s_ * comp[j]
check("post hoc: where Route A found the image of a line under W_11, W_13 or W_143 to be + or - another of the 28, mtft's model gives the same line and sign", al_ok and nimg >= 20, f"{nimg} images")
print(f"\n{sum(res)} of {len(res)} checks passed")
sys.exit(0 if all(res) else 1)
