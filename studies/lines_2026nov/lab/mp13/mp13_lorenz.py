#!/usr/bin/env python3
"""
mp13_lorenz.py  --  Route C for MP13 (test T-G9): level 1, through Ghys's Lorenz template.

Modular knots and links are Lorenz knots and links (Ghys): a primitive hyperbolic A of
SL(2, Z) with positive trace is conjugate to a positive word in U = [[1, 1], [0, 1]] and
V = [[1, 0], [1, 1]], and its periodic orbit is the Lorenz orbit with the word read in
x = U, y = V.  On the template every crossing lies where the left strip passes the
right one, and two orbit strands cross there exactly when the strand that starts on the
left (first letter x) lands to the right of the one that starts on the right; the
linking number of two orbits is half the number of their mutual crossings.

The cat orbit bounds in PSL(2, Z)\\PSL(2, R) (psi = 0), and M = Gamma_0(143)\\PSL(2, R) covers
it with 168 sheets, so (MP13, P4) the totals of crossings on the curve found by Route A
must be 168 times the sums of these linking numbers.  Nothing here uses geometry.

Run:  python3 mp13_lorenz.py mp13_numbers.json
"""
import sys, json
from fractions import Fraction as Fr

R13 = json.load(open(sys.argv[1]))["numbers"]
res = []


def check(name, cond, detail=""):
    res.append(bool(cond))
    print(("ok   " if cond else "FAIL ") + name + (f"   [{detail}]" if detail else ""), flush=True)


def mm(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))


inv = lambda A: ((A[1][1], -A[0][1]), (-A[1][0], A[0][0]))
U, V = ((1, 1), (0, 1)), ((1, 0), (1, 1))
GENS = [((0, -1), (1, 0)), ((1, 1), (0, 1)), ((1, -1), (0, 1))]


def positive_word(A):
    """A conjugate of A with entries >= 0 (breadth-first over conjugation by S, T, T^-1), peeled into U and V."""
    if A[0][0] + A[1][1] < 0:
        A = tuple(tuple(-v for v in r) for r in A)
    frontier, seen = [A], {A}
    for _ in range(40):
        for B in frontier:
            if min(B[0] + B[1]) >= 0:
                (a, b), (c, d) = B
                w = ""
                while (a, b, c, d) != (1, 0, 0, 1):
                    if a >= c and b >= d:
                        a, b, w = a - c, b - d, w + "x"
                    else:
                        c, d, w = c - a, d - b, w + "y"
                    assert min(a, b, c, d) >= 0
                P = ((1, 0), (0, 1))
                for ch in w:
                    P = mm(P, U if ch == "x" else V)
                assert P == B
                return w
        nxt = []
        for B in frontier:
            for h in GENS:
                C = mm(mm(h, B), inv(h))
                if C not in seen:
                    seen.add(C)
                    nxt.append(C)
        frontier = nxt
    raise RuntimeError("no positive conjugate found")


def canon(w):                                          # the cyclic word, by its least rotation
    return min(w[i:] + w[:i] for i in range(len(w)))


def less(u, v):                                        # u^infinity < v^infinity, lexicographically with x < y
    n = 2 * (len(u) + len(v))
    return (u * (n // len(u) + 1))[:n] < (v * (n // len(v) + 1))[:n]


def lorenz_lk(a, b):
    """Half the crossings between the strands of the Lorenz orbits a and b."""
    sa = [a[i:] + a[:i] for i in range(len(a))]
    sb = [b[i:] + b[:i] for i in range(len(b))]
    shift = lambda s: s[1:] + s[0]
    cr = 0
    for P, Q in ((sa, sb), (sb, sa)):
        for p in P:
            if p[0] != "x":
                continue
            for q in Q:
                if q[0] == "y" and less(shift(q), shift(p)):
                    cr += 1
    assert cr % 2 == 0
    return cr // 2


MOT = {"cat": ((2, 1), (1, 1)), "s10": ((84, -47), (143, -80)), "pell": ((12, 1), (143, 12))}
W = {k: canon(positive_word(A)) for k, A in MOT.items()}
Wi = {k: canon(positive_word(inv(A))) for k, A in MOT.items()}
psi = {"cat": 0, "s10": -1, "pell": -21}
check("the positive words: cat xy, Sphereland 10 xyy, Pell x y^22, and of the inverses xy, xxy, x^22 y; #x - #y is psi (Ghys)",
      (W["cat"], W["s10"], W["pell"]) == ("xy", "xyy", "x" + "y" * 22) and (Wi["cat"], Wi["s10"], Wi["pell"]) == ("xy", "xxy", "x" * 22 + "y")
      and all(W[k].count("x") - W[k].count("y") == psi[k] for k in MOT), f"{W}; inverses {Wi}")
# the template on small cases: the two fixed orbits are unlinked, and the order-2 orbit links each order-3 one once
check("the template on small orbits: Lk(x, y) = 0, Lk(xy, xxy) = Lk(xy, xyy) = 1, and Lk is symmetric", lorenz_lk("x", "y") == 0 and lorenz_lk("xy", "xxy") == 1 and lorenz_lk("xy", "xyy") == 1 and lorenz_lk("xyy", "xy") == 1)
y22, x22 = "x" + "y" * 22, "x" * 22 + "y"
vals = {"Lk(xy, xyy)": lorenz_lk("xy", "xyy"), "Lk(xy, xxy)": lorenz_lk("xy", "xxy"), "Lk(xy, xy^22)": lorenz_lk("xy", y22), "Lk(xy, x^22y)": lorenz_lk("xy", x22),
        "Lk(xyy, xy^22)": lorenz_lk("xyy", y22), "Lk(xxy, x^22y)": lorenz_lk("xxy", x22), "Lk(xyy, x^22y)": lorenz_lk("xyy", x22), "Lk(xxy, xy^22)": lorenz_lk("xxy", y22)}
want = [1, 1, 1, 1, 2, 2, 1, 1]
check("P4 the pre-registered template values 1, 1, 1, 1, 2, 2, 1, 1", list(vals.values()) == want, str(vals))
pred = {"cat_s10": 168 * (vals["Lk(xy, xyy)"] + vals["Lk(xy, xxy)"]), "cat_pell": 168 * (vals["Lk(xy, xy^22)"] + vals["Lk(xy, x^22y)"]),
        "s10_pell": 168 * (vals["Lk(xyy, xy^22)"] + vals["Lk(xyy, x^22y)"] + vals["Lk(xxy, xy^22)"] + vals["Lk(xxy, x^22y)"])}
found = {k: R13["totals"][k] for k in pred}
check("P4 the transfer: Route A's totals of crossings on the curve are 168 times the template's sums", pred == found, f"template x 168: {pred}; Route A: {found}")
lk = R13["totals"]["Lk_cat_B"]
check("P4 in M: Lk(all cat lifts, all lifts of B) = -168 for B = Sphereland 10 and Pell, and its size is 168 Lk(xy, word of B)", lk == {"s10": -168, "pell": -168} and 168 * vals["Lk(xy, xyy)"] == 168 and 168 * vals["Lk(xy, xy^22)"] == 168, str(lk))
print(f"\n{sum(res)} of {len(res)} checks passed")
sys.exit(0 if all(res) else 1)
