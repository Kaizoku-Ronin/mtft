"""mtft.surface.smflux — multi-stack flux models on X0(143): search, anomaly ledger, Higgs sectors, purity
(SM-01, v0.30.0).  Model class (declared): gauge group U(3)_c x U(2)_L x U(1)^k, one line-bundle flux L_x of
degree m_x per stack, fermions in all bifundamentals (x, y-bar) with internal bundle S0 (x) L_x (x) L_y^{-1};
net left-handed multiplicity = m_x − m_y (index theorem).  Hypercharge Y = sum y_x Q_x, y_x in Z/6, defined up
to an overall shift (fix y_c = 1/6).

Results (EXACT, integer arithmetic): no four-stack model gives three chiral SM families (telescoping sum rule);
five-stack models exist — 30 with |m| <= 6 and y in [−1, 1] — all with m_L = −3 and >= 6 vector-like
doublet pairs.  Right-handed neutrinos are allowed (net count reported).  Realisation on X0(143): odd degrees
from W13-fixed CM points (W11 must break); purity of each chiral sector with S0 from the sign structure of the
unit u (two fixed points of each sign) and its cusp orders.
"""
from __future__ import annotations

import itertools
from fractions import Fraction as Fr
from typing import Dict, List

import numpy as np

SM = {("3", "2", Fr(1, 6)): 3, ("3bar", "1", Fr(-2, 3)): 3, ("3bar", "1", Fr(1, 3)): 3, ("1", "2", Fr(-1, 2)): 3, ("1", "1", Fr(1)): 3}
NU = ("1", "1", Fr(0))


def conj(rep):
    col, wk, Y = rep; return ({"3": "3bar", "3bar": "3", "1": "1"}[col], wk, -Y)


def canon(rep):
    c = conj(rep); return (rep, 1) if str(rep) <= str(c) else (c, -1)


def pair_rep(x, y):
    col = "3" if x == "c" else ("3bar" if y == "c" else "1"); wk = "2" if "L" in (x, y) else "1"; return col, wk


def search(stacks: List[str], R: int, ylist: Dict[str, List[Fr]]) -> List[Dict]:
    pairs = [(x, y) for i, x in enumerate(stacks) for y in stacks[i + 1:]]
    ms = np.array(list(itertools.product(range(-R, R + 1), repeat=len(stacks) - 1)))
    M = np.column_stack([np.zeros(len(ms), int), ms]); idx = {s: i for i, s in enumerate(stacks)}
    target = {}
    for rep, k in SM.items():
        crep, sgn = canon(rep); target[crep] = target.get(crep, 0) + sgn * k
    nuc, _ = canon(NU); sols = []
    for yv in itertools.product(*[ylist[s] for s in stacks]):
        y = dict(zip(stacks, yv)); net = {}
        for x, z in pairs:
            col, wk = pair_rep(x, z); crep, sgn = canon((col, wk, y[x] - y[z])); net[crep] = net.get(crep, 0) + sgn * (M[:, idx[x]] - M[:, idx[z]])
        if any(c not in net for c in target): continue
        good = np.ones(len(M), bool)
        for crep, v in net.items():
            if crep == nuc: continue
            good &= v == target.get(crep, 0)
        for i in np.where(good)[0]:
            tot = {}
            for x, z in pairs:
                col, wk = pair_rep(x, z); crep, _ = canon((col, wk, y[x] - y[z])); tot[crep] = tot.get(crep, 0) + abs(int(M[i, idx[x]] - M[i, idx[z]]))
            nvl = sum((tot[c] - abs(int(net[c][i]))) // 2 for c in tot)
            sols.append({"degrees": dict(zip(stacks, M[i].tolist())), "y": {k: v for k, v in y.items()}, "vector_like_pairs": nvl, "net_nu_c": int(net[nuc][i]) if nuc in net else 0})
    return sols


def rep_table(model: Dict) -> List[Dict]:
    stacks = list(model["degrees"]); m = model["degrees"]; y = model["y"]; rows = []
    for i, x in enumerate(stacks):
        for z in stacks[i + 1:]:
            col, wk = pair_rep(x, z); n = m[x] - m[z]; rep = (col, wk, y[x] - y[z])
            if n < 0: rep, n = conj(rep), -n
            if n: rows.append({"sector": f"({x},{z}bar)" if m[x] - m[z] > 0 else f"({x}bar,{z})", "rep": rep, "multiplicity": n, "bundle_degree": m[x] - m[z], "charges": {s: (1 if s == x else -1 if s == z else 0) * (1 if m[x] - m[z] > 0 else -1) for s in stacks}})
    return rows


def anomaly_ledger(model: Dict) -> Dict:
    """Mixed anomalies SU(3)^2 U(1)_x, SU(2)^2 U(1)_x, grav U(1)_x for each stack U(1) and for Y; anomaly-free
    combinations of the stack U(1)s (rational kernel of the mixed-anomaly matrix)."""
    rows = rep_table(model); stacks = list(model["degrees"]); y = model["y"]
    dimc = {"3": 3, "3bar": 3, "1": 1}; dimw = {"2": 2, "1": 1}
    def A(qfun):
        a3 = sum(r["multiplicity"] * qfun(r) * dimw[r["rep"][1]] * (Fr(1, 2) if r["rep"][0] != "1" else 0) for r in rows)
        a2 = sum(r["multiplicity"] * qfun(r) * dimc[r["rep"][0]] * (Fr(1, 2) if r["rep"][1] == "2" else 0) for r in rows)
        ag = sum(r["multiplicity"] * qfun(r) * dimc[r["rep"][0]] * dimw[r["rep"][1]] for r in rows)
        return {"SU3^2": a3, "SU2^2": a2, "grav": ag}
    led = {x: A(lambda r, x=x: r["charges"][x]) for x in stacks}
    led["Y"] = A(lambda r: r["rep"][2])
    # anomaly-free rational combinations of stack U(1)s (mixed with SU(3)^2 and SU(2)^2)
    from sympy import Matrix, Rational
    Mat = Matrix([[Rational(led[x]["SU3^2"]) for x in stacks], [Rational(led[x]["SU2^2"]) for x in stacks], [Rational(led[x]["grav"]) for x in stacks]])
    ker = Mat.nullspace()
    return {"ledger": {k: {kk: str(vv) for kk, vv in v.items()} for k, v in led.items()}, "anomaly_free_combinations": [[str(c) for c in v] for v in ker],
            "Y_in_anomaly_free_span": bool(Matrix([[Rational(y[s]) for s in stacks]]).T in [Matrix(v) for v in ker]) or Matrix.hstack(*ker).rank() == Matrix.hstack(*ker, Matrix([Rational(y[s]) for s in stacks])).rank() if ker else False}


def purity_with_S0(degree_difference: int, kind: str) -> Dict:
    """Brill–Noether purity of S0 (x) O(D) for the standard realisations on X0(143):
    kind 'cm': D = ± (P1+P2+P3) (three W13-fixed points), degree ±3; 'cm2': D = ± 2(P1+P2+P3), degree ±6;
    'trivial': D = 0 (h^0 = h^1 = 2, two vector-like pairs); 'generic0': generic degree-0 twist (0, 0).
    Uses: H^0(S0) = span{1, u}, u(P) = ±i/sqrt 13 with two of each sign among the four fixed points."""
    d = degree_difference
    if kind == "cm" and abs(d) == 3:
        return {"h0": 3 if d > 0 else 0, "h1": 0 if d > 0 else 3, "pure": True, "reason": "no combination a+bu vanishes at three fixed points (both signs present)"}
    if kind == "cm2" and abs(d) == 6:
        return {"h0": 6 if d > 0 else 0, "h1": 0 if d > 0 else 6, "pure": True, "reason": "a+bu cannot vanish to order 2 at three points"}
    if kind == "trivial" and d == 0: return {"h0": 2, "h1": 2, "pure": False, "reason": "S0 itself: two vector-like pairs (constant and u)"}
    if kind == "generic0" and d == 0: return {"h0": 0, "h1": 0, "pure": True, "reason": "generic point of the Jacobian is off the theta divisor of S0"}
    return {"h0": None, "h1": None, "pure": None, "reason": "not covered"}


# ------------------------------------------------ SM-14: hypercharge normalisation of a stack model
def hypercharge_normalisation(N, y):
    """For stacks U(N_x) with one 6D gauge coupling g and hypercharge Y = sum_x y_x Q_x (Q_x the U(1) charge, +1 on the
    fundamental), 1/g_Y^2 = (2/g^2) sum_x N_x y_x^2 =: k_Y / g^2.  Returns k_Y and sin^2 theta_W = 1/(1 + k_Y) at the
    compactification scale, with g_3 = g_2 = sqrt(k_Y) g_Y.  M1/M2: N = (3, 2, 1, 1, 1), y = (1/6, 0, -1/2, 1/2, -1/2)
    give k_Y = 5/3 and sin^2 theta_W = 3/8."""
    from fractions import Fraction
    kY = 2 * sum(Fraction(int(n)) * Fraction(yy) ** 2 for n, yy in zip(N, y))
    return {"k_Y": kY, "sin2_thetaW": 1 / (1 + kY), "g3_over_gY": kY ** Fraction(1, 2) if kY.denominator == 1 and int(kY) ** 0.5 == int(int(kY) ** 0.5) else None}


# ------------------------------------------------ review V0311 §2/§4: exact Abelian anomaly polynomial, B-L, Majorana obstruction
def abelian_anomaly_polynomial(N, m, q):
    """Cubic U(1) anomaly of the bifundamental spectrum: sum_{a<b} (m_a - m_b) N_a N_b (q_a - q_b)^3 for a charge vector q on the
    stacks (exact rationals).  Vanishes identically on span{common phase, Y, B-L} for M1 (review V0311 §2, verified)."""
    from fractions import Fraction
    names = list(N); return sum(Fraction(m[a] - m[b]) * N[a] * N[b] * (Fraction(q[a]) - Fraction(q[b])) ** 3 for i, a in enumerate(names) for b in names[i + 1:])

def mixed_nonabelian_anomalies(N, m, q):
    """U(1)_q – SU(N_y)^2 anomalies for every non-Abelian stack y: sum_z N_z (m_y - m_z)(q_y - q_z) (exact rationals)."""
    from fractions import Fraction
    return {y: sum(N[z] * Fraction(m[y] - m[z]) * (Fraction(q[y]) - Fraction(q[z])) for z in N if z != y) for y in N if N[y] >= 2}

M1_STACKS = {"N": {"c": 3, "L": 2, "a": 1, "b": 1, "d": 1}, "m": {"c": 0, "L": -3, "a": 3, "b": 3, "d": 0},
             "Y": {"c": "1/6", "L": 0, "a": "-1/2", "b": "1/2", "d": "-1/2"}, "B-L": {"c": "1/3", "L": 0, "a": 0, "b": 0, "d": -1}}

def majorana_obstruction(BL=None):
    """nu^c = (a, d-bar) carries B-L = BL[a] - BL[d] = +1, so nu^c nu^c has charge +2: a bare Majorana mass is forbidden while
    B-L is unbroken.  A seesaw needs a B-L = -2 scalar (or B-L breaking by another mechanism) — outside the present class."""
    from fractions import Fraction
    BL = BL or M1_STACKS["B-L"]; c = Fraction(BL["a"]) - Fraction(BL["d"]); return {"nu_c_B_minus_L": c, "majorana_bilinear_charge": 2 * c, "bare_majorana_allowed": 2 * c == 0}
