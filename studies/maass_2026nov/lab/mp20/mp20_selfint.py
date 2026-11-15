#!/usr/bin/env python3
"""MP20, a hidden self-interaction?  Plan: ../MP20_Self_Interaction_Plan.md (saved at 07:59 New York time on 14 November 2026,
before this file existed).

The participation ratio I = vol * int u^4 / (int u^2)^2 of a Maass form over the Ford domain F* of Gamma_0*(143) (area 14 pi):
to first order a cubic self-interaction -Lap u + g |u|^2 u = lambda u shifts lambda by kappa * I.

F* = { |x| <= 1/2, outside every isometric circle }: the circles have centre -Q d/(N c) and radius sqrt(Q)/(N c), Q | 143,
c >= 1, gcd(Q d, (N/Q) c) = 1.  For each y the allowed x are intervals; Gauss-Legendre on each, a midpoint rule in log y.
A form is evaluated by its own expansion, u = sum c_n sqrt(y) K_iR(2 pi f_n y) trig(2 pi f_n x), with Hejhal's coefficients
(f_n = n, or n - 1/2 for the eta-twisted forms).  Route B is a Monte Carlo over the same domain, with its own evaluation loop.

Usage: python3 -I mp20_selfint.py run [Ny]      -> mp20_numbers.json   (Ny log-spaced heights, default 1600; the first run)
       python3 -I mp20_selfint.py exact [n]     -> mp20_numbers.json   (panels between the breakpoints of the domain, n nodes each)

The first run (midpoint in log y) missed the second clause of V13: the orthogonality came out 1.4e-4 against a registered 1e-5,
because the row integral has a square-root singularity at the top of every circle and a kink at every vertex.  `exact` puts
panel ends at those heights and removes the square root with y = y_b - (y_b - y_a) tau^2."""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
spec = importlib.util.spec_from_file_location("mp18_twisted", os.path.join(HERE, "..", "mp18", "mp18_twisted.py"))
tw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tw)
m = tw.m
N = 143
VOL = 14 * math.pi
YMIN, YMAX = 0.0100, 4.0

# (name, twisted, parity, eps, R, lattice partner)
FORMS = [("0-+: odd (-,-)", False, 1, (-1, -1), 0.3738518157),
         ("0++*: even (-,-)", False, 0, (-1, -1), 0.4449803610),
         ("0-+*: odd (+,+)", False, 1, (1, 1), 0.5555514070),
         ("twisted even (+,+)", True, 0, (1, 1), 0.4800862975),
         ("twisted odd (-,-)", True, 1, (-1, -1), 0.4627415219),
         ("twisted even (+,-)", True, 0, (1, -1), 0.5851199462),
         ("twisted odd (+,-)", True, 1, (1, -1), 0.4988429343),
         ("check: even (-,-), second form", False, 0, (-1, -1), 1.0820932551)]


def circles(rmin=YMIN):
    out = []
    for Q in (1, 11, 13, 143):
        Nq, c = N // Q, 1
        while math.sqrt(Q) / (N * c) > rmin:
            r = math.sqrt(Q) / (N * c)
            for d in range(-N * c, N * c + 1):
                if math.gcd(Q * d, Nq * c) != 1:
                    continue
                xc = -Q * d / (N * c)
                if -r < xc < 0.5 + r:
                    out.append((xc, r))
            c += 1
    return sorted(set(out))


def intervals(y, circ):
    """The x in [0, 1/2] outside every circle at height y, as a list of intervals."""
    cut = sorted((xc - math.sqrt(r * r - y * y), xc + math.sqrt(r * r - y * y)) for xc, r in circ if r > y)
    out, a = [], 0.0
    for lo, hi in cut:
        if lo > a:
            out.append((a, min(lo, 0.5)))
        a = max(a, hi)
        if a >= 0.5:
            break
    if a < 0.5:
        out.append((a, 0.5))
    return [(p, q) for p, q in out if q - p > 1e-13]


_GL = {}


def gl(n):
    if n not in _GL:
        _GL[n] = np.polynomial.legendre.leggauss(n)
    return _GL[n]


def breakpoints(circ):
    """Heights where the set of allowed x changes: tops of circles and vertices (pairwise intersections), kept when the point
    lies in the closed half-strip and strictly inside no other circle."""
    C = np.array(circ)

    def free(x, y, skip):
        if not (-1e-12 <= x <= 0.5 + 1e-12):
            return False
        for k, (xc, r) in enumerate(circ):
            if k not in skip and (x - xc) ** 2 + y * y < r * r * (1 - 1e-10):
                return False
        return True
    tops = [r for k, (xc, r) in enumerate(circ) if free(xc, r, (k,))]
    verts = []
    for i in range(len(circ)):
        for j in range(i + 1, len(circ)):
            (x1, r1), (x2, r2) = circ[i], circ[j]
            d = x2 - x1
            if d < 1e-14 or d >= r1 + r2 or d <= abs(r1 - r2):
                continue
            x = x1 + (d * d + r1 * r1 - r2 * r2) / (2 * d)
            y2 = r1 * r1 - (x - x1) ** 2
            if y2 > 0 and free(x, math.sqrt(y2), (i, j)):
                verts.append((x, math.sqrt(y2)))
    return sorted(tops), sorted(verts, key=lambda v: v[1])


def exact_rows(circ, n=16, n_top=40, panels_top=12):
    """(y, weight in dy / y^2 for both halves, intervals) with panel ends at every breakpoint."""
    tops, verts = breakpoints(circ)
    ys = sorted(set([round(t, 14) for t in tops] + [round(v[1], 14) for v in verts]))
    bp = [ys[0]]
    for v in ys[1:]:
        if v - bp[-1] > 1e-11:
            bp.append(v)
    t0, w0 = np.polynomial.legendre.leggauss(n)
    tau, wt = 0.5 * (t0 + 1), 0.5 * w0
    rows = []
    for ya, yb in zip(bp[:-1], bp[1:]):
        for t, w in zip(tau, wt):
            y = yb - (yb - ya) * t * t
            rows.append((y, 2.0 * 2 * (yb - ya) * t * w / (y * y)))
    t1, w1 = np.polynomial.legendre.leggauss(n_top)
    edges = np.linspace(math.log(bp[-1]), math.log(YMAX), panels_top + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        for t, w in zip(t1, w1):
            y = math.exp(0.5 * (b - a) * t + 0.5 * (a + b))
            rows.append((y, 2.0 * 0.5 * (b - a) * w / y))
    return rows, {"tops": len(tops), "vertices": len(verts), "lowest vertex": verts[0], "highest top": tops[-1], "panels": len(bp) - 1,
                  "vertex list": verts}


def coefficient_vectors():
    Y0, M, hcs = tw.setup(0.96, 0.88)
    out = []
    for name, twisted, par, eps, R in FORMS:
        fr = np.arange(1, M + 1) - (0.5 if twisted else 0.0)
        c1, c2 = tw.coeffs(R, hcs, fr, par, eps, twist=twisted)
        out.append({"name": name, "twisted": twisted, "parity": par, "eps": eps, "R": R, "f": fr, "c": c1,
                    "agree": float(np.max(np.abs(c1[1:20] - c2[1:20])))})
        print(f"  {name}: R = {R}, horocycles agree to {out[-1]['agree']:.1e}", flush=True)
    return out, M


def run(Ny, rows=None, forms_M=None, xfac=1.0):
    circ = circles() if rows is None else circles(0.004)
    forms, M = forms_M if forms_M else coefficient_vectors()
    nf = len(forms)
    s_edges = np.linspace(math.log(YMIN), math.log(YMAX), max(Ny, 1) + 1)
    ds = s_edges[1] - s_edges[0]
    if rows is None:
        rows = [(math.exp(0.5 * (s_edges[k] + s_edges[k + 1])), None) for k in range(Ny)]
    Ny = len(rows)
    area = 0.0
    N2, N4 = np.zeros(nf), np.zeros(nf)
    cross = {(1, 7): 0.0}                                        # the two even (-,-) forms: orthogonality
    groups = {}
    for i, f in enumerate(forms):
        groups.setdefault((f["twisted"], f["parity"]), []).append(i)
    for k in range(Ny):
        y, wrow = rows[k]
        iv = intervals(y, circ)
        if not iv:
            continue
        Meff = min(M, int(30.0 / (2 * math.pi * y)) + 2)
        xs, ws = [], []
        for a, b in iv:
            x0, w0 = gl(int(xfac * (2 * math.pi * Meff * (b - a) + 30)))
            xs.append(0.5 * (b - a) * x0 + 0.5 * (a + b))
            ws.append(0.5 * (b - a) * w0)
        x, w = np.concatenate(xs), np.concatenate(ws)
        wy = 2.0 * ds / y if wrow is None else wrow              # both halves; dx dy / y^2 with dy = y ds
        area += wy * w.sum()
        vals = {}
        for (twisted, par), idx in groups.items():
            fr = forms[idx[0]]["f"][:Meff]
            Tm = (np.cos if par == 0 else np.sin)(2 * np.pi * np.outer(fr, x))
            for i in idx:
                amp = forms[i]["c"][:Meff] * m.Kir(forms[i]["R"], 2 * np.pi * fr * y)
                u = math.sqrt(y) * (amp @ Tm)
                vals[i] = u
                N2[i] += wy * float(np.sum(w * u * u))
                N4[i] += wy * float(np.sum(w * u ** 4))
        cross[(1, 7)] += wy * float(np.sum(w * vals[1] * vals[7]))
        if k % 400 == 0:
            print(f"  y = {y:.5f} ({k}/{Ny}), {len(iv)} intervals, {x.size} nodes, {time.time() - T0:.0f} s", flush=True)
    area += 1.0 / YMAX                                           # the strip above YMAX, both halves
    I = VOL * N4 / N2 ** 2
    return {"Ny": Ny, "circles": len(circ), "area": area, "area / 14 pi": area / VOL, "N2": N2.tolist(), "N4": N4.tolist(), "I": I.tolist(),
            "orthogonality of the two even (-,-) forms": cross[(1, 7)] / math.sqrt(N2[1] * N2[7]),
            "names": [f["name"] for f in forms], "agree": [f["agree"] for f in forms]}, forms, circ


def monte_carlo(forms, circ, which, n_acc=30000, seed=20):
    """Route B: points uniform for dx dy / y^2 in the strip above YMIN, kept if outside every circle; u from a plain loop."""
    rng = np.random.default_rng(seed)
    C = np.array(circ)
    pts = []
    while len(pts) < n_acc:
        yv = YMIN / rng.random(20000)
        xv = rng.random(20000) - 0.5
        keep = np.ones(20000, dtype=bool)
        for xc, r in C:
            keep &= (np.abs(xv) - xc) ** 2 + yv ** 2 >= r * r
        pts += list(zip(np.abs(xv[keep]), yv[keep]))
    pts = np.array(pts[:n_acc])
    pts = pts[np.argsort(-pts[:, 1])]
    out = {}
    for i in which:
        f = forms[i]
        u = np.zeros(n_acc)
        for a in range(0, n_acc, 500):
            xb, yb = pts[a:a + 500, 0], pts[a:a + 500, 1]
            Meff = min(len(f["c"]), int(30.0 / (2 * math.pi * yb.min())) + 2)
            fr = f["f"][:Meff]
            K = m.Kir(f["R"], 2 * np.pi * np.outer(yb, fr))
            trig = (np.cos if f["parity"] == 0 else np.sin)(2 * np.pi * np.outer(xb, fr))
            u[a:a + 500] = np.sqrt(yb) * np.sum(K * trig * f["c"][:Meff][None, :], axis=1)
        m2, m4 = np.mean(u ** 2), np.mean(u ** 4)
        # delta-method error of m4 / m2^2
        g = u ** 4 / m2 ** 2 - 2 * m4 * u ** 2 / m2 ** 3
        out[f["name"]] = {"I": float(m4 / m2 ** 2), "standard error": float(np.std(g) / math.sqrt(n_acc))}
    return out


def analysis(res):
    I = dict(zip(res["names"], res["I"]))
    need = json.load(open(os.path.join(HERE, "..", "needed_shifts.json")))
    lat = {"0++": (3.405, 0.021), "0++*": (5.855, 0.041), "0-+": (5.276, 0.045), "0-+*": (7.29, 0.13)}
    ks = ["0-+", "0++*", "0-+*"]
    names = ["0-+: odd (-,-)", "0++*: even (-,-)", "0-+*: odd (+,+)"]
    m0, s0 = lat["0++"]
    r = [lat[k][0] / m0 for k in ks]
    cov = np.array([[r[i] * r[j] * ((s0 / m0) ** 2 + ((lat[ks[i]][1] / lat[ks[i]][0]) ** 2 if i == j else 0)) for j in range(3)] for i in range(3)]) / 16
    d = np.array([need[k]["needed shift in lambda"] for k in ks])
    Iv = np.array([I[n] for n in names])
    Ci = np.linalg.inv(cov)
    kap = float(Iv @ Ci @ d / (Iv @ Ci @ Iv))
    chi = float((d - kap * Iv) @ Ci @ (d - kap * Iv))
    p1 = {"kappa_plus": kap, "chi2": chi, "dof": 2, "p": float(chi2.sf(chi, 2)), "pass": bool(chi2.sf(chi, 2) > 0.05),
          "I": dict(zip(names, Iv.tolist())), "needed shifts": d.tolist(), "fitted shifts": (kap * Iv).tolist(),
          "chi2 with no interaction": float(d @ Ci @ d), "needed I ratio odd(-,-)/even(-,-)": float(d[0] / d[1]), "actual I ratio": float(Iv[0] / Iv[1])}
    # P2
    lam = lambda R: 0.25 + R * R
    A_even, A_odd = (9.32 / 3.405, need["A1+-"]["error"]), (9.78 / 3.405, need["A1--"]["error"])
    Ie, Io, Ie2, Io2 = I["twisted even (+,+)"], I["twisted odd (-,-)"], I["twisted even (+,-)"], I["twisted odd (+,-)"]
    km = need["A1+-"]["needed shift in lambda"] / Ie
    skm = need["A1+-"]["its error"] / Ie
    pred = {"twisted odd (-,-)": 4 * (lam(0.4627415219) + km * Io), "twisted odd (+,-)": 4 * (lam(0.4988429343) + km * Io2),
            "twisted even (+,-)": 4 * (lam(0.5851199462) + km * Ie2)}
    lowest_odd = min(pred["twisted odd (-,-)"], pred["twisted odd (+,-)"])
    Ilow = Io if pred["twisted odd (-,-)"] <= pred["twisted odd (+,-)"] else Io2
    sig_low = math.hypot(A_odd[1], 4 * skm * Ilow)
    n_low = (lowest_odd - A_odd[0]) / sig_low
    other_odd = max(pred["twisted odd (-,-)"], pred["twisted odd (+,-)"])
    n_even2 = (pred["twisted even (+,-)"] - A_even[0]) / math.hypot(A_even[1], 4 * skm * Ie2)
    p2 = {"kappa_minus": km, "its error": skm, "shifted ratios": pred, "lowest odd against A1--": {"value": lowest_odd, "sigma": n_low},
          "second even against A1+- (may not lie 2 sigma below)": {"value": pred["twisted even (+,-)"], "sigma": n_even2},
          "second odd": other_odd, "pass": bool(abs(n_low) <= 2 and n_even2 >= -2)}
    f1 = {"I": {n: I[n] for n in res["names"][:7]}, "pass": bool(all(1.5 <= I[n] <= 4.5 for n in res["names"][:7]))}
    return {"P1": p1, "P2": p2, "F1": f1}


def main_exact(n):
    circ = circles(0.004)                                        # every circle that could reach the domain, and more
    t_small, v_small = breakpoints(circles())
    t_all, v_all = breakpoints(circ)
    assert len(t_small) == len(t_all) and len(v_small) == len(v_all) and abs(v_small[0][1] - v_all[0][1]) < 1e-15, "smaller circles reach the domain"
    fm = coefficient_vectors()
    rows, info = exact_rows(circ, n)
    print(f"  {info['tops']} tops, {info['vertices']} vertices, {info['panels']} panels, {len(rows)} heights; lowest vertex {info['lowest vertex']}", flush=True)
    res, forms, _ = run(0, rows, fm)
    res["domain"] = info
    print(f"area / 14 pi - 1 = {res['area / 14 pi'] - 1:.2e}; orthogonality {res['orthogonality of the two even (-,-) forms']:.2e}; {time.time() - T0:.0f} s")
    for nm, v in zip(res["names"], res["I"]):
        print(f"  I = {v:.7f}   {nm}")
    out = {"quadrature": res, "analysis": analysis(res)}
    rows2, _ = exact_rows(circ, n + 8, 56, 12)
    finer, _, _ = run(0, rows2, fm, xfac=1.5)
    out["finer quadrature (n + 8 nodes a panel, half as many again in x)"] = {"I": finer["I"], "area / 14 pi - 1": finer["area / 14 pi"] - 1,
                                                                             "orthogonality": finer["orthogonality of the two even (-,-) forms"]}
    print("finer: max change of I", float(np.max(np.abs(np.array(finer["I"]) - np.array(res["I"])))), "orthogonality", finer["orthogonality of the two even (-,-) forms"])
    out["Route B (Monte Carlo)"] = monte_carlo(forms, circ, (0, 1, 2, 3, 4, 5, 6), n_acc=60000)
    print("Route B:", {k: (round(v["I"], 4), round(v["standard error"], 4)) for k, v in out["Route B (Monte Carlo)"].items()})
    out["first run (midpoint in log y, 1600 heights)"] = json.load(open(os.path.join(HERE, "mp20_numbers_first_run.json")))["quadrature"]
    out["seconds"] = time.time() - T0
    json.dump(out, open(os.path.join(HERE, "mp20_numbers.json"), "w"), indent=1)
    a = out["analysis"]
    print("P1:", {k: a["P1"][k] for k in ("kappa_plus", "chi2", "p", "pass", "needed I ratio odd(-,-)/even(-,-)", "actual I ratio")})
    print("P2:", {k: a["P2"][k] for k in ("kappa_minus", "shifted ratios", "lowest odd against A1--", "second even against A1+- (may not lie 2 sigma below)", "pass")})
    print("F1:", a["F1"]["pass"], f"({time.time() - T0:.0f} s)")


if __name__ == "__main__":
    if sys.argv[1] == "exact":
        main_exact(int(sys.argv[2]) if len(sys.argv) > 2 else 16)
        sys.exit()
    Ny = int(sys.argv[2]) if len(sys.argv) > 2 else 1600
    res, forms, circ = run(Ny)
    print(f"area / 14 pi = {res['area / 14 pi']:.7f}; orthogonality {res['orthogonality of the two even (-,-) forms']:.2e}; {time.time() - T0:.0f} s")
    for n, v in zip(res["names"], res["I"]):
        print(f"  I = {v:.5f}   {n}")
    out = {"quadrature": res, "analysis": analysis(res)}
    if Ny >= 1600:
        half, _, _ = run(Ny // 2)
        out["quadrature at half the heights"] = {"I": half["I"], "area / 14 pi": half["area / 14 pi"]}
        out["Route B (Monte Carlo)"] = monte_carlo(forms, circ, (0, 1, 3, 4))
        print("Route B:", {k: (round(v["I"], 4), round(v["standard error"], 4)) for k, v in out["Route B (Monte Carlo)"].items()})
    out["seconds"] = time.time() - T0
    json.dump(out, open(os.path.join(HERE, "mp20_numbers.json"), "w"), indent=1)
    a = out["analysis"]
    print("P1:", {k: a["P1"][k] for k in ("kappa_plus", "chi2", "p", "pass", "needed I ratio odd(-,-)/even(-,-)", "actual I ratio")})
    print("P2:", {k: a["P2"][k] for k in ("kappa_minus", "shifted ratios", "lowest odd against A1--", "second even against A1+- (may not lie 2 sigma below)", "pass")})
    print("F1:", a["F1"]["pass"], f"({time.time() - T0:.0f} s)")
