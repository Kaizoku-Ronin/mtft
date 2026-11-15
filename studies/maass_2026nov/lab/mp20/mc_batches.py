"""Route B, error from the scatter of batches: 40 batches of 10,000 points for the even (-,-) form and the lowest odd form.
The ratio I = <u^4>/<u^2>^2 is taken on the pooled means, its error by a jackknife over the batches."""
import json, math, sys, time
import numpy as np
sys.path.insert(0, "/home/claude/mp/mp20")
import mp20_selfint as s
m = s.m
t0 = time.time()
circ = s.circles(0.004)
C = np.array(circ)
forms, M = s.coefficient_vectors()
q = json.load(open("mp20_numbers.json"))["quadrature"]
rng = np.random.default_rng(777)
NB, NP = 40, 10000
which = (0, 1, 3, 4)
m2 = {i: [] for i in which}; m4 = {i: [] for i in which}; top = {i: [] for i in which}
for b in range(NB):
    pts = []
    while len(pts) < NP:
        yv = s.YMIN * 0.98 / rng.random(20000)                   # the strip above 0.0098, below the lowest vertex
        xv = rng.random(20000) - 0.5
        keep = np.ones(20000, dtype=bool)
        for xc, r in C:
            keep &= (np.abs(xv) - xc) ** 2 + yv ** 2 >= r * r
        pts += list(zip(np.abs(xv[keep]), yv[keep]))
    pts = np.array(pts[:NP]); pts = pts[np.argsort(-pts[:, 1])]
    for i in which:
        f = forms[i]; u = np.zeros(NP)
        for a in range(0, NP, 500):
            xb, yb = pts[a:a + 500, 0], pts[a:a + 500, 1]
            Meff = min(len(f["c"]), int(30.0 / (2 * math.pi * yb.min())) + 2)
            fr = f["f"][:Meff]
            K = m.Kir(f["R"], 2 * np.pi * np.outer(yb, fr))
            trig = (np.cos if f["parity"] == 0 else np.sin)(2 * np.pi * np.outer(xb, fr))
            u[a:a + 500] = np.sqrt(yb) * np.sum(K * trig * f["c"][:Meff][None, :], axis=1)
        m2[i].append(np.mean(u ** 2)); m4[i].append(np.mean(u ** 4)); top[i].append(float(np.sort(u ** 4)[-10:].sum() / np.sum(u ** 4)))
    if b % 10 == 9:
        print(f"  batch {b + 1}, {time.time() - t0:.0f} s", flush=True)
out = {}
for i in which:
    a2, a4 = np.array(m2[i]), np.array(m4[i])
    I = a4.mean() / a2.mean() ** 2
    jk = np.array([np.delete(a4, k).mean() / np.delete(a2, k).mean() ** 2 for k in range(NB)])
    se = math.sqrt((NB - 1) / NB * np.sum((jk - jk.mean()) ** 2))
    per = a4 / a2 ** 2
    Iq = q["I"][i]
    out[forms[i]["name"]] = {"I": float(I), "jackknife error": se, "quadrature": Iq, "sigma": float((I - Iq) / se), "scatter of the 40 single-batch ratios": float(per.std(ddof=1)),
                             "share of sum u^4 in the 10 largest points of a batch, mean": float(np.mean(top[i]))}
    print(forms[i]["name"], {k: round(v, 5) for k, v in out[forms[i]["name"]].items()})
json.dump(out, open("mc_batches.json", "w"), indent=1)
