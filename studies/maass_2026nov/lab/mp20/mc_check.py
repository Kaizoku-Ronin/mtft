"""Route B again with other seeds and more points: is the 3-sigma gap on the even (-,-) form a fluctuation?"""
import json, math, sys, time
import numpy as np
sys.path.insert(0, "/home/claude/mp/mp20")
import mp20_selfint as s
t0 = time.time()
circ = s.circles(0.004)
forms, M = s.coefficient_vectors()
q = json.load(open("mp20_numbers.json"))["quadrature"]
out = {}
for seed in (101, 202, 303):
    r = s.monte_carlo(forms, circ, (0, 1, 2), n_acc=120000, seed=seed)
    out[seed] = r
    print(seed, {k: (round(v["I"], 4), round(v["standard error"], 4), round((v["I"] - q["I"][q["names"].index(k)]) / v["standard error"], 2)) for k, v in r.items()}, f"{time.time() - t0:.0f} s", flush=True)
json.dump(out, open("mc_check.json", "w"), indent=1)
