"""MP23, Q5 (the plan's caution): is the frozen frame of the holomorphic class unique?  The fit of Q2 at degree <= 2 had full rank.
Here the same fit at degree <= 4 (15 monomials, 120 real unknowns): its rank, and the solution."""
import json, sys
import numpy as np
sys.argv = [sys.argv[0], "none"]
sys.path.insert(0, "/home/claude/mp/mp23")
import mp23_tg1 as m
MON = m.monos(4)
c, res, rank, sv, M, y = m.fit(MON, m.rand_points(40, 23), -1)
c = c.reshape(len(MON), 4)
big = {str(mono): [[float(v.real), float(v.imag)] for v in row] for mono, row in zip(MON, c) if np.abs(row).max() > 1e-9}
out = {"degree": 4, "real unknowns": 2 * c.size, "rank": rank, "residual / |W|": res, "singular values, largest and three smallest": [float(sv[0])] + [float(v) for v in sv[-3:]],
       "monomials with a coefficient above 1e-9": big, "largest coefficient of degree 2 to 4": float(np.abs(c[3:]).max())}
print(json.dumps(out, indent=1))
d = json.load(open("/home/claude/mp/mp23/mp23_numbers.json"))
d["Q5: the fit at degree <= 4"] = out
d.get("checks", {}).pop("run of sections none", None)
json.dump(d, open("/home/claude/mp/mp23/mp23_numbers.json", "w"), indent=1)
