"""MP22, after the run (not predicted): why MP11's search lost the line at t = 0.8.  MP11 looked, from each return, for an upward
crossing of the plane through (0.6, -0.3, 0.5) across the field there, within an arc of 60.  Follow the line for 100 circuits and
list the gaps, in arc length in space, between successive upward crossings of that plane."""
import json, math, sys
import numpy as np
from scipy.integrate import solve_ivp
sys.argv = [sys.argv[0], "none"]
sys.path.insert(0, "/home/claude/mp/mp22")
import mp22_tg8 as m
tv = 0.8
u = m.field(tv, "B")
n0 = u(m.X0H)
ev = lambda s, X: (X - m.X0H) @ n0
ev.direction = 1
so = solve_ivp(lambda s, X: u(X), (0, 1500.0), m.X0H, method="DOP853", rtol=1e-11, atol=1e-13, events=ev)
s = np.concatenate([[0.0], so.t_events[0]])
gaps = np.diff(s)
d = np.linalg.norm(so.y_events[0] - m.X0H, axis=1)
mod = m.Model(tv)
out = {"arc length followed": float(so.t[-1]), "upward crossings": int(gaps.size), "longest gap": float(gaps.max()), "gaps above 60": int((gaps > 60).sum()),
       "typical gap (median)": float(np.median(gaps)), "first ten gaps": [float(v) for v in gaps[:14]], "distance of the crossings from the point, least and most": [float(d.min()), float(d.max())],
       "drift of I_2 / range": float(np.ptp(mod.I2_space(so.y.T, "B")) / mod.I2_range)}
print(json.dumps(out, indent=1))
allnum = json.load(open("/home/claude/mp/mp22/mp22_numbers.json"))
allnum["after the run: MP11's plane at t = 0.8"] = out
json.dump(allnum, open("/home/claude/mp/mp22/mp22_numbers.json", "w"), indent=1)
