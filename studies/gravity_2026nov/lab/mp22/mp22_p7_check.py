"""MP22, P7, after its first run printed MISS on the count of crossings (66, not 65): where are the crossings of y = 0?"""
import json, math, sys
import numpy as np
from scipy.integrate import solve_ivp
sys.argv = [sys.argv[0], "none"]
sys.path.insert(0, "/home/claude/mp/mp22")
import mp22_tg8 as m
d = json.load(open("/home/claude/mp/mp22/mp22_numbers.json"))
tv = 0.4; mod = m.Model(tv); K = d["P7"]["K with W = 65/64"]
_, _, (X1, X2, X3) = m.periods(K, mod.eps)
c = float((X1 + X2) / 2); Pi = lambda X: K + X - mod.g * X * (X - 0.75) ** 2; sn = math.sqrt(1 - c)
z = mod.f * math.sqrt(c) + 0j; p = -math.sqrt(-Pi(c)) / sn + 1j * math.sqrt(Pi(1.0)) / sn
X0 = m.chart(tv, *[np.array(v) for v in mod.labels_of(z, p, "B")])
emu = lambda l, Y: float(mod.zp(*m.labels(tv, Y), "B")[0].imag); emu.direction = 1
so = solve_ivp(m.rhs_line(tv, "B"), (0, 70 * 2 * math.pi), X0, method="DOP853", rtol=1e-12, atol=1e-13, events=emu)
te = so.t_events[0]; L = d["P7"]["on the torus W = 65/64"]["label length"]
out = {"y at the start": emu(0, X0), "first three crossings of y = 0, label length": [float(v) for v in te[:3]], "crossings with label length above 1e-6, up to the closure": int(np.sum((te > 1e-6) & (te < L + 1e-6))),
       "crossings counted by the first run": int(np.sum(te < L + 1e-6)), "label length at the closure": L, "the crossing nearest the closure, minus the closure": float(te[np.argmin(abs(te - L))] - L)}
print(json.dumps(out, indent=1))
d["P7"]["after the miss: the count of crossings"] = out
json.dump(d, open("/home/claude/mp/mp22/mp22_numbers.json", "w"), indent=1)
