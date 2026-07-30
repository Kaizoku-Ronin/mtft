# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
import math, numpy as np, importlib.util, os
spec = importlib.util.spec_from_file_location("c", os.path.join(os.path.dirname(os.path.abspath(__file__)), "ep_census.py"))
c = importlib.util.module_from_spec(spec); spec.loader.exec_module(c)
out = open("census2.out", "w", buffering=1)

out.write("== F1: nb-stability WHERE IT MATTERS (r <= 4.8) ==\n")
for r in (3.0, 4.0, 4.8, 9.0):
    row = []
    for nb in (30, 40, 50):
        g, B = c.internal(nb=nb)
        row.append(c.count(r, g, B, npts=700).real)
    out.write(f"  r={r:4.1f}: nb30={row[0]:9.4f}  nb40={row[1]:9.4f}  "
              f"nb50={row[2]:9.4f}   drift={max(row)-min(row):8.4f}\n")

out.write("\n== F2: the staircase — off-integer counts LOCALIZE EPs ==\n")
g, B = c.internal(nb=40)
prev = None
for r in np.arange(1.0, 4.90, 0.10):
    v = c.count(r, g, B, npts=700).real
    frac = abs(v - round(v))
    flag = "  <-- EP near contour" if frac > 0.01 else ""
    if prev is not None and round(v) != round(prev):
        flag += f"  [JUMP +{round(v)-round(prev)}: EP pair modulus in ({r-0.1:.2f},{r:.2f})]"
    out.write(f"  r={r:4.2f}: {v:9.4f}{flag}\n")
    prev = v
out.write("done\n"); out.close()
