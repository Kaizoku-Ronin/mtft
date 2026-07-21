#!/usr/bin/env python3
"""
Generate viz/hero_stiffness.png and viz/stiffness_navigator.html
=================================================================
The Yang-Mills stiffness landscape mu_N(y) over gauge groups SU(N),
rendered from the arithmetic weights w_n = sum_{d|n} (log d)/d.

Every point on this surface is computed from integers — zero free
parameters. The ridge structure at small y is confinement; the
uniform positivity IS the mass gap statement.

Reproducible: python viz/make_hero.py   (regenerates both assets)
Requires: numpy, matplotlib (PNG); the HTML is dependency-free plotly-CDN.
"""
import math
import numpy as np

# ── arithmetic weights (the whole input) ──────────────────────
NMAX = 4000
w = np.zeros(NMAX + 1)
for d in range(2, NMAX + 1):
    w[d::d] += math.log(d) / d

def mu(N, y):
    n = np.arange(2, NMAX + 1)
    return np.sum(n * n * w[2:] * np.exp(-2 * np.pi * y * n)
                  * (1 - np.cos(2 * np.pi * n / N)))

Y_C = 0.18174  # canonical confinement depth (CriticalDepths.y_conf)

ys = np.linspace(0.02, 0.50, 160)
Ns = np.arange(2, 17)
Z = np.array([[math.log10(max(mu(N, y), 1e-300)) for y in ys] for N in Ns])

# ── PNG (GitHub-dark Mathematical Brutalism) ──────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

BG = "#0d1117"
fig = plt.figure(figsize=(16, 9), facecolor=BG)
ax = fig.add_subplot(111, projection="3d", facecolor=BG)
Yg, Ng = np.meshgrid(ys, Ns)
surf = ax.plot_surface(Yg, Ng, Z, cmap=cm.magma, rstride=1, cstride=2,
                       linewidth=0, antialiased=True, alpha=0.97)
# confinement-depth trace
zc = [math.log10(mu(N, Y_C)) for N in Ns]
ax.plot([Y_C] * len(Ns), Ns, zc, color="#58e6d9", lw=2.5, zorder=10)
ax.text(Y_C, Ns[-1] + 0.4, zc[-1], "y_c = 0.18174", color="#58e6d9",
        fontsize=11, family="monospace")

ax.set_xlabel("modular depth  y", color="#e6edf3", labelpad=10)
ax.set_ylabel("gauge group  SU(N)", color="#e6edf3", labelpad=10)
ax.set_zlabel("log10  mu_N(y)", color="#e6edf3", labelpad=8)
for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
    pane.set_facecolor(BG); pane.set_edgecolor("#21262d")
ax.tick_params(colors="#8b949e")
ax.grid(False)
ax.view_init(elev=27, azim=-58)
fig.text(0.5, 0.945, "THE STIFFNESS LANDSCAPE", ha="center",
         color="#e6edf3", fontsize=22, family="monospace", weight="bold")
fig.text(0.5, 0.905,
         "mu_N(y) = sum n^2 w_n e^(-2 pi y n) (1 - cos 2 pi n/N)   —   "
         "zero free parameters, computed from the integers",
         ha="center", color="#8b949e", fontsize=11, family="monospace")
fig.text(0.5, 0.045, "MTFT  ·  X0(143)  ·  pip install mtft",
         ha="center", color="#484f58", fontsize=10, family="monospace")
plt.savefig("hero_stiffness.png", dpi=140, facecolor=BG,
            bbox_inches="tight", pad_inches=0.35)
print("hero_stiffness.png written")

# ── interactive HTML (plotly CDN, self-contained data) ────────
ys_h = np.linspace(0.02, 0.50, 90)
Zh = [[round(math.log10(max(mu(N, y), 1e-300)), 4) for y in ys_h] for N in Ns]
html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>MTFT — Stiffness Navigator</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>body{{margin:0;background:#0d1117;font-family:ui-monospace,monospace}}
#hdr{{color:#e6edf3;text-align:center;padding:14px 0 2px;font-size:20px;font-weight:bold}}
#sub{{color:#8b949e;text-align:center;font-size:12px;padding-bottom:6px}}
#plot{{width:100vw;height:86vh}}</style></head><body>
<div id="hdr">THE STIFFNESS LANDSCAPE &mu;<sub>N</sub>(y)</div>
<div id="sub">drag to rotate &middot; scroll to zoom &middot; zero free parameters &middot; MTFT / X&#8320;(143)</div>
<div id="plot"></div><script>
const ys={list(map(float, np.round(ys_h,4)))},Ns={list(map(int,Ns))},Z={Zh};
Plotly.newPlot('plot',[{{type:'surface',x:ys,y:Ns,z:Z,colorscale:'Magma',
 contours:{{z:{{show:true,usecolormap:true,project:{{z:true}}}}}},
 colorbar:{{title:'log10 mu',tickfont:{{color:'#8b949e'}},titlefont:{{color:'#8b949e'}}}}}},
 {{type:'scatter3d',mode:'lines',x:Array(Ns.length).fill({Y_C}),y:Ns,
   z:Ns.map((n,i)=>Z[i][{int(np.argmin(np.abs(ys_h-Y_C)))}]),
   line:{{color:'#58e6d9',width:6}},name:'y_c = 0.18174'}}],
 {{paper_bgcolor:'#0d1117',scene:{{xaxis:{{title:'modular depth y',color:'#8b949e',gridcolor:'#21262d'}},
   yaxis:{{title:'SU(N)',color:'#8b949e',gridcolor:'#21262d'}},
   zaxis:{{title:'log10 mu_N',color:'#8b949e',gridcolor:'#21262d'}},
   camera:{{eye:{{x:-1.6,y:-1.5,z:0.9}}}}}},showlegend:false,
   margin:{{l:0,r:0,t:0,b:0}}}});</script></body></html>"""
with open("stiffness_navigator.html", "w") as f:
    f.write(html)
print("stiffness_navigator.html written (%d KB)" % (len(html)//1024))
