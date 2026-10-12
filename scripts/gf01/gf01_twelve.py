"""Twelve and phi: the polygon units, the two shortest loops of the modular surface, and which loop is shortest at each torus level."""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Arc
from matplotlib.lines import Line2D

BLUE, ORANGE, AQUA, INK, INK2, GRID, SURF, PALE = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#c9c8c3", "#fcfcfb", "#dddcd7"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.edgecolor": GRID, "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})
PHI = (1 + 5 ** 0.5) / 2

fig = plt.figure(figsize=(14.6, 6.4))
axA = fig.add_axes([0.01, 0.12, 0.33, 0.72])
axB = fig.add_axes([0.36, 0.10, 0.33, 0.76])
axC = fig.add_axes([0.735, 0.22, 0.25, 0.60])

# ---------------------------------------------------------------- A: the polygon units
axA.set_aspect("equal"); axA.axis("off")
polys = [(5, 2, "φ = (1+√5)/2 = 1.618", "Q(√5)", ORANGE, -2.6), (8, 3, "1+√2 = 2.414", "Q(√2)", INK2, 0.0), (12, 5, "2+√3 = 3.732", "Q(√3)", BLUE, 2.6)]
for n, k, lab, field, col, cx in polys:
    ang = np.pi / 2 + 2 * np.pi * np.arange(n) / n
    P = np.c_[cx + np.cos(ang), np.sin(ang)]
    axA.add_patch(Polygon(P, closed=True, fill=False, ec=GRID, lw=1.2))
    axA.plot([P[0, 0], P[1, 0]], [P[0, 1], P[1, 1]], color=INK, lw=2.6, solid_capstyle="round")          # the side
    axA.plot([P[0, 0], P[k, 0]], [P[0, 1], P[k, 1]], color=col, lw=2.6, solid_capstyle="round")          # the chord
    axA.text(cx, -1.22, f"{n}-gon: chord ÷ side\n{lab}\nthe unit of {field}", ha="center", va="top", fontsize=8.8, color=INK)
axA.set_xlim(-3.9, 3.9); axA.set_ylim(-2.35, 1.3)
fig.text(0.012, 0.95, "A. Each polygon has its own unit; the 12-gon's is 2+√3, not φ", fontsize=11.5, va="top")
fig.text(0.012, 0.115, "sin 15° = (√6 − √2)/4 = 0.259; the golden sine is sin 18° = (φ − 1)/2 = 0.309. The dodecagon's lengths live in Q(√2, √3),\n"
         "whose quadratic subfields are Q(√2), Q(√3) and Q(√6): φ is not among them.", fontsize=8.6, color=INK2, va="top")

# ---------------------------------------------------------------- B: the two shortest loops of the modular surface
axB.set_aspect("equal"); axB.axis("off")
th = np.linspace(np.pi / 3, 2 * np.pi / 3, 200)
axB.fill_between(np.cos(th), np.sin(th), 2.05, color="#efeeea", lw=0, zorder=0)
axB.plot(np.cos(th), np.sin(th), color=ORANGE, lw=1.4, zorder=2)                              # |tau| = 1, theta = pi class
for x in (-0.5, 0.5):
    axB.plot([x, x], [3 ** 0.5 / 2, 2.05], color=ORANGE, lw=1.4, zorder=2)                       # Re tau = +-1/2, theta = pi class
axB.plot([0, 0], [1, 2.05], color=BLUE, lw=1.4, ls=(0, (4, 3)), zorder=2)                        # Re tau = 0, theta = 0 class

def arc(cx, r, x0, x1, **kw):
    xs = np.linspace(x0, x1, 300); ys = np.sqrt(np.maximum(r * r - (xs - cx) ** 2, 0)); axB.plot(xs, ys, **kw)

R5 = 5 ** 0.5 / 2; R12 = 3 ** 0.5
arc(0.5, R5, 0.0, 0.5, color=ORANGE, lw=3.2, zorder=4); arc(-0.5, R5, -0.5, 0.0, color=ORANGE, lw=3.2, zorder=4)            # golden loop
arc(1.0, R12, -0.5, 0.5, color=BLUE, lw=2.6, zorder=3); arc(0.0, R12, -0.5, 0.5, color=BLUE, lw=2.6, zorder=3)              # twelve loop
arc(-1.0, R12, -0.5, 0.5, color=BLUE, lw=2.6, zorder=3)
pts = [((0, 1), "i\n(the half-turn)", (0.07, 0.86), BLUE), ((0.5, R5), "(1+i√5)/2\ndisc −20", (0.56, 1.03), ORANGE),
       ((-0.5, 3 ** 0.5 / 2), "ρ (the 1/3-turn)", (-0.93, 0.72), BLUE), ((0, 2 ** 0.5), "i√2  disc −8", (0.09, 1.33), BLUE),
       ((0, 3 ** 0.5), "i√3  disc −12", (0.09, 1.80), BLUE), ((0.5, 11 ** 0.5 / 2), "(1+i√11)/2\ndisc −11", (0.56, 1.60), BLUE)]
for (x, y), lab, (tx, ty), col in pts:
    axB.scatter([x], [y], s=34, color=INK, zorder=6, edgecolor=SURF, lw=1.0)
    axB.text(tx, ty, lab, fontsize=8.2, color=INK, va="center", zorder=7)
axB.text(0.0, 2.2, "the modular tile; θ = π mirrors in orange, the θ = 0 mirror dashed blue", ha="center", fontsize=8.6, color=INK2)
axB.set_xlim(-1.25, 1.25); axB.set_ylim(0.55, 2.3)
axB.legend(handles=[Line2D([], [], color=ORANGE, lw=3.2, label="the golden loop: trace 3, disc 5, length 4 log φ = 1.925"),
                    Line2D([], [], color=BLUE, lw=2.6, label="the twelve loop: trace 4, disc 12, length 2 log(2+√3) = 2.634")],
           loc="lower center", bbox_to_anchor=(0.5, -0.12), frameon=False, fontsize=8.4, handlelength=2.0)
fig.text(0.36, 0.95, "B. The two shortest closed loops of the modular surface X(1)", fontsize=11.5, va="top")

# ---------------------------------------------------------------- C: systole length at the thirteen torus levels
levels = [57, 58, 65, 77, 82, 91, 118, 123, 141, 142, 143, 145, 155]
def systole_trace(N):
    for t in range(3, 60):
        if any((a * (t - a) - 1) % N == 0 for a in range(N)):
            return t
tr = [systole_trace(N) for N in levels]; L = [2 * math.acosh(t / 2) for t in tr]
cols = [ORANGE if t == 3 else BLUE if t == 4 else PALE for t in tr]
y = np.arange(len(levels))[::-1]
axC.barh(y, L, color=cols, height=0.68, edgecolor=SURF, lw=2, zorder=2)
for yi, N, t, l in zip(y, levels, tr, L):
    axC.text(l + 0.06, yi, f"trace {t}, disc {t*t-4}", va="center", fontsize=8.2, color=INK2)
axC.set_yticks(y); axC.set_yticklabels([str(N) for N in levels], fontsize=9)
for lab in axC.get_yticklabels():
    if lab.get_text() == "143":
        lab.set_fontweight("bold")
axC.set_xlim(0, 6.9); axC.set_xlabel("length of the shortest closed loop on X0(N)", fontsize=9)
for sp_ in ("right", "top"):
    axC.spines[sp_].set_visible(False)
axC.tick_params(axis="y", length=0)
fig.text(0.735, 0.95, "C. Which loop is shortest at each torus level", fontsize=11.5, va="top")
fig.text(0.735, 0.095, "The golden loop closes at N only if 5 is a square mod every\nprime of N; (5/13) = −1, so at 143 the shortest loop is the\n"
         "twelve loop. Property tried after looking (X17, not evidence).", fontsize=8.6, color=INK2, va="top")
fig.savefig("docs/SM/figs/gf01_twelve.png", dpi=190)
print("ok", list(zip(levels, tr)))
