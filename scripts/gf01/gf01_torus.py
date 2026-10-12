"""The force torus: couplings as its shape (electroweak Pythagoras), and the climb of the 8-torus's sheets into U(8), SO(16) and E8."""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Arc, FancyArrowPatch

BLUE, ORANGE, AQUA, INK, INK2, GRID, SURF, PALE = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#c9c8c3", "#fcfcfb", "#dddcd7"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.edgecolor": GRID, "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})

# ---------------------------------------------------------------- inputs (MS-bar at m_Z) and derived lengths
ainv, s2 = 127.95, 0.2312
e = math.sqrt(4 * math.pi / ainv); g = e / math.sqrt(s2); gp = e / math.sqrt(1 - s2)
X, Y = 1 / g, 1 / gp                     # one unit of T3 is 1/g long, one unit of Y is 1/g' long
thW = math.degrees(math.atan2(X, Y))     # angle of the photon diagonal from the Y axis
Yps = X * math.sqrt(5 / 3)               # same width, Pati-Salam equal couplings: g'^2 = (3/5) g^2
thPS = math.degrees(math.atan2(X, Yps))
assert abs(math.hypot(X, Y) - 1 / e) < 1e-12 and abs(math.sin(math.radians(thPS)) ** 2 - 3 / 8) < 1e-12

fig = plt.figure(figsize=(14.6, 6.6))
axD = fig.add_axes([0.02, 0.17, 0.43, 0.70])
axE1 = fig.add_axes([0.54, 0.50, 0.44, 0.34])
axE2 = fig.add_axes([0.54, 0.13, 0.44, 0.24])

# ---------------------------------------------------------------- D: the electroweak torus
axD.set_aspect("equal"); axD.axis("off")
axD.add_patch(Rectangle((0, 0), X, Y, fill=False, ec=BLUE, lw=1.8, zorder=2))
axD.add_patch(Rectangle((0, 0), X, Yps, fill=False, ec=INK2, lw=1.2, ls=(0, (4, 3)), zorder=2))
axD.plot([0, X], [0, Yps], color=INK2, lw=1.2, ls=(0, (4, 3)), zorder=2)
axD.add_patch(FancyArrowPatch((0, 0), (X, Y), arrowstyle="-|>", mutation_scale=16, color=ORANGE, lw=2.6, zorder=4))
zx, zy = Y, -X; n = math.hypot(zx, zy); zx, zy = 1.35 * zx / n, 1.35 * zy / n
axD.add_patch(FancyArrowPatch((0, 0), (zx, zy), arrowstyle="-|>", mutation_scale=12, color=INK2, lw=1.3, ls=(0, (3, 2)), zorder=3))
axD.add_patch(Arc((0, 0), 1.25, 1.25, theta1=90 - thW, theta2=90, color=INK, lw=1.0, zorder=3))
axD.text(-0.08, 0.66, f"θ_W = {thW:.1f}°", fontsize=9.5, ha="right", va="center")
axD.text(X * 0.70, -0.07, f"T₃: 1/g = {X:.2f}", ha="center", va="top", fontsize=9, color=INK2)
axD.text(-0.1, Y * 0.62, f"Y: 1/g′ = {Y:.2f}", ha="right", va="center", fontsize=9, color=INK2)
axD.text(-1.7, 3.12, "one unit of weak isospin T₃ is 1/g long, one unit of hypercharge Y is 1/g′ long (at m_Z)", ha="left", va="bottom", fontsize=9, color=INK2)
axD.text(X + 0.14, Y - 0.08, "the photon: the diagonal Q = T₃ + Y\n"
         f"its length 1/e = {1/e:.2f} = √({X:.2f}² + {Y:.2f}²)\n"
         f"its tilt is the weak angle: sin²θ_W = {s2} at m_Z", ha="left", va="top", fontsize=9)
axD.text(X + 0.14, Yps - 0.02, "same width, Pati–Salam's equal couplings\n(g₄ = g_L = g_R at matching):\nsin²θ_W = 3/8, θ_W = 37.8°",
         ha="left", va="top", fontsize=9, color=INK2)
axD.text(zx + 0.12, zy - 0.02, "the Z: the perpendicular direction,\nthe one the Higgs breaks", ha="left", va="center", fontsize=9, color=INK2)
axD.set_xlim(-1.75, 4.9); axD.set_ylim(-1.2, 3.4)
fig.text(0.012, 0.95, "D. The force torus has a shape, and the shape is the coupling", fontsize=11.5, va="top")
import textwrap
capD = ("Lengths are measured with the couplings, so the photon's strength is Pythagoras: 1/e² = 1/g² + 1/g′². GC-01's matching rule "
        "α_Y⁻¹ = α_R⁻¹ + ⅔ α₄⁻¹ is the same theorem on the Pati–Salam torus, hypercharge being the diagonal Y = T₃R + √⅔ T₁₅. With its "
        "θ-angle each coupling is a point τ = θ/2π + 4πi/g² of the upper half-plane, the arena of TR-01's mirrors.")
fig.text(0.012, 0.13, "\n".join(textwrap.wrap(capD, 98)), fontsize=8.8, color=INK2, va="top")

# ---------------------------------------------------------------- E1: the rungs of the 8-torus
ks = np.arange(9); counts = [math.comb(8, int(k)) for k in ks]
cols = [AQUA if k % 2 == 0 else PALE for k in ks]
axE1.bar(ks, counts, width=0.72, color=cols, edgecolor=SURF, lw=2, zorder=2)
for k, c in zip(ks, counts):
    axE1.text(k, c + 2.0, str(c), ha="center", va="bottom", fontsize=9, color=INK2)
axE1.annotate("the 8 circles:\nthe dots", xy=(1, 4), xytext=(0.2, 46), fontsize=8.6, color=INK, ha="center",
              arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
axE1.annotate("28 sheets = pairs of\ncircles: the lines", xy=(1.68, 14), xytext=(1.25, 68), fontsize=8.6, color=INK, ha="center",
              arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
axE1.set_xticks(ks); axE1.set_xlim(-0.6, 8.6); axE1.set_ylim(0, 86); axE1.set_yticks([])
axE1.set_xlabel("k: sheets spanned by k of the eight circles  (the rung Λᵏ(8))", fontsize=9)
for sp in ("left", "right", "top"):
    axE1.spines[sp].set_visible(False)
from matplotlib.patches import Patch
axE1.legend(handles=[Patch(color=AQUA, label="even k: 1 + 28 + 70 + 28 + 1 = 128"), Patch(color=PALE, label="odd k: 8 + 56 + 56 + 8 = 128")],
            loc="upper right", frameon=False, fontsize=8.6)
fig.text(0.535, 0.95, "E. The climb: dots, lines and sheets of the 8-torus, and where it ends", fontsize=11.5, va="top")

# ---------------------------------------------------------------- E2: E8 as dots + lines + even sheets
segs = [(8, INK2, "8", "circles"), (56, BLUE, "56", "lines eᵢ − eⱼ"), (56, ORANGE, "56", "lines ±(eᵢ + eⱼ)"), (128, AQUA, "128", "even sheets ½(±1, …, ±1)")]
x0 = 0
for w, c, num, lab in segs:
    axE2.barh(0, w, left=x0, height=0.5, color=c, edgecolor=SURF, lw=2, zorder=2)
    axE2.text(x0 + w / 2, -0.42, f"{num}\n{lab}" if w > 10 else num, ha="center", va="top", fontsize=8.6, color=INK)
    x0 += w
for (end, y, lab) in [(64, 0.45, "U(8) = 64: our parent's forces and particles"), (120, 0.95, "SO(16) = 120: + Λ²(8) and its conjugate, PV-01's particle lines"),
                      (248, 1.45, "E8 = 248: + the 128 even sheets")]:
    axE2.plot([0, end], [y, y], color=INK2, lw=1.0)
    axE2.plot([end, end], [y - 0.08, y], color=INK2, lw=1.0)
    axE2.plot([0, 0], [y - 0.08, y], color=INK2, lw=1.0)
    axE2.text(1, y + 0.05, lab, va="bottom", fontsize=8.8, color=INK)
axE2.set_xlim(-2, 260); axE2.set_ylim(-1.25, 1.85); axE2.axis("off")
capE = ("The rungs are Λᵏ(8), the k-sheets of the 8-torus as the moves see them. E8's 248 is all of the dots, all 112 lines and the 128 even "
        "sheets, and its root lattice is the only even self-dual lattice in eight dimensions. Exact and standard; not evidence that the parent is E8.")
fig.text(0.535, 0.075, "\n".join(textwrap.wrap(capE, 92)), fontsize=8.8, color=INK2, va="top")

fig.savefig("docs/SM/figs/gf01_torus.png", dpi=190)
print("ok", round(X, 4), round(Y, 4), round(1 / e, 4), round(thW, 2), round(thPS, 2))
