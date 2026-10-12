"""Forces first: the dots-and-lines picture of the parent U(8), the strong force's lattice Z[omega], and the ladder in string language."""
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch

BLUE, ORANGE, AQUA, INK, INK2, GRID, SURF = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#c9c8c3", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": INK, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "savefig.facecolor": SURF})

fig = plt.figure(figsize=(14.6, 6.5))
axA = fig.add_axes([0.005, 0.10, 0.40, 0.78])
axB = fig.add_axes([0.415, 0.16, 0.25, 0.68])
axC = fig.add_axes([0.68, 0.03, 0.315, 0.85])
for ax in (axA, axB, axC):
    ax.set_aspect("equal"); ax.axis("off")

# ---------------------------------------------------------------- A: eight dots, grouped 4 + 2 + 2
groups = {"SU(4)": [(-2.0, 0.9), (-1.0, 1.25), (-1.0, -0.15), (-2.0, -0.5)],
          "SU(2)_L": [(1.55, 1.75), (2.45, 1.25)],
          "SU(2)_R": [(1.55, -1.05), (2.45, -0.55)]}
labels_dot = {"SU(4)": ["c", "c", "c", "A1"], "SU(2)_L": ["L", "L"], "SU(2)_R": ["B", "A2"]}
pts = [(g, p) for g, ps in groups.items() for p in ps]
for (g1, p1), (g2, p2) in itertools.combinations(pts, 2):
    if g1 == g2:
        axA.plot([p1[0], p2[0]], [p1[1], p2[1]], color=BLUE, lw=2.6, zorder=3, solid_capstyle="round")
    else:
        axA.plot([p1[0], p2[0]], [p1[1], p2[1]], color=ORANGE, lw=1.0, alpha=0.55, zorder=1)
for g, ps in groups.items():
    for p, lab in zip(ps, labels_dot[g]):
        axA.add_patch(Circle(p, 0.2, fill=False, ec=BLUE, lw=1.4, zorder=4))
        axA.add_patch(Circle(p, 0.085, color=INK, zorder=5))
        axA.text(p[0] + 0.24, p[1] + 0.2, lab, fontsize=8.5, color=INK2, zorder=6)
axA.text(-1.5, 2.05, "SU(4)\ncolour (c × 3) + lepton (A1)", ha="center", fontsize=9.5)
axA.text(2.0, 2.35, "SU(2)_L  (weak)", ha="center", fontsize=9.5)
axA.text(2.0, -1.65, "SU(2)_R", ha="center", fontsize=9.5)
axA.text(0.25, 0.95, "(4, 2, 1)\nleft-handed\nquarks + leptons", fontsize=8.5, color=INK, ha="center",
         bbox=dict(boxstyle="round,pad=0.2", fc=SURF, ec="none", alpha=0.9), zorder=7)
axA.text(0.25, -0.75, "(4̄, 1, 2)\nright-handed\nquarks + leptons", fontsize=8.5, color=INK, ha="center",
         bbox=dict(boxstyle="round,pad=0.2", fc=SURF, ec="none", alpha=0.9), zorder=7)
axA.text(2.75, 0.35, "(1, 2, 2)\nHiggs", fontsize=8.5, color=INK, ha="left",
         bbox=dict(boxstyle="round,pad=0.2", fc=SURF, ec="none", alpha=0.9), zorder=7)
axA.set_xlim(-2.8, 3.9); axA.set_ylim(-2.75, 2.75)
axA.plot([], [], color=BLUE, lw=2.6, label="line inside a group: a force carrier (8 lines × 2 orientations = 16, plus 8 circles = 24)")
axA.plot([], [], color=ORANGE, lw=1.4, alpha=0.7, label="line between groups: a particle (20 lines × 2 = 40 states: one Pati–Salam family + the Higgs)")
axA.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.06), frameon=False, fontsize=8.6, handlelength=2.2)
fig.text(0.012, 0.018, "The programme's working reading: every field is a U(8)-adjoint block. PV-01's non-supersymmetric parent puts the\n"
         "particles on the companion lines eᵢ + eⱼ of Λ²(8) instead — which reading holds is open question I.1.", fontsize=8, color=INK2)
fig.text(0.012, 0.95, "A. Eight dots (the parent's U(8)): choosing the groups chooses the forces,\n    and the lines between groups are the particles",
         fontsize=11, va="top")

# ---------------------------------------------------------------- B: the strong force on three dots: Z[omega]
w = np.exp(2j * np.pi / 3)
lat = [a + b * w for a in range(-4, 5) for b in range(-4, 5)]
lat = [z for z in lat if abs(z) < 2.25]
axB.scatter([z.real for z in lat], [z.imag for z in lat], s=7, color=GRID, zorder=1)
roots = [(1 - w) * u for u in (1, -1, w, -w, w * w, -w * w)]
hexa = sorted(roots, key=lambda z: np.angle(z))
hx = [z.real for z in hexa] + [hexa[0].real]; hy = [z.imag for z in hexa] + [hexa[0].imag]
axB.plot(hx, hy, color=BLUE, lw=0.8, alpha=0.5, zorder=2)
axB.scatter([z.real for z in roots], [z.imag for z in roots], s=70, color=BLUE, zorder=4)
axB.scatter([0], [0], s=70, color=BLUE, zorder=4); axB.scatter([0.13], [0], s=70, color=BLUE, zorder=4, ec=SURF)
quarks = [1, w, w * w]
axB.plot([z.real for z in quarks] + [1], [z.imag for z in quarks] + [0], color=ORANGE, lw=1.2, zorder=3)
axB.scatter([z.real for z in quarks], [z.imag for z in quarks], s=80, color=ORANGE, zorder=5, marker="^")
axB.scatter([-z.real for z in quarks], [-z.imag for z in quarks], s=80, facecolor=SURF, edgecolor=ORANGE, lw=1.6, zorder=5, marker="v")
axB.text(1.08, 0.05, "1", fontsize=9); axB.text(w.real - 0.32, w.imag + 0.05, "ω", fontsize=9); axB.text((w * w).real - 0.32, (w * w).imag - 0.2, "ω²", fontsize=9)
axB.set_xlim(-2.3, 2.3); axB.set_ylim(-2.3, 2.3)
axB.text(0, -2.55, "quarks ▲ at 1, ω, ω² (triality 1); antiquarks ▽ (triality 2);\n8 gluons ● on the hexagon and at the centre (triality 0).\n"
         "Only triality 0 escapes confinement: Z[ω]/(1 − ω) = Z/3.\nCharge ≡ ⅔ × triality (mod 1), so every free particle\ncarries a whole charge.", ha="center", va="top", fontsize=8.6, color=INK2)
fig.text(0.418, 0.95, "B. Three dots: the strong force's charges\n    live on the Eisenstein lattice Z[ω]", fontsize=11, va="top")

# ---------------------------------------------------------------- C: the ladder in string language
ys = [3.35, 2.0, 0.55, -1.05]
axC.add_patch(Circle((0, ys[0]), 0.09, color=INK)); axC.add_patch(Circle((0, ys[0]), 0.22, fill=False, ec=BLUE, lw=1.4))
axC.add_patch(Circle((-0.45, ys[1]), 0.09, color=INK)); axC.add_patch(Circle((0.45, ys[1]), 0.09, color=INK))
axC.plot([-0.45, 0.45], [ys[1], ys[1]], color=BLUE, lw=2.6)
axC.add_patch(Circle((0, ys[2]), 0.32, fill=False, ec=AQUA, lw=2.4))
axC.add_patch(Ellipse((0, ys[3]), 1.25, 0.62, fill=False, ec=INK, lw=1.8))
axC.add_patch(Ellipse((0, ys[3] + 0.04), 0.5, 0.17, fill=False, ec=INK, lw=1.4))
import textwrap
txt = ["a dot: a brane, carrying one circle, U(1)",
       "two dots and a line: a string between them, a force carrier. Drawn on top of each other the line is massless and U(1)×U(1) grows to U(2)",
       "a loop with no ends: a closed string, the graviton. Gravity is not a line between dots",
       "a loop swept round a loop: a torus. Summing over its shapes τ (the modular tile) works only for even self-dual charge lattices: in the heterotic string the torus picks the forces: E8×E8 or Spin(32)/Z2"]
txt = ["\n".join(textwrap.wrap(t, 45)) for t in txt]
for y, t in zip(ys, txt):
    axC.text(0.75, y, t, va="center", fontsize=8.6, linespacing=1.25)
for y0, y1 in zip(ys[:-1], ys[1:]):
    axC.add_patch(FancyArrowPatch((0, y0 - 0.38), (0, y1 + 0.42), arrowstyle="-|>", mutation_scale=10, color=INK2, lw=1.0))
axC.set_xlim(-0.8, 5.4); axC.set_ylim(-2.05, 3.85)
fig.text(0.685, 0.95, "C. Your ladder, in string language\n    (known physics, not an MTFT result)", fontsize=11, va="top")
fig.savefig("docs/SM/figs/gf01_dots.png", dpi=190)
print("ok")
