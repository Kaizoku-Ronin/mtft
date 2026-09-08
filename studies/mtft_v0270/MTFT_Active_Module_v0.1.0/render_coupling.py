"""Render registered coupling data without fitting or changing the data.

Usage: python render_coupling.py
Inputs are resolved beside this script. All plotted fractions use the first
registered run (period-input dps=50, frame seed=143). Matrix arithmetic remains
float64; the period-input precision is not an operator error certificate.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np


ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT / "coupling_results.json").read_text())
kernel = json.loads((ROOT / "kernel_check_results.json").read_text())
run = data["runs"][0]
names = data["primary"] + data["replication"]
keys = ["active_to_active", "active_to_fixed", "fixed_to_active", "fixed_to_fixed"]
labels = ["Active → active", "Active → fixed", "Fixed → active", "Fixed → fixed"]
colors = ["#197C88", "#E9A04A", "#BE5B3A", "#536781"]
values = np.array([[run["operators"][n]["block_squared_fractions"][k]
                    for k in keys] for n in names])
if not np.allclose(values.sum(axis=1), 1, atol=1e-12, rtol=0):
    raise ValueError("Squared block fractions do not sum to one.")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelcolor": "#344455", "text.color": "#20303F",
    "xtick.color": "#526170", "ytick.color": "#344455",
    "svg.fonttype": "none", "savefig.facecolor": "white",
})
fig = plt.figure(figsize=(15.8, 10.1), facecolor="white")
fig.text(.055, .95, "Arithmetic coupling across the 8 + 5 split",
         fontsize=24, weight="bold", va="top")
fig.text(.055, .905, "MTFT v0.26.2 · X₀(143) · registered operator sample",
         fontsize=12.5, color="#526170")
fig.text(.055, .872,
         "Active: 8 complex = 16 real dimensions     |     Fixed: 5 complex = 10 real dimensions",
         fontsize=12.5)

ax = fig.add_axes([.085, .20, .60, .57])
y = np.arange(len(names), dtype=float)
y[8:] += .55
left = np.zeros(len(names))
for j, color in enumerate(colors):
    widths = 100 * values[:, j]
    ax.barh(y, widths, left=left, color=color, edgecolor="white", linewidth=.9,
            height=.65, zorder=3)
    for yi, l, w in zip(y, left, widths):
        if w >= 5:
            ax.text(l + w/2, yi, f"{w:.1f}", ha="center", va="center",
                    fontsize=10, color="#20303F" if j == 1 else "white")
    left += widths
display = {n: ("STAR" if n == "STAR" else rf"${n[0]}_{{{n[1:]}}}$") for n in names}
ax.set_yticks(y, [display[n] for n in names], fontsize=13)
ax.set_ylim(y[-1] + .6, -.8)
ax.set_xlim(0, 100)
ax.set_xticks(np.arange(0, 101, 20))
ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
ax.set_xlabel("Share of the operator’s squared Frobenius norm", labelpad=12)
ax.grid(axis="x", color="#E5E9ED", linewidth=.8, zorder=0)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_color("#CBD3DA")
ax.tick_params(axis="y", length=0, pad=10)
ax.axhline(7.77, color="#CCD4DB", linewidth=.8)
ax.text(0, 8.00, "REPLICATION OPERATORS", va="center", fontsize=8.5,
        color="#697988")
for yi, row in zip(y, values):
    mixed = 100 * (row[1] + row[2])
    label = "< 10⁻²⁰ %" if mixed < 1e-20 else f"{mixed:.1f}%"
    ax.text(102.3, yi, label, ha="left", va="center", fontsize=10.5,
            color="#687989" if mixed < 1e-20 else "#8F432B", clip_on=False)
ax.text(102.3, -.77, "Cross-sector\nshare", fontsize=10, weight="bold",
        ha="left", va="bottom", clip_on=False)

fig.legend([Patch(facecolor=c) for c in colors], labels, loc="upper left",
           bbox_to_anchor=(.078, .83), ncol=4, frameon=False, fontsize=10.5,
           columnspacing=1.8, handlelength=1.6)

# The hull is computed in the 26-dimensional real Hodge frame. The second
# value is dim(span(V_active, T2 V_active)), not the rank of T2 itself.
hull = data["hulls"]["T2_only"]["growth"]
ah = fig.add_axes([.805, .49, .155, .265])
ah.bar([0, 1], hull[:2], color=[colors[0], colors[3]], width=.56)
ah.set_ylim(0, 29)
ah.set_yticks([0, 8, 16, 26])
ah.set_ylabel("Real dimension", fontsize=10)
ah.set_xticks([0, 1], ["Active", "Active +\nT₂(active)"], fontsize=10)
ah.grid(axis="y", color="#E5E9ED", linewidth=.8)
ah.set_axisbelow(True)
ah.spines["left"].set_visible(False)
ah.spines["bottom"].set_color("#CBD3DA")
ah.tick_params(axis="y", length=0)
for i, v in enumerate(hull[:2]):
    ah.text(i, v + .55, str(v), ha="center", va="bottom", fontsize=14, weight="bold")
ah.set_title("T₂ alone spans\nthe full space", loc="left", fontsize=13,
             weight="bold", pad=17)

fig.text(.805, .397, "Reading the map", fontsize=12, weight="bold")
fig.text(.805, .375,
         "STAR is consistent with\nsector preservation.\n\n"
         "W₁₄₃ mixes less than W₁₁\nand W₁₃, but its mixing\nis numerically resolved.",
         fontsize=10.5, color="#526170", va="top", linespacing=1.45)

gap = kernel["Cij_svd"]["gap"]
fig.text(.055, .113,
         "Squared Frobenius norm shares in a Hodge orthonormal frame; these are not transition probabilities.",
         fontsize=11, weight="bold")
fig.text(.055, .080,
         f"Direct seed kernel: five complex null directions; retained/discarded singular-value gap ≈ {gap:.1e}.\n"
         "Numerical diagnostic (float64). This figure does not establish exact invariance, a spacetime signature, or particle dynamics.",
         fontsize=10, color="#526170", linespacing=1.5)

for extension in ("png", "svg"):
    path = ROOT / f"Arithmetic_Coupling_Map.{extension}"
    fig.savefig(path, dpi=180)
    print(path)
plt.close(fig)
