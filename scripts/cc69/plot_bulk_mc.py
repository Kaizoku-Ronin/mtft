"""Hysteresis scans of the repaired SU(2) action -> figure + summary.  usage: python plot_bulk_mc.py out17 fig.png"""
import glob
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D, OUT = sys.argv[1], sys.argv[2]
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.titlesize": 10, "axes.titleweight": "bold",
                     "axes.titlecolor": INK, "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6})
runs = sorted((json.load(open(f)) for f in glob.glob(os.path.join(D, "mc_k*.json"))), key=lambda r: r["kappa"])
fig, axs = plt.subplots(1, 2, figsize=(11.4, 4.4), layout="constrained")
ax = axs[0]
summary = {}
for i, r in enumerate(runs):
    up = sorted(r["hot_up"], key=lambda q: q["beta"]); dn = sorted(r["cold_down"], key=lambda q: q["beta"])
    bu = np.array([q["beta"] for q in up]); pu = np.array([q["P"] for q in up]); pd = np.array([q["P"] for q in dn])
    gap = pd - pu; j = int(np.argmax(np.abs(gap)))
    dP = np.diff(pu); jmax = int(np.argmax(dP))
    summary[r["kappa"]] = {"beta_A": r["beta_A_equivalent"], "max_hysteresis": float(gap[j]), "at_beta": float(bu[j]),
                           "largest_step_up": float(dP[jmax]), "step_at": float((bu[jmax] + bu[jmax + 1]) / 2),
                           "P_beta0_hot": float(pu[0]), "P_beta0_cold": float(pd[0])}
    ax.plot(bu, pu, "-", color=C[i], lw=1.5, label=f"κ = {r['kappa']:g}  (β_A ≈ {r['beta_A_equivalent']:.2f})")
    ax.plot(bu, pd, "--", color=C[i], lw=1.2)
    ax.fill_between(bu, pu, pd, color=C[i], alpha=0.18, lw=0)
ax.set(xlabel="β  (fundamental coupling)", ylabel="plaquette ⟨½ Tr U_p⟩", xlim=(0, 3), ylim=(-0.02, 0.9))
ax.text(0.1, 0.84, "solid: β up from a hot start;  dashed: β down from a cold start;  band: hysteresis", color=INK2, fontsize=8)
ax.legend(loc="lower right", fontsize=7.5, frameon=False)
ax.set_title("A  SU(2), 6⁴ lattice: the tower opens a first-order bulk transition")
ax = axs[1]
ks = np.array(sorted(summary)); hy = np.array([summary[k]["max_hysteresis"] for k in ks]); st = np.array([summary[k]["largest_step_up"] for k in ks])
ax.plot(ks, hy, "o-", color=C[0], lw=1.5, ms=6, mec=SURF, mew=1.4, label="largest hysteresis gap in ⟨P⟩")
ax.plot(ks, st, "s-", color=C[1], lw=1.5, ms=6, mec=SURF, mew=1.4, label="largest single-step jump of ⟨P⟩ (Δβ = 0.125)")
ax.axvline(1, color=AXIS, lw=0.9); ax.text(1.4, 0.3, "package default κ = 1", color=INK2, fontsize=8, rotation=90, va="center")
ax.set(xlabel="κ  (tower strength)", ylabel="size of the discontinuity", xlim=(-2, max(ks) + 5), ylim=(0, max(0.35, hy.max() * 1.15)))
ax.legend(loc="upper left", fontsize=7.5, frameon=False)
ax.set_title("B  Where the transition appears")
fig.savefig(OUT, dpi=170)
print(json.dumps(summary, indent=1))
