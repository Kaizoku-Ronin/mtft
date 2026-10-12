"""PG-01 figure: the two energies of the primon gas, the shape well, the volume runaway, the Wilson-line torus, the toy run's phase portrait."""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpmath import mp, zeta, log as mlog, diff as mdiff

from mtft.matter import vacuum_atlas as VA
from mtft.matter.path_a import _lattice

BLUE, ORANGE, AQUA, INK, INK2, GRID, SURF, PALE = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#c9c8c3", "#fcfcfb", "#dddcd7"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.5, "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.edgecolor": GRID, "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})
mp.dps = 20

fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.6))
fig.subplots_adjust(left=0.05, right=0.985, top=0.9, bottom=0.08, wspace=0.3, hspace=0.42)
(axA, axB, axC), (axD, axE, axF) = axes
for ax in axes.flat:
    for sp_ in ("right", "top"):
        ax.spines[sp_].set_visible(False)
    ax.grid(color=GRID, lw=0.5, alpha=0.6)

# ---------------------------------------------------------------- A: the two energies
betas = np.linspace(1.03, 6.0, 160)
UD = [float(-mdiff(lambda x: mlog(-zeta(x) * zeta(x + 1, derivative=1)), b)) for b in betas]
UU = [float(-zeta(b, derivative=1) / zeta(b)) for b in betas]
axA.plot(betas, UD, color=BLUE, lw=2.0, label="marked ensemble ⟨log n⟩_D, Z_D = −ζ(β)ζ′(β+1)")
axA.plot(betas, UU, color=ORANGE, lw=2.0, label="unmarked gas −ζ′/ζ = Σ Λ(n) n^{−β}")
axA.axvline(1.0, color=INK2, lw=1.0, ls=(0, (4, 3)))
axA.set_yscale("log"); axA.set_xlim(0.9, 6.0); axA.set_xlabel("inverse temperature β"); axA.set_ylabel("mean energy")
axA.text(1.12, 0.03, "Hagedorn wall β = 1: the weights\nare wₙ = wₙ(β = 1)", fontsize=8.5, color=INK2, va="bottom")
axA.legend(loc="upper right", frameon=False, fontsize=8)
axA.set_title("A. Two energies of the primon gas", loc="left", fontsize=10.5)

# ---------------------------------------------------------------- B: the shape well
sw = VA.shape_well(); phis = float(np.log(float(sw["r_star"])))
ph = np.linspace(phis - 2.2, phis + 2.2, 300)
axB.plot(ph, np.cosh(ph - phis), color=BLUE, lw=2.2)
axB.scatter([phis], [1.0], s=40, color=INK, zorder=5)
axB.annotate("floor: the self-dual point r* = √(19/11)\ncurvature exactly 1 in units of E*", xy=(phis, 1.0), xytext=(phis + 0.25, 3.0), fontsize=8.5,
             arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
axB.set_xlabel("φ = log r,  r = A_X / A_E"); axB.set_ylabel("E(r) / E*"); axB.set_ylim(0.6, 5.2)
axB.set_title("B. The shape modulus: E/E* = cosh(φ − φ*)", loc="left", fontsize=10.5)

# ---------------------------------------------------------------- C: the volume runaway
L = np.linspace(-0.5, 3.0, 300)
V = np.exp(-2 * L) + np.exp(-1.5 * L)
axC.plot(L, V, color=ORANGE, lw=2.2)
axC.set_xlabel("log Vol"); axC.set_ylabel("V(Vol) ∝ Vol⁻² + Vol⁻³ᐟ²  (shape fixed)")
axC.text(0.9, 1.3, "every term falls with the volume:\nno finite critical point — the flow\nruns to decompactification (PV-01)", fontsize=8.5, color=INK2)
axC.set_title("C. The volume has no critical point", loc="left", fontsize=10.5)

# ---------------------------------------------------------------- D: the Wilson-line torus
lat = _lattice((0, 1), 40)
t1 = np.linspace(-0.5, 0.5, 121); t2 = np.linspace(-0.5, 0.5, 121); T1, T2 = np.meshgrid(t1, t2)
Vw = np.zeros_like(T1)
for m, n, r2, odd in lat:
    if odd:
        Vw += -2 * r2 ** -4 * np.cos(2 * np.pi * (T1 * m + T2 * n))
Vw = Vw / abs(Vw).max()
im = axD.pcolormesh(T1, T2, Vw, cmap="Blues_r", shading="auto")
axD.scatter([0], [0], s=50, color=INK, zorder=5); axD.scatter([0], [0.5], s=50, color=ORANGE, zorder=5); axD.scatter([0], [-0.5], s=50, color=ORANGE, zorder=5)
axD.text(0.04, 0.02, "jet vacuum θ = 0: sink", fontsize=8.5, color=INK)
axD.text(0.04, 0.44, "θ_η: source", fontsize=8.5, color=ORANGE)
axD.set_xlabel("θ₁ (periods)"); axD.set_ylabel("θ₂ (periods)"); axD.set_aspect("equal"); axD.grid(False)
axD.set_title("D. The Wilson-line potential on its torus (PA-01)", loc="left", fontsize=10.5)

# ---------------------------------------------------------------- E: the toy run's phase portrait
for H, col, lab in ((0.2, BLUE, "H = 0.2 m: spiral attract (Type 5)"), (1.0, ORANGE, "H = 1.0 m: node (overdamped)")):
    run = VA.toy_run(H, t_end=50.0)
    axE.plot(run["phi"] - run["phi_star"], run["phidot"], color=col, lw=1.4, label=lab)
axE.scatter([0], [0], s=40, color=INK, zorder=5)
axE.axhline(0, color=GRID, lw=0.8); axE.axvline(0, color=GRID, lw=0.8)
axE.set_xlabel("φ − φ*"); axE.set_ylabel("dφ/dt"); axE.legend(loc="upper right", frameon=False, fontsize=8)
axE.set_title("E. Toy run with Hubble friction (DIAGNOSTIC)", loc="left", fontsize=10.5)

# ---------------------------------------------------------------- F: time series, shape captured and volume running
run = VA.toy_run(0.4, mode="radiation", t_end=80.0); vc = VA.volume_companion(t_end=80.0)
axF.plot(run["t"] - run["t"][0], run["phi"] - run["phi_star"], color=BLUE, lw=1.6, label="shape φ − φ* (radiation era H = 1/2t): captured")
axF.plot(vc["t"], vc["log_vol"], color=ORANGE, lw=1.6, label="log Vol (same friction): runs away, ~ (2/3) log t")
axF.axhline(0, color=GRID, lw=0.8)
axF.set_xlabel("time (units of 1/m)"); axF.legend(loc="center right", bbox_to_anchor=(1.0, 0.56), frameon=False, fontsize=8)
axF.set_title("F. Shape captured, volume not (DIAGNOSTIC)", loc="left", fontsize=10.5)

fig.text(0.05, 0.955, "PG-01 — energy and the two clocks: the weights are the divisor-gas energy at the Hagedorn wall; the vacuum's shape has a sink, its volume has none",
         fontsize=11.5, va="center")
fig.savefig("docs/SM/figs/pg01_energy_clocks.png", dpi=170)
print("ok")
