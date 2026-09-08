"""Render the diagnostic geometry beside exact packaged-matrix sign data."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import numpy as np


ROOT = Path(__file__).resolve().parent


def main():
    diagnostic = json.loads((ROOT / "block_coupling_results.json").read_text())
    exact = json.loads((ROOT / "exact_blocks_results.json").read_text())
    run = diagnostic["runs"][0]
    names = diagnostic["block_order"]
    assert names == ["ell", "old", "q4", "q6"]
    operators = ["W11", "W143"]
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 12,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.labelcolor": "#27334A", "text.color": "#18273F",
        "xtick.color": "#27334A", "ytick.color": "#27334A",
        "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(18, 8.5), facecolor="#FAFCFF")
    grid = fig.add_gridspec(1, 3, left=.055, right=.925, bottom=.24,
                           top=.72, width_ratios=[1.05, 1, 1], wspace=.32)
    fig.text(.055, .94, "Arithmetic blocks and the local active sector", fontsize=25,
             fontweight="bold", va="top")
    fig.text(.055, .879,
             "The tested arithmetic operators preserve every Hecke block exactly. "
             "The local active projector P crosses those blocks.", fontsize=14)

    ax = fig.add_subplot(grid[0, 0])
    fractions = np.array([run["geometry"]["blocks"][b]["active_fraction_of_block"]
                          for b in names])
    x = np.arange(4)
    ax.bar(x, fractions * 100, width=.75, color="#138C93", label="Active overlap")
    ax.bar(x, (1-fractions) * 100, bottom=fractions * 100, width=.75,
           color="#DCE3EC", label="Fixed overlap")
    for i, f in enumerate(fractions):
        ax.text(i, f*50, f"{f*100:.2f}%", ha="center", va="center",
                color="white", fontsize=12, fontweight="bold")
    ax.set_xticks(x, [f"{b}\ndim {exact['blocks'][b]['dimension_complex']}" for b in names])
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 25, 50, 75, 100], ["0%", "25%", "50%", "75%", "100%"])
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#E4E9F0", linewidth=.8)
    ax.set_title("Overlap within each Hecke block", loc="left", fontsize=15, pad=34)
    ax.text(0, 1.035, "Trace fraction; dimensions below are complex", transform=ax.transAxes,
            fontsize=11)
    ax.legend(loc="lower left", bbox_to_anchor=(-.02, -.22), ncol=2,
              frameon=False, fontsize=11, columnspacing=1.1, handlelength=1.3)

    norm = Normalize(0, 17)
    cmap = plt.get_cmap("YlOrBr")
    image = None
    for j, op in enumerate(operators):
        ax = fig.add_subplot(grid[0, j+1])
        info = run["operators"][op]
        matrix = np.array(info["commutator_block_squared_fractions"]) * 100
        labels = []
        for b in names:
            signs = exact["operators"][op]["block_data"][b]
            sign = ("+" if signs["minus_eigenspace_dimension"] == 0 else
                    "−" if signs["plus_eigenspace_dimension"] == 0 else "±")
            labels.append(f"{b} ({sign})")
        image = ax.imshow(matrix, cmap=cmap, norm=norm)
        ax.set_xticks(x, labels, fontsize=11)
        ax.set_yticks(x, labels, fontsize=11)
        ax.tick_params(length=0)
        for row in range(4):
            for col in range(4):
                ax.text(col, row, f"{matrix[row,col]:.2f}", ha="center", va="center",
                        fontsize=12, color="white" if matrix[row,col] > 10 else "#18273F")
        for pos in [.5, 1.5, 2.5]:
            ax.axhline(pos, color="white", linewidth=2)
            ax.axvline(pos, color="white", linewidth=2)
        total = info["commutator_norm_squared_over_operator_squared"] * 100
        ax.set_title(f"{op}: commutator contributions", loc="left", fontsize=15, pad=34)
        ax.text(0, 1.045, f"Sum: {total:.2f}%  |  Coupling rank: {info['coupling_rank_real']} real",
                transform=ax.transAxes, fontsize=11)
        ax.text(.5, -.18, "Exact operator signs shown in parentheses", transform=ax.transAxes,
                ha="center", fontsize=10.5)
        if op == "W143":
            for row, col in [(2, 3), (3, 2)]:
                ax.add_patch(Rectangle((col-.48, row-.48), .96, .96,
                                       fill=False, ec="#138C93", lw=2.5))
    color_ax = fig.add_axes([.945, .30, .012, .37])
    cb = fig.colorbar(image, cax=color_ax, ticks=[0, 4, 8, 12, 16])
    cb.set_label("Common scale: % of ‖A‖²F", fontsize=11, labelpad=10)

    fig.text(.42, .82, r"Heatmap entry: $100\,\|E_i[A,P]E_j\|_F^2\,/\,\|A\|_F^2$",
             fontsize=14)
    fig.text(.055, .105,
             "Why W143 couples less: q4 and q6 both have sign −, so their commutator blocks vanish exactly (teal outlines).",
             fontsize=13, fontweight="bold")
    fig.text(.055, .063,
             "Overlap is neither an intersection dimension nor a probability. Heatmaps describe [A,P], not arithmetic transitions between Hecke blocks.",
             fontsize=11)
    fig.text(.055, .025,
             "All plotted magnitudes and coupling ranks are numerical diagnostics in a metric-orthonormal frame (seed 143). "
             "Signs and block preservation are exact on packaged matrices.\n"
             "Entries shown as 0.00 are rounded; only the stated sign cancellations are asserted exact here. MTFT v0.26.2.",
             fontsize=10, color="#55627A", linespacing=1.35)
    fig.savefig(ROOT / "Hecke_Block_Map.png", dpi=160, facecolor=fig.get_facecolor())
    fig.savefig(ROOT / "Hecke_Block_Map.svg", facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    main()
