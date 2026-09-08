#!/usr/bin/env python3
"""Render the frozen T2 memory-approximation diagnostics from saved results."""
from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, FixedLocator
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.directory
    result = json.loads((root / "approximation_results.json").read_text())
    selection = json.loads((root / "selection_frozen.json").read_text())
    pointwise = np.load(root / "approximation_pointwise.npz")
    models = result["models"]
    families = ["svd", "poles", "snapshots"]
    labels = ["Coupling SVD", "Selected poles", "Response snapshots"]
    colors = ["#126D85", "#C26731", "#8060A1"]
    assert all(not models[f"{family}_r{rank}"]["heldout_budget_passes"]["0.1"]
               for family in families for rank in range(5))
    assert all(selection["training_choices"]["0.01"][family]["rank"] == 5
               for family in families)

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.titlesize": 13, "axes.titleweight": "bold",
        "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#7F8B93", "text.color": "#172C3A",
        "axes.labelcolor": "#172C3A", "xtick.color": "#415765",
        "ytick.color": "#415765", "svg.fonttype": "none",
    })
    fig, axs = plt.subplots(2, 2, figsize=(15.5, 10.8))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=.075, right=.96, bottom=.14, top=.825,
                        wspace=.28, hspace=.53)
    fig.text(.075, .955, r"$T_2$ memory approximation: the tested reductions miss the budget",
             fontsize=20, weight="bold", va="top")
    fig.text(.075, .906,
             "None of the three families at r = 0–4 passes even the 10% joint gate on the held-out grids.",
             fontsize=12.5)
    fig.text(.075, .875,
             "The frozen 1% selections retain all five memory coordinates (r = 5): a full control, with no compression.",
             fontsize=11.5, color="#4B6473")

    # Panel A: fixed rank-four comparisons. Percentages are relative Frobenius errors.
    ax = axs[0, 0]
    ax.set_title("A  Four retained coordinates still miss the gates", loc="left", pad=14)
    x = np.arange(3)
    response = np.array([models[f"{f}_r4"]["spectral_maximum_relative_error"] for f in families])*100
    temporal = np.array([models[f"{f}_r4"]["time_maximum_relative_error"] for f in families])*100
    a = ax.bar(x-.17, response, width=.30, color="#126D85", label="Response error")
    b = ax.bar(x+.17, temporal, width=.30, color="#91BBC4", label="Evolution error")
    for bars in [a, b]:
        for bar in bars:
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()*1.12,
                    f"{bar.get_height():.2f}%", ha="center", va="bottom", fontsize=9)
    for budget, style in [(1,"-"),(5,":"),(10,"--")]:
        ax.axhline(budget, color="#BB4654", linestyle=style, linewidth=1.2, zorder=0)
        ax.text(2.54, budget, f"{budget}%", va="center", fontsize=9, color="#A33342")
    ax.set_yscale("log")
    ax.set_ylim(.65, 700)
    ax.set_xlim(-.55, 2.84)
    ax.set_xticks(x, ["Coupling\nSVD", "Selected\npoles", "Response\nsnapshots"])
    ax.set_ylabel("Maximum relative error (%)")
    ax.yaxis.set_major_locator(FixedLocator([1, 10, 100]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.legend(loc="upper left", frameon=False, fontsize=9)

    # Panel B: running maximum makes the registered time-prefix comparison explicit.
    ax = axs[0, 1]
    ax.set_title("B  A short interval can hide the later error", loc="left", pad=14)
    times = pointwise["times"]
    err = np.maximum.accumulate(pointwise["time_relative_svd_r4"])*100
    keep = times > 0
    ax.plot(times[keep], err[keep], color=colors[0], linewidth=2.4,
            label="SVD, r = 4: maximum so far")
    ax.axhline(1, color="#BB4654", linewidth=1.2)
    ax.axhline(10, color="#BB4654", linewidth=1.2, linestyle="--")
    for t, key, offset in [(.25,"0.25",(18,16)), (1.,"1",(15,-30)),
                            (2*np.pi,str(2*np.pi),(-85,10))]:
        val = models["svd_r4"]["time_prefixes"][key]["maximum_relative_amplitude_error"]*100
        ax.scatter([t], [val], color=colors[0], zorder=4, s=32)
        ax.annotate(f"{val:.4g}%", (t,val), xytext=offset,
                    textcoords="offset points", fontsize=10,
                    arrowprops={"arrowstyle":"-", "color":"#80909A"})
    ax.set_yscale("log")
    ax.set_ylim(.00005, 95)
    ax.set_xlim(0, 6.7)
    ax.set_xticks([0,1,np.pi,2*np.pi], ["0","1",r"$\pi$",r"$2\pi$"])
    ax.set_xlabel("Dimensionless evolution parameter t")
    ax.set_ylabel("Maximum relative amplitude error (%)")
    ax.legend(loc="lower right", fontsize=9, frameon=False)
    ax.grid(axis="y", alpha=.16)

    # Panel C: the weak discarded direct channel has much larger indirect mixing.
    ax = axs[1, 0]
    ax.set_title("C  The discarded direction still mixes indirectly", loc="left", pad=14)
    svd = models["svd_r4"]
    norms = [svd["discarded_direct_coupling_norm"], svd["retained_discarded_D_mixing_norm"]]
    ax.barh([1,0], norms, height=.36, color=["#91BBC4", "#126D85"])
    for yy, val in zip([1,0], norms):
        ax.text(val+.025, yy, f"{val:.5f}", va="center", fontsize=11, weight="bold")
    ax.set_yticks([])
    ax.text(.01, 1.28, "Direct coupling to active sector", fontsize=10)
    ax.text(.01, .28, "Mixing with retained memory coordinates", fontsize=10)
    ax.set_xlim(0,1.17)
    ax.set_ylim(-.6,1.8)
    ax.set_xlabel("Operator norm (dimensionless); SVD, r = 4")
    ax.text(.01, 1.6, f"Indirect mixing is about {norms[1]/norms[0]:.1f}× the discarded direct coupling.",
            fontsize=10, color="#4B6473", va="center")
    ax.grid(axis="x", alpha=.16)
    ax.set_axisbelow(True)

    # Panel D: response error is still above 10% at the broadest tested smoothing.
    ax = axs[1, 1]
    ax.set_title("D  Broader spectral smoothing helps, but is insufficient", loc="left", pad=14)
    etas = np.array([.05,.2,.7])
    for family, label, color in zip(families,labels,colors):
        vals = [models[f"{family}_r4"]["spectral_error_by_eta"][str(e)]*100 for e in etas]
        ax.plot(etas, vals, "o-", color=color, label=label, linewidth=2, markersize=5)
    ax.axhline(10, color="#BB4654", linestyle="--", linewidth=1.2)
    ax.text(.053, 10.7, "10% gate", color="#A33342", fontsize=9)
    svd_broad = models["svd_r4"]["spectral_error_by_eta"]["0.7"]*100
    ax.annotate(f"SVD: {svd_broad:.2f}%", (.7,svd_broad), xytext=(.22,23),
                textcoords="data", fontsize=10,
                arrowprops={"arrowstyle":"-", "color":"#80909A"})
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(7, 650)
    ax.set_xlim(.042,.88)
    ax.set_xticks(etas,["0.05","0.2","0.7"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.yaxis.set_major_locator(FixedLocator([10,100]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel(r"Imaginary part $\eta$ of $z=x+i\eta$ (dimensionless)")
    ax.set_ylabel("Maximum response error (%), r = 4")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.grid(axis="y", alpha=.16)

    full_max = max(max(models[f"{f}_r5"]["spectral_maximum_relative_error"],
                       models[f"{f}_r5"]["time_maximum_relative_error"]) for f in families)
    fig.text(.075,.068,
             f"Finite-grid diagnostics: {result['heldout_count']:,} held-out spectral points and "
             f"{result['time_count']} evolution points. All r = 5 controls have relative errors < 10⁻¹³.",
             fontsize=10.5, color="#415765")
    fig.text(.075,.040,
             "Bases were frozen using 183 training spectral points. Errors are relative Frobenius norms; "
             "this establishes no physical clock or global approximation bound.",
             fontsize=10, color="#415765")
    assert full_max < 1e-13
    fig.savefig(root / "T2_Memory_Approximation.png", dpi=160, facecolor="white")
    fig.savefig(root / "T2_Memory_Approximation.svg", facecolor="white")
    plt.close(fig)
    print("Rendered T2_Memory_Approximation.png and .svg from frozen saved results.")


if __name__ == "__main__":
    main()
