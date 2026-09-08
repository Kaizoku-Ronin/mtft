#!/usr/bin/env python3
"""Render the frozen T2 experiment from its recorded JSON, without rerunning it."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parent


def main():
    data = json.loads((ROOT / "t2_effective_operator_results.json").read_text())
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.titlesize": 13, "axes.labelsize": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#9ba5af", "axes.labelcolor": "#283644",
        "text.color": "#1e2d3a", "xtick.color": "#485866", "ytick.color": "#485866",
        "svg.fonttype": "none", "savefig.facecolor": "white",
    })
    teal, orange, purple, ink = "#007f86", "#d76a20", "#6b55a1", "#233744"
    fig, axs = plt.subplots(2, 2, figsize=(15.6, 11.6))
    fig.subplots_adjust(left=.079, right=.965, top=.79, bottom=.155, wspace=.23, hspace=.59)
    fig.text(.079, .948, r"$T_2$: coupled dynamics beyond independent pairs", fontsize=24, weight="bold")
    fig.text(.079, .899, "8 active + 5 complementary complex dimensions  •  all five complementary directions couple", fontsize=12.5, color="#526472")
    fig.text(.079, .867, "The singular-vector pairs do not form invariant blocks. Three input directions are initially uncoupled, but develop delayed transfer.", fontsize=11.5, color="#526472")

    ax = axs[0, 0]
    rows = sorted((r for r in data["regular_response"] if r["z_imag"] == .05), key=lambda r: r["z_real"])
    x = np.array([r["z_real"] for r in rows])
    full = np.array([r["average_spectral_response_full"] for r in rows])
    naive = np.array([r["average_spectral_response_naive"] for r in rows])
    ax.set_title("A  Removing feedback changes the response", loc="left", pad=29, weight="bold")
    ax.text(0, 1.035, r"Average over the active eight; $z=x+0.05i$", transform=ax.transAxes, color="#526472", fontsize=10.5)
    ax.plot(x, full, color=teal, lw=2.1, label="Full projected response")
    ax.plot(x, naive, color=orange, ls="--", lw=1.8, label="Active block alone")
    ax.set_xlim(-3, 3)
    ax.set_ylim(0, max(full.max(), naive.max()) * 1.28)
    ax.set_xlabel("Dimensionless spectral parameter x")
    ax.set_ylabel(r"$-\mathrm{Im}\,\mathrm{Tr}(G)/(8\pi)$")
    ax.legend(loc="upper left", frameon=False, fontsize=10)

    ax = axs[0, 1]
    rows = data["dimensionless_dynamics"]
    t = np.array([r["t"] for r in rows])
    avg = np.array([r["average_transfer"] for r in rows])
    dark = np.array([r["dark_average_transfer"] for r in rows])
    ax.set_title("B  Initially uncoupled inputs also transfer", loc="left", pad=29, weight="bold")
    ax.text(0, 1.035, r"Chosen evolution $U(t)=\exp(-itT_2)$; averages over orthonormal inputs", transform=ax.transAxes, color="#526472", fontsize=10)
    ax.plot(t, avg, color=teal, lw=2.2, label="All 8 active inputs")
    ax.plot(t, dark, color=purple, lw=2.2, label="3 initially uncoupled inputs")
    for key, color, label, shift in [
        ("sampled_max_average_transfer", teal, "37.15%", (12, 9)),
        ("sampled_max_dark_average_transfer", purple, "25.80%", (-14, -23)),
    ]:
        item = data["summary"][key]
        ax.scatter([item["t"]], [item["value"]], color=color, s=24, zorder=5)
        ax.annotate(label, xy=(item["t"], item["value"]), xytext=shift,
                    textcoords="offset points", color=color, fontsize=10,
                    ha="left" if shift[0] > 0 else "right")
    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(0, .495)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi], ["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
    ax.set_yticks([0, .1, .2, .3, .4])
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.set_xlabel("Chosen dimensionless evolution parameter t")
    ax.set_ylabel("Average complementary squared norm fraction")
    ax.legend(loc="upper right", frameon=False, fontsize=9.5)
    ax.text(.985, .035, "Labels mark sampled maxima", transform=ax.transAxes, ha="right", fontsize=9, color="#526472")

    ax = axs[1, 0]
    rows = data["short_time_controls"]
    st = np.array([r["t"] for r in rows])
    c2 = data["structure"]["ordinary_short_time_coefficient"]
    c4 = data["structure"]["dark_short_time_coefficient"]
    ordinary_ratio = np.array([r["ordinary_average_over_t_squared"] for r in rows]) / c2
    dark_ratio = np.array([r["dark_average_over_t_fourth"] for r in rows]) / c4
    ax.set_title("C  Delayed transfer has a different leading order", loc="left", pad=29, weight="bold")
    ax.text(0, 1.035, "Analytic expansions; exponents are derived, not fitted", transform=ax.transAxes, color="#526472", fontsize=10.5)
    ax.axhline(1, color="#8c9aa4", ls=":", lw=1.1)
    ax.semilogx(st, ordinary_ratio, color=teal, lw=2.1, marker="o", ms=5,
                label=r"Active average / $(c_2t^2)$")
    ax.semilogx(st, dark_ratio, color=purple, lw=2.1, marker="s", ms=4.5,
                label=r"Initially uncoupled average / $(c_4t^4)$")
    ax.set_xlim(.105, .0029)
    ax.set_ylim(.9875, 1.0012)
    ax.set_xlabel("Short-time parameter t (decreasing →)")
    ax.set_ylabel("Measured / predicted leading term")
    ax.set_xticks([.1, .025, .00625], ["0.1", "0.025", "0.00625"])
    ax.minorticks_off()
    ax.set_yticks([.988, .992, .996, 1])
    ax.legend(loc="lower right", frameon=False, fontsize=9.5)
    ax.text(.49, .60, f"c₂ = {c2:.9f}\nc₄ = {c4:.9f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=10.5, color=ink,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .9, "pad": 2.5})

    ax = axs[1, 1]
    rows = data["pole_stress"]
    eps = sorted({r["z_imag"] for r in rows}, reverse=True)
    raw = [max(r["effective_relative_error"] for r in rows if r["z_imag"] == e) for e in eps]
    stable = [max(r["spectral_relative_error"] for r in rows if r["z_imag"] == e) for e in eps]
    ax.set_title("D  Near a compression pole, use the full spectrum", loc="left", pad=29, weight="bold")
    ax.text(0, 1.035, r"$z=d_j+i\varepsilon$; maximum error over all five complementary poles", transform=ax.transAxes, color="#526472", fontsize=10)
    ax.loglog(eps, raw, color=orange, lw=2.1, marker="o", ms=4, label="Raw Schur evaluation")
    ax.loglog(eps, stable, color=teal, lw=2.1, marker="s", ms=3.8, label="Full spectral reconstruction")
    ax.set_xlim(1e-2, 1e-12)
    ax.set_ylim(1e-15, 3e-3)
    ax.set_xticks([1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12])
    ax.set_yticks([1e-15, 1e-12, 1e-9, 1e-6, 1e-3])
    ax.set_xlabel(r"Imaginary offset $\varepsilon$")
    ax.set_ylabel("Relative error against full projected inverse")
    ax.legend(loc="upper left", frameon=False, fontsize=9.5)
    ax.annotate(r"$1.64\times10^{-4}$", xy=(eps[-1], raw[-1]), xytext=(-10, 12),
                textcoords="offset points", ha="right", fontsize=10.5, color=orange)
    ax.text(.97, .17, r"Spectral: $\leq 1.49\times10^{-14}$", transform=ax.transAxes,
            ha="right", fontsize=10.5, color=teal)

    for ax in axs.flat:
        ax.grid(axis="y", which="major", color="#e5e9ec", lw=.75, zorder=0)
        ax.tick_params(axis="both", labelsize=10)

    summary = data["summary"]
    fig.text(.079, .092,
             f"Maximum feedback error ≈ {summary['max_regular_effective_relative_error']:.2e} on the {summary['regular_count']:,}-point regular grid.",
             fontsize=11.5, color=ink)
    fig.text(.079, .058, "Numerical matrix diagnostics: T₂ and t have no identified physical energy or clock. All reported transfer maxima are sampled values.",
             fontsize=10.5, color="#526472")
    fig.text(.079, .029, "MTFT v0.26.2  •  X₀(143)  •  Source: t2_effective_operator_results.json  •  Frozen protocol; no fitted scale or physical interpretation",
             fontsize=9.5, color="#667581")
    for extension in ("png", "svg"):
        fig.savefig(ROOT / f"T2_Effective_Operator.{extension}", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
