#!/usr/bin/env python3
"""Render the W143 reduction diagnostics from the recorded JSON only."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parent


def main():
    data = json.loads((ROOT / "effective_operator_results.json").read_text())
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10.5,
        "axes.titlesize": 12, "axes.labelsize": 10.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#9ba5af", "axes.labelcolor": "#283644",
        "text.color": "#1e2d3a", "xtick.color": "#485866", "ytick.color": "#485866",
        "svg.fonttype": "none", "savefig.facecolor": "white",
    })
    teal, orange, purple = "#007f86", "#d76a20", "#6b55a1"
    fig, axes = plt.subplots(1, 3, figsize=(16.6, 6.5))
    fig.subplots_adjust(left=.055, right=.976, top=.735, bottom=.258, wspace=.32)
    fig.text(.055, .935, r"$W_{143}$: retaining feedback from the complementary sector", fontsize=22, weight="bold")
    fig.text(.055, .885, "13 complex dimensions  •  two coupled 2 × 2 blocks  •  nine uncoupled −1 directions", fontsize=12, color="#526472")
    fig.text(.055, .837, "Numerical geometry and exact involution identities; dimensionless matrix experiments", fontsize=11, color="#526472")

    rows = sorted((r for r in data["regular_response"] if r["z_imag"] == .05), key=lambda r: r["z_real"])
    x = np.array([r["z_real"] for r in rows])
    full = np.array([r["modes"][0]["spectral_response_full"] for r in rows])
    eff = np.array([-r["modes"][0]["effective_imag"] / np.pi for r in rows])
    naive = np.array([r["modes"][0]["spectral_response_naive"] for r in rows])
    ax = axes[0]
    ax.set_title("A  Restoring the projected response", loc="left", pad=27, weight="bold")
    ax.text(0, 1.035, r"Mode 1; $z=x+0.05i$", transform=ax.transAxes, color="#526472", fontsize=10)
    ax.semilogy(x, full, color="#233744", lw=2.4, label="Full operator")
    ax.semilogy(x, eff, color=teal, lw=0, marker="o", markersize=4.0,
                markerfacecolor="white", markeredgewidth=1.2, markevery=12,
                label="With feedback (overlaps)")
    ax.semilogy(x, naive, color=orange, ls="--", lw=1.9, label="Active block alone")
    for pole in [-1, 1]:
        ax.axvline(pole, lw=.9, color="#97a3ac", ls=":", zorder=0)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(.0015, 10)
    ax.set_xlabel("Spectral parameter x")
    ax.set_ylabel(r"Projected response $-\mathrm{Im}\,G/\pi$")
    ax.set_xticks([-1.5, -1, 0, 1, 1.5])
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="white", framealpha=.93, fontsize=9)
    ax.text(.97, .06, "Dotted lines: full poles ±1", transform=ax.transAxes, ha="right", color="#526472", fontsize=9)

    ax = axes[1]
    ax.set_title("B  Transfer to the complement", loc="left", pad=27, weight="bold")
    ax.text(0, 1.035, r"Chosen control $U(t)=\exp(-itW_{143})$", transform=ax.transAxes, color="#526472", fontsize=10)
    ctl = data["dimensionless_unitary_control"]
    t = np.array(ctl["grid_t"])
    for j, color in enumerate([teal, purple]):
        f = np.array(ctl["mode_complement_norm_fractions"][j])
        peak = data["modes"][j]["maximum_toy_complement_norm_fraction"]
        ax.plot(t, f, color=color, lw=2.2, label=f"Mode {j+1}: peak {100*peak:.4f}%")
        ax.axhline(peak, color=color, lw=.85, alpha=.5, ls=":")
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(0, .15)
    ax.set_xlabel("Chosen dimensionless parameter t")
    ax.set_ylabel("Squared norm fraction in complement")
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi], ["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.set_yticks([0, .04, .08, .12])
    ax.legend(loc="upper center", fontsize=9, frameon=False)

    ax = axes[2]
    ax.set_title("C  Stability near compression poles", loc="left", pad=27, weight="bold")
    ax.text(0, 1.035, r"$z=d_j+i\varepsilon$; maximum over the two poles", transform=ax.transAxes, color="#526472", fontsize=10)
    stresses = data["pole_stress"]
    eps = sorted({r["z_imag"] for r in stresses}, reverse=True)
    raw = [max(r["effective_relative_error"] for r in stresses if r["z_imag"] == e) for e in eps]
    continued = [max(r["involution_relative_error"] for r in stresses if r["z_imag"] == e) for e in eps]
    ax.loglog(eps, raw, color=orange, lw=2.1, marker="o", markersize=4, label="Raw Schur evaluation")
    ax.loglog(eps, continued, color=teal, lw=2.1, marker="s", markersize=3.6, label="Continued involution formula")
    ax.set_xlim(1e-2, 1e-12)
    ax.set_ylim(1e-15, 1e-3)
    ax.set_xticks([1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12])
    ax.set_xlabel(r"Imaginary offset $\varepsilon$")
    ax.set_ylabel("Relative error against full projected inverse")
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    ax.annotate(r"$6.59\times10^{-5}$", xy=(eps[-1], raw[-1]), xytext=(-10, 13),
                textcoords="offset points", ha="right", fontsize=10, color=orange)
    ax.text(.97, .215, r"Continued: $\leq 6.55\times10^{-14}$", transform=ax.transAxes, ha="right", fontsize=10, color=teal)
    for ax in axes:
        ax.grid(axis="y", which="major", color="#e5e9ec", lw=.7, zorder=0)
        ax.tick_params(axis="both", labelsize=9)

    s = data["summary"]
    fig.text(.055, .156,
             f"Regular-grid feedback error ≤ {s['max_regular_effective_relative_error']:.2e} across {s['regular_point_count']} points.  "
             f"Schur denominator condition number reaches {s['max_stress_Schur_condition']:.2e} in the pole stress test.",
             fontsize=11, color="#263e4e")
    fig.text(.055, .102,
             "These are matrix diagnostics: ±1 are operator eigenvalues, not identified physical energies; t has no identified physical clock.",
             fontsize=10.5, color="#526472")
    fig.text(.055, .059,
             "MTFT v0.26.2  •  X₀(143)  •  Source: effective_operator_results.json  •  No continuum or particle interpretation inferred",
             fontsize=9.5, color="#667581")
    for suffix in ("png", "svg"):
        fig.savefig(ROOT / f"W143_Effective_Operator.{suffix}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
