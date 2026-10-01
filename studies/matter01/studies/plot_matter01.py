"""Plot the archived Matter-01 reference and conditional-reduction results.

Dependencies: numpy and matplotlib. Paths resolve relative to this study bundle.
Cross-section curves evaluate the independently checked analytic Born oracle;
markers and cohomology bounds are read from the archived JSON calculations.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, ScalarFormatter


ROOT = Path(__file__).resolve().parents[1]
# Load the independent, math-only oracle directly. Importing the entire MTFT
# package would bring in unrelated research dependencies for a simple figure.
REPO = ROOT / "repo_candidate" if (ROOT / "repo_candidate").exists() else ROOT.parents[1]
_oracle_path = REPO / "src" / "mtft" / "matter" / "oracle.py"
_spec = importlib.util.spec_from_file_location("matter01_figure_oracle", _oracle_path)
assert _spec is not None and _spec.loader is not None
oracle = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(oracle)


def read_result(name: str) -> dict:
    return json.loads((ROOT / "results" / name).read_text())


def style_axis(ax, title: str):
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=14, color="#183047")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#b4c2ce")
    ax.tick_params(colors="#41576b", labelsize=10)
    ax.grid(True, axis="y", color="#e1e7ed", linewidth=.7)
    ax.set_axisbelow(True)


def plot(output: Path) -> Path:
    archive = read_result("qed_oracle_checks.json")
    vacuum = read_result("vacuum_audit.json")
    reduction = read_result("u1_reduction_checks.json")
    inputs = archive["inputs"]
    mass = inputs["electron_mass_gev"]
    alpha = 1 / inputs["alpha_inverse"]
    nb = inputs["gev2_to_nb"]
    assert math.isclose(mass, oracle.MASS_GEV, rel_tol=0, abs_tol=1e-20)
    assert math.isclose(alpha, oracle.ALPHA, rel_tol=0, abs_tol=1e-18)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.labelcolor": "#41576b", "figure.facecolor": "#f7f9fb", "axes.facecolor": "#ffffff", "savefig.facecolor": "#f7f9fb"})
    fig, axs = plt.subplots(2, 2, figsize=(13.8, 10.7))
    fig.subplots_adjust(left=.085, right=.96, top=.85, bottom=.18, hspace=.64, wspace=.26)
    fig.suptitle("Matter-01 · Start with U(1)", x=.085, y=.975, ha="left", fontsize=25, fontweight="bold", color="#183047")
    fig.text(.085, .933, "Reference QED, a conditional gauge–geometry bridge, and a spectrum audit", fontsize=14, color="#536a7d")
    fig.text(.085, .899, f"INPUTS  electron mass = {mass*1000:.11f} MeV   ·   α⁻¹ = {inputs['alpha_inverse']:.9f}   ·   unpolarized tree level", fontsize=10.5, color="#87631d")

    ax = axs[0, 0]
    style_axis(ax, "A  Photon angles: finite-mass event cross section")
    cosine = np.linspace(-1, 1, 1201)
    colors = ["#237fb3", "#bd5f36", "#167b63"]
    for beta, color in zip((.2, .7, .95), colors):
        y = np.array([oracle.dsigma_dcosine(beta, float(c), mass, alpha) * nb / 1e6 for c in cosine])
        ax.plot(cosine, y, color=color, lw=2.2, label=rf"$\beta={beta:.2f}$")
    ax.set_yscale("log")
    ax.set_xlim(-1, 1)
    ax.set_xticks([-1, -.5, 0, .5, 1])
    ax.set_xlabel(r"Photon direction $\cos\theta$ in the CM frame")
    ax.set_ylabel(r"$d\sigma_{\rm event}/d\cos\theta$  [mb]")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=1, loc="upper center", ncol=3, fontsize=10)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.text(0, -.25, "The ½ identical-photon factor is included over the full sphere.\nForward/backward peaks remain finite when the electron mass is kept.", transform=ax.transAxes, fontsize=9.5, color="#536a7d", va="top")

    ax = axs[0, 1]
    style_axis(ax, "B  Total event cross section above threshold")
    energies = np.geomspace(2 * mass / math.sqrt(1-.1**2), .02, 700)
    betas = np.sqrt(1-(2*mass/energies)**2)
    totals = np.array([oracle.sigma_total(float(b), mass, alpha) * nb / 1e6 for b in betas])
    ax.loglog(energies*1000, totals, lw=2.2, color="#237fb3", label="Massive Born formula")
    samples = archive["total_grid"]
    ax.scatter([s["sqrt_s_gev"]*1000 for s in samples], [s["event_sigma_nb"]/1e6 for s in samples], marker="o", s=36, edgecolors="white", linewidths=.8, color="#bd5f36", zorder=4, label="Archived independent checks")
    threshold = 2*mass*1000
    ax.axvline(threshold, color="#9aaabb", linestyle="--", lw=1)
    ax.set_xlim(threshold*.98, 21)
    ax.set_xticks([1.1, 2, 5, 10, 20])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.set_xlabel(r"Total CM energy $\sqrt{s}$  [MeV]")
    ax.set_ylabel(r"Event cross section $\sigma$  [mb]")
    ax.legend(frameon=False, fontsize=9.5, loc="upper right")
    ax.text(0, -.25, f"Pair threshold: 2mₑ = {threshold:.6f} MeV.  1 mb = 10⁶ nb.\nMeasured inputs; curves are reference QED, not MTFT predictions.", transform=ax.transAxes, fontsize=9.5, color="#536a7d", va="top")

    ax = axs[1, 0]
    style_axis(ax, "C  Conditional C3X gauge–breathing reduction")
    x = np.linspace(-.6, .6, 401)
    Z = np.exp(x)
    g = np.exp(-x/2)
    ax.plot(x, Z, lw=2.2, color="#237fb3", label=r"$Z(\chi)=e^{\chi/M_P}$")
    ax.plot(x, g, lw=2.2, color="#167b63", label=r"$g_{\rm eff}/g_0=e^{-\chi/(2M_P)}$")
    ax.axhline(1, color="#9aaabb", linestyle="--", lw=1)
    ax.axvline(0, color="#dde5ec", lw=1)
    witness = reduction["C3X_breathing"]
    exponent = witness["canonical_exponent"]
    assert exponent == 1.0 and witness["internal_dimension"] == 2
    ax.scatter([.1, .1], [witness["Z"], witness["g_eff_over_g0"]], c=["#237fb3", "#167b63"], s=38, edgecolors="white", linewidths=.8, zorder=4)
    ax.set_xlim(-.6, .6)
    ax.set_ylim(.43, 1.98)
    ax.set_xlabel(r"Uniform breathing field $\chi/M_P$")
    ax.set_ylabel("Dimensionless coefficient / coupling ratio")
    ax.legend(frameon=False, loc="upper left", fontsize=10)
    ax.text(0, -.26, "d = 2 internal surface; uniform classical ansatz; constant χ background.\nThis plots a conditional formula, not coupled gauge–matter evolution.", transform=ax.transAxes, fontsize=9.5, color="#536a7d", va="top")

    ax = axs[1, 1]
    style_axis(ax, "D  Net indices can hide large paired cohomology")
    selected = {r["orientation"]: r for r in vacuum["witness"]["RR_cohomology_lower_bounds"]}
    control = {r["orientation"]: r for r in vacuum["collinear_comparison_control"]["RR_cohomology_lower_bounds"]}
    rows = [selected["c/L"], control["c/L"], selected["L/B"], control["L/B"]]
    lower1 = [int(r["h1_L_lower_bound"]) for r in rows]
    lower2 = [int(r["h1_dual_lower_bound"]) for r in rows]
    pos = np.arange(4)
    width = .31
    bars1 = ax.bar(pos-width/2, lower1, width=width, color="#237fb3", label=r"$h^1(\mathcal{L})$ lower bound")
    bars2 = ax.bar(pos+width/2, lower2, width=width, color="#bd5f36", label=r"$h^1(\mathcal{L}^{\vee})$ lower bound")
    for bars in (bars1, bars2):
        for bar in bars:
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()*1.11, str(int(bar.get_height())), ha="center", va="bottom", fontsize=9, color="#183047")
    ax.set_yscale("log")
    ax.set_ylim(1, 4500)
    ax.set_xticks(pos, ["Selected 31\nc/L", "Control 36\nc/L", "Selected 31\nL/B", "Control 36\nL/B"])
    ax.set_ylabel("Riemann–Roch cohomology lower bound")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.legend(frameon=False, fontsize=9.2, loc="upper left")
    ax.text(0, -.26, "31 / 36 are index-based extras scores, not particle counts.\nBars are lower bounds, not full dimensions or a physical massless spectrum.", transform=ax.transAxes, fontsize=9.5, color="#536a7d", va="top")

    fig.text(.085, .042, "Sources: archived Matter-01 oracle, U(1) reduction, and selected-vacuum audits.\nReference QED does not derive photon masslessness or electron masses in the MTFT vacuum.", fontsize=9, color="#536a7d", linespacing=1.6)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "viz" / "Matter01_U1_Overview.png")
    args = parser.parse_args()
    print(plot(args.output))
