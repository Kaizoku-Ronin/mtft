"""
MTFT Visualization Helpers
===========================

Optional plotting functions for 2D and 3D visualization.
Works with matplotlib (2D/3D), plotly (interactive 3D), or
exports raw numpy arrays for any other library.

Install extras:  pip install mtft[full]  (adds matplotlib + scipy)
For interactive:  pip install plotly

All functions return raw data (numpy arrays / dicts) by default.
Pass plot=True to render with matplotlib.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from mtft.arithmetic import weight_array, S1, stiffness_S, center_stiffness, mass_gap_stiffness
from mtft.constants import (
    PI, FEIGENBAUM_DELTA, CriticalDepths, GAUGE, T_INF, LEPTONS,
)


# ═══════════════════════════════════════════════════════════════
#  Stiffness Landscape (2D)
# ═══════════════════════════════════════════════════════════════

def stiffness_landscape(
    y_range: tuple = (0.01, 0.5),
    n_points: int = 200,
    N_values: tuple = (2, 3, 4, 5),
    plot: bool = False,
) -> dict:
    """
    Holonomy stiffness μ_N(y) across gauge groups and modular depth.

    Returns dict with 'y', and 'mu_N' for each N.
    Pass plot=True to render with matplotlib.
    """
    y_arr = np.linspace(y_range[0], y_range[1], n_points)
    data = {"y": y_arr}

    for N in N_values:
        mu = np.array([mass_gap_stiffness(y, N=N, n_max=200) for y in y_arr])
        data[f"mu_{N}"] = mu

    if plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))
        for N in N_values:
            ax.semilogy(y_arr, np.maximum(data[f"mu_{N}"], 1e-20), label=f"SU({N})")

        ax.axvline(CriticalDepths.y_conf, color="red", ls="--", alpha=0.7, label="y_conf")
        ax.axvline(CriticalDepths.y_spec, color="blue", ls="--", alpha=0.7, label="y_spec")
        ax.set_xlabel("Modular depth y")
        ax.set_ylabel("μ_N(y)  [mass gap stiffness]")
        ax.set_title("MTFT Yang-Mills Confinement Landscape")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        return {**data, "fig": fig}

    return data


# ═══════════════════════════════════════════════════════════════
#  Hosotani Potential (2D)
# ═══════════════════════════════════════════════════════════════

def hosotani_potential_plot(
    fermion_fraction: float = 0.4,
    n_points: int = 500,
    plot: bool = False,
) -> dict:
    """
    One-loop effective potential V_eff(θ) for the holonomy angle.
    """
    from mtft.hosotani import HosotaniPotential

    hp = HosotaniPotential(fermion_fraction=fermion_fraction)
    theta = np.linspace(0.01, PI - 0.01, n_points)
    V = hp(theta)
    theta0 = hp.find_vacuum()
    data = {"theta": theta, "V": V, "theta0": theta0}

    if plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(theta, V, "b-", lw=2)
        ax.axvline(theta0, color="red", ls="--", label=f"vacuum θ₀={theta0:.3f}")
        ax.set_xlabel("θ_H (holonomy angle)")
        ax.set_ylabel("V_eff(θ)")
        ax.set_title("Hosotani Effective Potential (Gauge-Higgs Unification)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        return {**data, "fig": fig}

    return data


# ═══════════════════════════════════════════════════════════════
#  Rotation Curve (2D)
# ═══════════════════════════════════════════════════════════════

def rotation_curve_plot(
    A: float = 1.0,
    r_range: tuple = (0.1, 50.0),
    n_points: int = 300,
    M_baryon: float = 0.5,
    plot: bool = False,
) -> dict:
    """
    Galaxy rotation curve from τ-vortex dark halo.

    Uses dimensionless units for demonstration.
    v²_∞ = 2πGA² gives the flat asymptote.
    """
    r = np.linspace(r_range[0], r_range[1], n_points)

    # Baryonic (Keplerian) contribution
    v_baryon = np.sqrt(M_baryon / r)

    # τ-vortex halo contribution: v² = 2πGA² (flat)
    v_tau = np.full_like(r, A)

    # Combined
    v_total = np.sqrt(v_baryon ** 2 + v_tau ** 2)

    data = {"r": r, "v_baryon": v_baryon, "v_tau": v_tau, "v_total": v_total}

    if plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(r, v_baryon, "b--", label="Baryonic (Keplerian)", alpha=0.7)
        ax.plot(r, v_tau, "r:", label="τ-vortex halo", alpha=0.7)
        ax.plot(r, v_total, "k-", lw=2, label="Total (MTFT)")
        ax.set_xlabel("Radius r")
        ax.set_ylabel("Circular velocity v(r)")
        ax.set_title("Flat Rotation Curve from τ-Vortex Dark Halo")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        return {**data, "fig": fig}

    return data


# ═══════════════════════════════════════════════════════════════
#  Koide Manifold (3D)
# ═══════════════════════════════════════════════════════════════

def koide_manifold_3d(
    n_phi: int = 200,
    n_r: int = 5,
    r_values: tuple = (0.5, 1.0, 2.0, 5.0, 10.0),
    plot: bool = False,
) -> dict:
    """
    The Koide manifold K ≅ S¹ × ℝ₊ in mass-root space (√m₁, √m₂, √m₃).

    Each circle at fixed r satisfies K = 2/3 exactly.
    The lepton point (e, μ, τ) sits on one specific circle.
    """
    from mtft.koide import koide_manifold_point

    phi_arr = np.linspace(0, 2 * PI, n_phi)
    curves = {}

    for r in r_values:
        pts = np.array([koide_manifold_point(phi, r) for phi in phi_arr])
        curves[f"r={r}"] = pts

    # Lepton point
    lepton_pt = np.array([
        math.sqrt(LEPTONS.e), math.sqrt(LEPTONS.mu), math.sqrt(LEPTONS.tau)
    ])

    data = {"phi": phi_arr, "curves": curves, "lepton_point": lepton_pt}

    if plot:
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection="3d")

        for label, pts in curves.items():
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], alpha=0.6, label=label)

        ax.scatter(*lepton_pt, color="red", s=100, zorder=5, label="(e, μ, τ)")

        ax.set_xlabel("√m₁")
        ax.set_ylabel("√m₂")
        ax.set_zlabel("√m₃")
        ax.set_title("Koide Manifold K ≅ S¹ × ℝ₊")
        ax.legend(fontsize=8)
        plt.tight_layout()
        return {**data, "fig": fig}

    return data


# ═══════════════════════════════════════════════════════════════
#  Burning Ship Fractal (2D image)
# ═══════════════════════════════════════════════════════════════

def burning_ship_plot(
    re_range: tuple = (-2.0, 1.5),
    im_range: tuple = (-2.0, 0.5),
    resolution: int = 600,
    max_iter: int = 200,
    plot: bool = False,
) -> dict:
    """
    Burning Ship fractal escape-time map.
    """
    from mtft.burning_ship import burning_ship_array

    img = burning_ship_array(re_range, im_range, resolution, max_iter)
    data = {
        "image": img,
        "re_range": re_range,
        "im_range": im_range,
        "resolution": resolution,
    }

    if plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 7))
        ax.imshow(
            img,
            extent=[re_range[0], re_range[1], im_range[0], im_range[1]],
            cmap="inferno",
            origin="lower",
            aspect="auto",
        )
        ax.set_xlabel("Re(c)")
        ax.set_ylabel("Im(c)")
        ax.set_title("Burning Ship Fractal — MTFT Fermion Vacuum")
        plt.tight_layout()
        return {**data, "fig": fig}

    return data


# ═══════════════════════════════════════════════════════════════
#  Verification Scorecard (bar chart)
# ═══════════════════════════════════════════════════════════════

def scorecard_plot(plot: bool = True) -> dict:
    """
    Bar chart of all MTFT predictions vs experiment.
    """
    from mtft.verify import full_report

    report = full_report()
    names = [p.name for p in report]
    errors = [min(p.error_percent, 50) for p in report]  # cap for display
    colors = ["green" if e < 1 else "orange" if e < 5 else "red" for e in errors]

    data = {"names": names, "errors": errors, "colors": colors}

    if plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(12, 7))
        y_pos = range(len(names))
        ax.barh(y_pos, errors, color=colors, edgecolor="gray", linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel("Error (%)")
        ax.set_title("MTFT Verification Scorecard — All Predictions vs PDG")
        ax.axvline(1.0, color="green", ls="--", alpha=0.5, label="1% threshold")
        ax.axvline(5.0, color="orange", ls="--", alpha=0.5, label="5% threshold")
        ax.legend()
        ax.invert_yaxis()
        plt.tight_layout()
        return {**data, "fig": fig}

    return data


# ═══════════════════════════════════════════════════════════════
#  Data Export (for any visualization library)
# ═══════════════════════════════════════════════════════════════

def export_for_plotly(data: dict) -> dict:
    """
    Convert any mtft.viz data dict to plotly-friendly format.

    NumPy arrays → lists, complex → (real, imag) tuples.
    """
    out = {}
    for k, v in data.items():
        if isinstance(v, np.ndarray):
            if np.iscomplexobj(v):
                out[k] = {"real": v.real.tolist(), "imag": v.imag.tolist()}
            else:
                out[k] = v.tolist()
        elif k == "fig":
            continue  # skip matplotlib figures
        else:
            out[k] = v
    return out
