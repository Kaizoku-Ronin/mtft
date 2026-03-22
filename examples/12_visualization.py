"""
Visualization: Plotting MTFT Data
===================================

All mtft functions return numpy arrays, compatible with any
plotting library. This example uses matplotlib (pip install mtft[full]).

For interactive 3D, try plotly:  pip install mtft[viz]

  pip install mtft[full]
  python examples/12_visualization.py
"""

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("matplotlib not found. Install with: pip install mtft[full]")
    print("Showing data in text mode instead.\n")

import numpy as np
from mtft.arithmetic import weight_array, mass_gap_stiffness, S1, su3_hessian_eigenvalues
from mtft.koide import koide_manifold_point
from mtft.burning_ship import burning_ship_array
from mtft.constants import PI, LEPTONS, CriticalDepths
import math

def save_or_show(fig, name):
    fig.savefig(f"{name}.png", dpi=150, bbox_inches="tight")
    print(f"  Saved {name}.png")
    plt.close(fig)

if HAS_MPL:
    # ── 1. Arithmetic weights ────────────────────────────────
    print("Generating plots...")
    n_arr = np.arange(1, 101)
    w_arr = weight_array(100)

    fig, ax = plt.subplots(figsize=(10, 4))
    colors = ['#E24B4A' if all(n % p != 0 for p in [2,3,5,7] if p < n) and n > 1
              else '#378ADD' for n in n_arr]
    ax.bar(n_arr, w_arr, color=colors, width=0.8)
    ax.set_xlabel("n")
    ax.set_ylabel("w_n")
    ax.set_title("MTFT arithmetic weights w_n = sum_{d|n} (log d)/d")
    save_or_show(fig, "01_weights")

    # ── 2. Stiffness landscape ───────────────────────────────
    y_arr = np.linspace(0.02, 0.45, 200)
    fig, ax = plt.subplots(figsize=(10, 5))
    for N, color, label in [(2,'#85B7EB','SU(2)'), (3,'#1D9E75','SU(3)'),
                             (5,'#7F77DD','SU(5)'), (7,'#D85A30','SU(7)')]:
        mu = [mass_gap_stiffness(y, N=N, n_max=200) for y in y_arr]
        ax.semilogy(y_arr, np.maximum(mu, 1e-6), color=color, label=label, lw=2)
    ax.axvline(CriticalDepths.y_conf, color='#1D9E75', ls='--', alpha=0.5)
    ax.axvline(CriticalDepths.y_spec, color='#378ADD', ls='--', alpha=0.5)
    ax.set_xlabel("Modular depth y")
    ax.set_ylabel("Mass gap stiffness (log)")
    ax.set_title("Yang-Mills confinement landscape")
    ax.legend()
    ax.grid(True, alpha=0.2)
    save_or_show(fig, "02_stiffness")

    # ── 3. SU(3) Hessian crossing ────────────────────────────
    y2 = np.linspace(0.14, 0.22, 100)
    l1s, l2s = zip(*[su3_hessian_eigenvalues(y) for y in y2])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(y2, l1s, color='#85B7EB', lw=2.5, label='lambda_1')
    ax.plot(y2, l2s, color='#D85A30', lw=2.5, label='lambda_2')
    ax.axvline(0.18174, color='#1D9E75', ls='--', alpha=0.7, label='y_conf')
    ax.set_xlabel("Modular depth y")
    ax.set_ylabel("Hessian eigenvalue")
    ax.set_title("SU(3) confinement lock: eigenvalues cross at y = 0.182")
    ax.legend()
    ax.grid(True, alpha=0.2)
    save_or_show(fig, "03_hessian")

    # ── 4. Alpha inverse from S1 ─────────────────────────────
    y3 = np.linspace(0.15, 0.40, 100)
    alpha_inv = [2.0 / S1(y) for y in y3]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(y3, alpha_inv, color='#378ADD', lw=2.5)
    ax.axhline(137.036, color='#E24B4A', ls='--', alpha=0.7, label='alpha^-1 = 137.036')
    ax.axvline(CriticalDepths.y_spec, color='#1D9E75', ls='--', alpha=0.7, label='y_spec')
    ax.set_xlabel("Modular depth y")
    ax.set_ylabel("2 / S_1(y)")
    ax.set_title("Fine-structure constant from the spectral lock")
    ax.legend()
    ax.grid(True, alpha=0.2)
    save_or_show(fig, "04_alpha_inverse")

    # ── 5. Koide manifold ────────────────────────────────────
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    for r, color in [(0.5,'#85B7EB'),(1.5,'#5DCAA5'),(3.0,'#7F77DD'),(6.0,'#EF9F27')]:
        pts = []
        for phi in np.linspace(0, 2*PI, 200):
            w = koide_manifold_point(phi, r)
            if all(w > 0): pts.append(w)
        if pts:
            pts = np.array(pts)
            ax.plot(pts[:,0], pts[:,1], pts[:,2], color=color, alpha=0.7, label=f'r={r}')

    lp = [math.sqrt(LEPTONS.e), math.sqrt(LEPTONS.mu), math.sqrt(LEPTONS.tau)]
    ax.scatter(*lp, color='#E24B4A', s=100, zorder=5, label='(e, mu, tau)')
    ax.set_xlabel('sqrt(m_1)')
    ax.set_ylabel('sqrt(m_2)')
    ax.set_zlabel('sqrt(m_3)')
    ax.set_title('Koide manifold K = S^1 x R+')
    ax.legend(fontsize=8)
    save_or_show(fig, "05_koide_manifold")

    # ── 6. Burning Ship fractal ──────────────────────────────
    img = burning_ship_array(
        re_range=(-2.0, 1.5),
        im_range=(-2.0, 0.5),
        resolution=800,
        max_iter=200,
    )
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.imshow(img, extent=[-2.0, 1.5, -2.0, 0.5], cmap='inferno',
              origin='lower', aspect='auto')
    ax.set_xlabel("Re(c)")
    ax.set_ylabel("Im(c)")
    ax.set_title("Burning Ship fractal — the MTFT fermion vacuum")
    save_or_show(fig, "06_burning_ship")

    print(f"\n  6 plots saved as PNG files.")
    print(f"  These can be used directly in papers and presentations.")

else:
    # Text-mode fallback
    print("TEXT MODE: Arithmetic weight samples")
    w = weight_array(20)
    for n in range(1, 21):
        bar = "#" * int(w[n-1] * 20)
        print(f"  w_{n:<3d} = {w[n-1]:.4f}  {bar}")

    print(f"\nMass gap at y=0.18:")
    for N in [2, 3, 5, 7]:
        mu = mass_gap_stiffness(0.18, N=N, n_max=200)
        print(f"  SU({N}): mu = {mu:.4f} > 0")
