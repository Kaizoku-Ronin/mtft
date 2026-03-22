"""
Dark Sector and Cosmology
=========================

MTFT replaces dark matter particles with tau-vortex halos:
  rho(r) = A^2 / (2 * r^2)

This gives flat rotation curves v_inf^2 = 2*pi*G*A^2 with
zero new particles. The Tully-Fisher relation emerges from
black hole entropy, not phenomenological fitting.

  pip install mtft
  python examples/05_dark_sector.py
"""

import numpy as np
from mtft.dark_sector import TauVortexHalo, rotation_curve, tully_fisher
from mtft.cosmology import FriedmannMTFT

# ── 1. Tau-vortex halo ───────────────────────────────────────
print("TAU-VORTEX DARK HALO")
print("=" * 55)
halo = TauVortexHalo(A=1e10)

print(f"  Vortex amplitude A = {halo.A:.2e}")
print(f"  Density profile: rho(r) = A^2 / (2 r^2)")
print()
print(f"  {'Radius':>10s}  {'rho(r)':>12s}")
print(f"  {'-'*24}")
for r in [0.1, 1.0, 10.0, 50.0, 100.0]:
    rho = halo.rho(r)
    print(f"  {r:10.2f}  {rho:12.2e}")

print(f"\n  Key: rho ~ 1/r^2 means M(r) ~ r  -> v(r) = const")
print(f"  No NFW profile. No dark matter particle. Just topology.")

# ── 2. Flat rotation curve ───────────────────────────────────
print(f"\nROTATION CURVE")
print("-" * 55)
r_arr = np.linspace(0.5, 50, 10)
v_bary_arr = np.sqrt(0.5 / r_arr)  # Keplerian
v_halo = 1.0  # constant for tau-vortex
for i, r in enumerate(r_arr):
    v_b = float(v_bary_arr[i])
    v_t = float(np.sqrt(v_b**2 + v_halo**2))
    print(f"  r={r:5.1f}  v_bary={v_b:.3f}  v_halo={v_halo:.3f}  v_total={v_t:.3f}")

print(f"\n  v_total flattens -> this IS the observed galaxy rotation curve")

# ── 3. Tully-Fisher relation ─────────────────────────────────
print(f"\nTULLY-FISHER: L ~ v^4 (from BH entropy)")
print("-" * 55)
for v_inf in [100, 150, 200, 250, 300]:
    L = tully_fisher(float(v_inf))
    print(f"  v_inf = {v_inf} km/s  ->  L ~ {L:.2e}")

print(f"\n  The exponent 4 comes from black hole entropy S ~ A ~ M^2,")
print(f"  not from baryonic physics fitting.")

# ── 4. Modified Friedmann cosmology ──────────────────────────
print(f"\nMODIFIED FRIEDMANN EQUATION WITH tau-FIELD")
print("=" * 55)
cosmo = FriedmannMTFT(
    H0=70.0,           # km/s/Mpc
    Omega_m=0.3,
    Omega_r=9e-5,
    Omega_tau=0.01,     # tau-field energy fraction
)

print(f"  H_0         = {cosmo.H0} km/s/Mpc")
print(f"  Omega_m     = {cosmo.Omega_m}")
print(f"  Omega_tau   = {cosmo.Omega_tau}")
print(f"  Omega_total = {cosmo.Omega_m + cosmo.Omega_r + cosmo.Omega_tau + cosmo.Omega_Lambda:.4f}")
print()

history = cosmo.expansion_history(a_min=0.01, a_max=1.0, n_points=200)
print(f"  {'Scale a':>8s}  {'H(a)/H0':>10s}  {'Redshift z':>10s}")
print(f"  {'-'*30}")
for i in range(0, len(history["a"]), 40):
    a = history["a"][i]
    H = history["H_over_H0"][i]
    z = 1.0 / a - 1
    print(f"  {a:8.3f}  {H:10.2f}  {z:10.1f}")

print(f"\n  The tau-field adds a small oscillatory component to H(a),")
print(f"  potentially detectable in BAO data.")
