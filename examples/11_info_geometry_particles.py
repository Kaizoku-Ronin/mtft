"""
Information Geometry and the Standard Model Catalog
====================================================

Information geometry connects the arithmetic measure on the
modular parameter to physical observables:
  - The Fisher-Rao metric measures distinguishability
  - Its curvature diverges at the Feigenbaum point (onset of chaos)
  - Zamolodchikov's c-theorem constrains RG flow monotonicity

The particle catalog encodes every SM particle with its
MTFT kappa-coupling, mass prediction, and hierarchy level.

  pip install mtft
  python examples/11_info_geometry_particles.py
"""

import numpy as np
from mtft.info_geometry import (
    fisher_rao_metric, ricci_scalar_logistic,
    lyapunov_exponent, ZamolodchikovFlow,
)
from mtft.particles import Particle, StandardModel
from mtft.hosotani import HosotaniPotential
from mtft.constants import FEIGENBAUM_DELTA, PI

# ── 1. Fisher-Rao metric ─────────────────────────────────────
print("FISHER-RAO METRIC ON THE MODULAR MEASURE")
print("=" * 55)
print(f"  The metric g(r) measures how fast the logistic map's")
print(f"  invariant measure changes with parameter r.")
print()
print(f"  {'r':>8s}  {'g(r)':>12s}")
print(f"  {'-'*22}")
for r in [2.5, 3.0, 3.2, 3.4, 3.5, 3.56, 3.57, 3.8]:
    g = fisher_rao_metric(r)
    print(f"  {r:8.2f}  {g:12.4f}")

print(f"\n  Peaks near r_inf = 3.5699 — the Feigenbaum accumulation")
print(f"  point. Maximum information content at onset of chaos.")

# ── 2. Logistic map bridge ───────────────────────────────────
print(f"\nLOGISTIC MAP AT THE EDGE OF CHAOS")
print("-" * 55)
r_values = [2.5, 3.0, 3.449, 3.544, 3.5699, 3.8, 4.0]
print(f"  {'r':>8s}  {'Lyapunov':>10s}  {'Ricci':>10s}  {'Status':>12s}")
for r in r_values:
    lam = lyapunov_exponent(r, n_iter=5000)
    R = ricci_scalar_logistic(r)
    status = "periodic" if lam < 0 else "edge" if abs(lam) < 0.01 else "chaotic"
    print(f"  {r:8.4f}  {lam:10.4f}  {R:10.2f}  {status:>12s}")

print(f"\n  The Feigenbaum point r_inf = 3.5699... is where the")
print(f"  Ricci curvature diverges — the onset of chaos.")
print(f"  delta = {FEIGENBAUM_DELTA:.10f} controls the approach.")

# ── 3. Zamolodchikov c-theorem (RG flow) ─────────────────────
print(f"\nZAMOLODCHIKOV c-THEOREM (RG flow monotonicity)")
print("-" * 55)
flow = ZamolodchikovFlow(kappa=1.0)
beta_val = flow.beta(1.0)
print(f"  beta function at g=1.0: {beta_val:.6f}")
print(f"  The c-theorem guarantees monotonic decrease along RG flow.")
print(f"  MTFT arithmetic weights satisfy this — the stiffness mu_N")
print(f"  decreases monotonically with increasing modular depth y.")

# ── 4. Hosotani mechanism ────────────────────────────────────
print(f"\nHOSOTANI GAUGE-HIGGS UNIFICATION")
print("=" * 55)
hp = HosotaniPotential(fermion_fraction=0.4)
theta0 = hp.find_vacuum()
print(f"  Fermion fraction: {hp.fermion_fraction}")
print(f"  Vacuum angle:     theta_0 = {theta0:.6f}")
print(f"  V(theta_0):       {hp.V_bosonic(theta0) + hp.V_fermionic(theta0):.6f}")
print(f"  V''(theta_0):     {hp.second_derivative(theta0):.4f} (> 0 = stable)")
print(f"\n  The Higgs is the 5th component of the gauge field.")
print(f"  EWSB = choosing theta_0 != 0 in the holonomy potential.")

# ── 5. Standard Model particle catalog ───────────────────────
print(f"\nSTANDARD MODEL PARTICLE CATALOG")
print("=" * 55)
sm = StandardModel()
particles = sm.all
print(f"  {len(particles)} particles loaded\n")

print(f"  {'Name':15s} {'Mass (GeV)':>12s} {'kappa':>8s} {'Type':>10s}")
print(f"  {'-'*48}")
for p in sorted(particles, key=lambda p: p.mass_GeV):
    print(f"  {p.name:15s} {p.mass_GeV:12.4f} {p.kappa:8.4f} {p.ptype.value:>10s}")

# ── 6. Mass hierarchy from kappa ─────────────────────────────
print(f"\nMASS HIERARCHY: kappa-couplings span 12 orders of magnitude")
print("-" * 55)
massive = [p for p in particles if p.mass_GeV > 0]
lightest = min(massive, key=lambda p: p.mass_GeV)
heaviest = max(massive, key=lambda p: p.mass_GeV)
ratio = heaviest.mass_GeV / lightest.mass_GeV
print(f"  Lightest: {lightest.name} ({lightest.mass_GeV:.6f} GeV)")
print(f"  Heaviest: {heaviest.name} ({heaviest.mass_GeV:.2f} GeV)")
print(f"  Ratio:    {ratio:.0f}")
print(f"\n  MTFT explains this hierarchy through the arithmetic")
print(f"  weights — highly composite n have large w_n, giving")
print(f"  strong coupling to the vacuum.")
