"""
The Burning Ship: Fermion Vacuum and Anisotropic Universality
=============================================================

The Burning Ship iteration z -> (|Re z| + i|Im z|)^2 + c
breaks analyticity through fold lines at x=0 and y=0.

MTFT maps:
  Mandelbrot (analytic)      -> gauge sector  -> Yang-Mills mass gap
  Burning Ship (non-analytic) -> fermion sector -> Koide formula

The fold lines create anisotropic Feigenbaum constants:
  delta_x ~ 7.30 (quark direction)
  delta_y ~ 3.38 (lepton direction)

The mass hierarchy between generations comes from:
  m_{n+1}/m_n ~ delta_eff^2

  pip install mtft
  python examples/07_burning_ship.py
"""

import math
from mtft.burning_ship import (
    burning_ship_iterate, burning_ship_array,
    sector, jacobian, jacobian_determinant, jacobian_trace,
    ANISOTROPIC, FOLD_PROXIMITY, satellite_mass_scale,
)
from mtft.constants import FEIGENBAUM_DELTA, DELTA_X, DELTA_Y

# ── 1. The iteration ─────────────────────────────────────────
print("BURNING SHIP: z -> (|Re z| + i|Im z|)^2 + c")
print("=" * 55)
test_points = [
    (0.0+0.0j, "origin (bounded)"),
    (-1.0+0.0j, "period-2 (bounded)"),
    (-1.756-0.028j, "deep ship (bounded)"),
    (2.0+0.0j, "outside (escapes)"),
]
for c, label in test_points:
    n = burning_ship_iterate(c, max_iter=200)
    status = "bounded" if n == 200 else f"escapes at n={n}"
    print(f"  c = {str(c):>14s}  {label:25s}  {status}")

# ── 2. Phase space sectors ───────────────────────────────────
print(f"\nPHASE SPACE SECTORS (sgn(x), sgn(y))")
print("-" * 55)
for z in [1+1j, -1+1j, 1-1j, -1-1j]:
    s = sector(z)
    J = jacobian(z)
    det = jacobian_determinant(z)
    tr = jacobian_trace(z)
    print(f"  z = {str(z):>8s}  sector = {s}  det(J) = {det:.1f}  tr(J) = {tr:.1f}")

print(f"\n  det(J) = 4|z|^2 is CONTINUOUS across folds")
print(f"  tr(J) is DISCONTINUOUS at x=0 (the fold line)")
print(f"  This discontinuity IS the origin of fermion mass splitting.")

# ── 3. Anisotropic Feigenbaum constants ──────────────────────
print(f"\nANISOTROPIC FEIGENBAUM UNIVERSALITY")
print("=" * 55)
A = ANISOTROPIC
print(f"  Standard delta (Mandelbrot): {A.delta:.10f}")
print(f"  delta_x (attractor fold):    {A.delta_x:.4f}  (quark direction)")
print(f"  delta_y (barrier fold):      {A.delta_y:.4f}  (lepton direction)")
print(f"  beta_x  (x-fold approach):   {A.beta_x:+.2f}")
print(f"  beta_y  (y-fold approach):   {A.beta_y:+.2f}")
print(f"  |beta_x| + |beta_y| = {A.fold_budget:.2f}  (should be 0.50)")
print()

print(f"  Generation mass ratios:")
print(f"    Leptons:  delta_y^2 = {A.generation_mass_ratio('lepton'):.1f}  (m_mu/m_e ~ 207)")
print(f"    Quarks:   delta_x^2 = {A.generation_mass_ratio('quark'):.1f}  (m_s/m_d ~ 20)")

# ── 4. Fold proximity data ───────────────────────────────────
print(f"\nFOLD PROXIMITY vs PERIOD (Paper 4, Table 2)")
print("-" * 55)
print(f"  {'Period':>6s}  {'min |x|':>10s}  {'min |y|':>10s}  {'x-cross':>8s}  {'y-cross':>8s}")
for period, (mx, my, xc, yc) in sorted(FOLD_PROXIMITY.items()):
    print(f"  {period:6d}  {mx:10.2e}  {my:10.2e}  {xc:8d}  {yc:8d}")

# ── 5. Three generations from satellites ─────────────────────
print(f"\nTHREE GENERATIONS FROM SELF-SIMILAR SATELLITES")
print("-" * 55)
for n in range(1, 6):
    scale = satellite_mass_scale(n, delta_eff=DELTA_Y)
    print(f"  Generation {n}: m_{n}/m_1 ~ {scale:.4e}")

print(f"\n  Generations 4+ fall below the EW scale -> not realized.")
print(f"  This is WHY there are exactly 3 generations of fermions.")
