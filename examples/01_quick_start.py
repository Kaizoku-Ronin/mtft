"""
MTFT Quick Start
================

Run this first:  pip install mtft
Then:            python examples/01_quick_start.py

This shows the core idea in 30 lines: fundamental constants
derived from the integers alone, with zero free parameters.
"""

import mtft

print(f"mtft v{mtft.__version__}\n")

# ── Every gauge coupling, from arithmetic ────────────────────
print("GAUGE SECTOR (derived from delta, T_inf, integers):")
print(f"  alpha^-1 = 2*pi*delta^2 + 1/(4*delta) - xi*T/delta^6")
print(f"           = {mtft.GAUGE.alpha_inv:.6f}  (CODATA: 137.035999)")
print(f"  alpha_s  = T_inf / 4          = {mtft.GAUGE.alpha_s:.4f}  (PDG: 0.1180)")
print(f"  sin^2 theta_W = 3/13          = {mtft.GAUGE.sin2_theta_W:.6f}  (PDG: 0.23122)")
print(f"  M_W/M_Z = 1/(2*Omega)         = {mtft.GAUGE.W_Z_ratio:.4f}  (PDG: 0.8814)")
print()

# ── The Higgs mass, from gamma and Omega ─────────────────────
print("HIGGS SECTOR:")
print(f"  m_H = v * gamma / (2*Omega)    = {mtft.HIGGS.m_H:.2f} GeV  (PDG: 125.25)")
print()

# ── The full scorecard ───────────────────────────────────────
from mtft.verify import full_report

report = full_report()
under_1 = sum(1 for p in report if p.error_percent < 1.0)
print(f"SCORECARD: {under_1} of {len(report)} predictions under 1% error")
print(f"           Zero free parameters. All from the integers.")
