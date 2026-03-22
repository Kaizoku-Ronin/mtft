"""
Radioactive Decay in Modular Time
===================================

MTFT predicts corrections to nuclear decay rates from the
tau-field gradient. The correction is tiny in lab conditions
but potentially observable near neutron stars.

The neutron lifetime anomaly (bottle vs beam measurements)
may be an MTFT observable.

  pip install mtft
  python examples/13_decay.py
"""

import math
from mtft.decay import ModularDecay, u238_example, neutron_beta_decay

# ── 1. The modular decay formula ─────────────────────────────
print("RADIOACTIVE DECAY IN MODULAR TIME")
print("=" * 55)
print("  Standard: Gamma = Gamma_0")
print("  MTFT:     Gamma = Gamma_0 / [1 + C * delta^2_eff / Theta(tau*)]")
print()
print("  where Theta(tau*) = Im(tau*) * exp(-2*pi*Im(tau*))")
print("  controls suppression at large modular depth.")
print()

# ── 2. Decay rate vs modular depth ───────────────────────────
print("DECAY CORRECTION vs TAU-FIELD DEPTH")
print("-" * 55)
Gamma_0 = 1e-10  # arbitrary base rate
print(f"  {'Im(tau*)':>10s}  {'Theta':>12s}  {'Correction':>12s}  {'T_half ratio':>14s}")
print(f"  {'-'*50}")
for y in [0.01, 0.05, 0.10, 0.15, 0.18, 0.25, 0.40]:
    d = ModularDecay(Gamma_0=Gamma_0, delta_eff=0.001, tau_star=1j*y)
    print(f"  {y:10.2f}  {d.theta_tau:12.4e}  {d.half_life_correction:12.4f}  "
          f"{'~1 (lab)' if d.half_life_correction < 1.001 else f'{d.half_life_correction:.4f}'}")

# ── 3. U-238 alpha decay ────────────────────────────────────
print(f"\nU-238 ALPHA DECAY (T_1/2 = 4.468 Gyr)")
print("-" * 55)
u = u238_example()
print(f"  Laboratory:")
print(f"    T_half = {u['lab']['T_half_s']:.4e} s")
print(f"    Correction factor: {u['lab']['correction_factor']:.6f}")
print(f"  Neutron star surface:")
print(f"    T_half = {u['neutron_star']['T_half_s']:.4e} s")
print(f"    Correction factor: {u['neutron_star']['correction_factor']:.4f}")
print(f"  Prediction: {u['prediction']}")

# ── 4. Neutron lifetime anomaly ──────────────────────────────
print(f"\nNEUTRON LIFETIME ANOMALY")
print("=" * 55)
n = neutron_beta_decay()
print(f"  Bottle measurement: T_1/2 = {n['T_half_bottle_s']:.2f} s")
print(f"  Beam measurement:   T_1/2 = {n['T_half_beam_s']:.1f} s")
print(f"  Discrepancy:        {n['discrepancy_percent']:.2f}%  (~4 sigma)")
print(f"  Q-value:            {n['Q_value_GeV']*1000:.1f} MeV")
print(f"\n  MTFT interpretation:")
print(f"  {n['mtft_interpretation']}")
