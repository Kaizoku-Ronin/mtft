"""
Modular Forms: eta, theta, and the Spectral Determinant
========================================================

The Dedekind eta function is the backbone of MTFT:
  eta(tau) = q^(1/24) * prod_{n>=1} (1 - q^n)

It connects Ramanujan's sum 1+2+3+... = -1/12 to particle masses
through the dimensional bridge formula:
  m_e / m_P = |eta(tau_c)|^{2*alpha^-1/pi}

  pip install mtft
  python examples/06_modular_forms.py
"""

import math
import numpy as np
from mtft.forms import dedekind_eta, jacobi_theta3, eisenstein_E2k, j_invariant
from mtft.modular import TauField, sl2z_transform
from mtft.dimensional_bridge import (
    electron_mass_from_eta,
    charge_from_feigenbaum,
    alpha_F_from_impedance,
    feigenbaum_product_lock,
    modular_impedance,
)
from mtft.constants import PI, CriticalDepths

# ── 1. The tau-field and SL(2,Z) ─────────────────────────────
print("THE MODULAR PARAMETER tau IN THE UPPER HALF-PLANE")
print("=" * 55)
tau = TauField(0.0 + 0.5j)
print(f"  tau = 0.5i (purely imaginary)")
print()

# SL(2,Z) transformations
S = np.array([[0, -1], [1, 0]])    # tau -> -1/tau
T = np.array([[1, 1], [0, 1]])     # tau -> tau + 1

tau_val = 0.3 + 0.5j
tau_S = sl2z_transform(tau_val, S)
tau_T = sl2z_transform(tau_val, T)
print(f"  tau       = {tau_val}")
print(f"  S(tau)    = -1/tau = {tau_S:.4f}")
print(f"  T(tau)    = tau+1  = {tau_T:.4f}")
print()

# ── 2. Dedekind eta at critical depths ───────────────────────
print("DEDEKIND ETA AT MTFT CRITICAL DEPTHS")
print("-" * 55)
for label, y in [
    ("y_s1 (skeleton 1)", CriticalDepths.y_s1),
    ("y_conf (confinement)", CriticalDepths.y_conf),
    ("y_s2 (skeleton 2)", CriticalDepths.y_s2),
    ("y_spec (EM extraction)", CriticalDepths.y_spec),
]:
    tau_c = 1j * y
    eta = dedekind_eta(tau_c)
    print(f"  {label:25s}:  |eta(i*{y:.4f})| = {abs(eta):.8f}")

# ── 3. Modular form zoo ──────────────────────────────────────
print(f"\nMODULAR FORM EVALUATIONS AT tau = i*0.5")
print("-" * 55)
tau_test = 0.5j
print(f"  eta(tau)      = {dedekind_eta(tau_test):.8f}")
print(f"  theta_3(0,tau) = {jacobi_theta3(0, tau_test):.8f}")
print(f"  E_2(tau)      = {eisenstein_E2k(tau_test, k=1):.6f}")
print(f"  E_4(tau)      = {eisenstein_E2k(tau_test, k=2):.6f}")
print(f"  E_6(tau)      = {eisenstein_E2k(tau_test, k=3):.6f}")
print(f"  j(tau)        = {j_invariant(tau_test):.2f}")
print(f"  j(i) (exact)  = {j_invariant(1j):.0f}  (should be 1728)")

# ── 4. The dimensional bridge ────────────────────────────────
print(f"\nTHE DIMENSIONAL BRIDGE: m_e from eta(tau)")
print("=" * 55)
result = electron_mass_from_eta()
print(f"  Formula: m_e/m_P = |eta(tau_c)|^(2*alpha^-1/pi)")
print(f"  y_c = {result['y_c']:.5f}")
print(f"  |eta(tau_c)| = {result['abs_eta']:.8f}")
print(f"  Exponent 2*alpha^-1/pi = {result['exponent_2alpha_inv_over_pi']:.4f}")
print(f"  ln(m_P/m_e) predicted = {result['ln_mP_over_me_predicted']:.4f}")
print(f"  ln(m_P/m_e) measured  = {result['ln_mP_over_me_measured']:.4f}")
print(f"  Error: {result['error_percent']:.2f}%")

# ── 5. Feigenbaum-EM identities ──────────────────────────────
print(f"\nFEIGENBAUM-ELECTROMAGNETISM IDENTITIES")
print("-" * 55)

ch = charge_from_feigenbaum()
print(f"  (i)   e = sqrt(2)/delta:")
print(f"        {ch['e_feigenbaum']:.6f} vs {ch['e_exact']:.6f}  ({ch['error_percent']:.3f}%)")

af = alpha_F_from_impedance()
print(f"  (ii)  alpha_F = exp(40*pi*alpha):")
print(f"        {af['alpha_F_predicted']:.6f} vs {af['alpha_F_actual']:.6f}  ({af['error_percent']:.3f}%)")

pl = feigenbaum_product_lock()
print(f"  (iii) delta^2 * ln(alpha_F) = 20:")
print(f"        {pl['delta2_ln_alphaF']:.6f}  ({pl['error_percent']:.4f}%)")

mi = modular_impedance()
print(f"  (iv)  (delta/alpha_F)^2 = 7/2:")
print(f"        {mi['ratio_squared']:.6f} vs 3.5  ({mi['error_percent']:.3f}%)")
