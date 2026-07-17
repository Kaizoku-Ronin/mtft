"""
Lepton Masses: Koide Theorem and the Tano Mass Formula
=======================================================

Two independent predictions of the charged lepton masses:

1. Koide formula K = 2/3 — predicts m_tau to 0.006%
2. Tano formula theta_0 = 2/delta^2 + 4*arg(a2_cx(f3))
   — predicts m_mu and m_tau from X_0(143) with zero free
   parameters beyond the electron mass anchor.

Together they show that lepton masses are determined by:
  - The Koide constraint (fixes the angle in mass-root space)
  - The Hecke eigenvalue of X_0(143) (fixes the phase)
  - The electron mass (sets the scale)
"""

import math
from mtft.koide import (
    koide_ratio, koide_leptons, koide_up_quarks, koide_down_quarks,
    koide_angle_leptons, predict_tau_mass, bisector_stability_K,
    koide_manifold_point, MAGIC_ANGLE,
)
from mtft.x0_143 import (
    tano_mass_predictions, koide_angle_tano, koide_angle_experimental,
    generation_count, A2_COMPLEX, ORBIT_DIMENSIONS, GENUS, INDEX,
    JACOBIAN, verify_complex_eigenvalue,
)
from mtft.constants import LEPTONS, FEIGENBAUM_DELTA

# ── 1. The Koide formula ─────────────────────────────────────
print("THE KOIDE FORMULA: K = (sum m_i) / (sum sqrt(m_i))^2 = 2/3")
print("=" * 60)
print(f"  K(e, mu, tau)   = {koide_leptons():.6f}  (exact: 0.666667)")
print(f"  K(u, c, t)      = {koide_up_quarks():.4f}  (conjectured: 2/3)")
print(f"  K(d, s, b)      = {koide_down_quarks():.4f}  (conjectured: 2/3)")
print()

# ── 2. Why K = 2/3: bisector stability ───────────────────────
print("WHY 2/3?  The Burning Ship fold geometry:")
print(f"  K = 1 / (3 * cos^2(45 deg)) = 1/(3 * 1/2) = {bisector_stability_K():.10f}")
print(f"  The 45 deg bisector is the unique RG fixed point")
print(f"  in the fold quadrant with delta_x != delta_y.")
print()

# ── 3. The Koide phase and magic angle ───────────────────────
print("THE KOIDE MANIFOLD K = S^1 x R+")
print(f"  Magic angle: arccos(1/sqrt(3)) = {math.degrees(MAGIC_ANGLE):.2f} deg")
print(f"  Lepton phase: {math.degrees(koide_angle_leptons()):.4f} deg")
print()

# Check: every point on the Koide manifold has K = 2/3
print("  Verification (K on manifold at various phases):")
for phi in [0.2, 0.5, 1.0, 2.0]:
    w = koide_manifold_point(phi, r=2.0)
    if all(w > 0):
        K = koide_ratio(w[0]**2, w[1]**2, w[2]**2)
        print(f"    phi={phi:.1f}, r=2.0:  K = {K:.12f}")
print()

# ── 4. Tau mass prediction from Koide ────────────────────────
print("TAU MASS FROM KOIDE (0.006% precision)")
print("-" * 60)
result = predict_tau_mass()
print(f"  Input:  m_e  = 0.51099895 MeV")
print(f"  Input:  m_mu = 105.6583755 MeV")
print(f"  Output: m_tau = {result['predicted_MeV']:.3f} MeV")
print(f"  PDG:    m_tau = {result['experimental_MeV']:.2f} MeV")
print(f"  Error:  {result['error_percent']:.4f}%")
print()

# ── 5. X_0(143) — the arithmetic stage ──────────────────────
print("THE MODULAR CURVE X_0(143)")
print("=" * 60)
print(f"  Level N = 143 = 11 x 13")
print(f"  Genus = {GENUS}")
print(f"  Index = {INDEX} = |PSL(2,7)| = dim SU(13)")
print(f"  Galois orbit dimensions: {ORBIT_DIMENSIONS}")
print(f"  Fermion generations: {generation_count()} (theorem, not assumption)")
print()

# ── 6. The "complex Hecke eigenvalue" — RETIRED (audit B11) ──────────────
print("THE PHANTOM EIGENVALUE (SUPERSEDED)")
print("-" * 60)
print(f"  Historical value a2_cx(f3) = {A2_COMPLEX}")
print(f"  This was believed to be the ONLY non-real datum in all of X_0(143).")
print(f"  AUDIT RESULT: it is a root of the WRONG polynomial.  The correct")
print(f"  T_2|f_3 charpoly x^6 - 10x^4 + 2x^3 + 24x^2 - 7x - 12 has 6 REAL")
print(f"  roots (trivial nebentypus => all Hecke eigenvalues are totally real).")

check = verify_complex_eigenvalue()
print(f"  all roots real:        {check['all_roots_real']}")
print(f"  max |root|:            {check['max_abs_root']:.4f} (Ramanujan 2*sqrt(2) = {check['ramanujan_bound']:.4f})")
print(f"  phantom min distance to a true root: {check['phantom_min_distance_to_true_root']:.4f}")
print()

# ── 7. The Tano mass formula — SUPERSEDED (audit B10) ────────────────────
print("THE TANO MASS FORMULA (Paper 26, Theorem 8.1) — SUPERSEDED")
print("=" * 60)
print(f"  theta_0 = 2/delta^2 + dim(f2) * arg(a2_cx(f3))")
print(f"          = {math.degrees(koide_angle_tano()):.3f} deg  (historical value)")
print(f"  Experimental: {math.degrees(koide_angle_experimental()):.3f} deg")
print(f"  Difference:   {abs(math.degrees(koide_angle_tano()) - math.degrees(koide_angle_experimental())):.3f} deg")
print(f"  NOTE: arg(a2_cx) uses the phantom value, so this is a numerical")
print(f"  coincidence, not arithmetic of X_0(143).  Retained for continuity.")
print()

pred = tano_mass_predictions()
print(f"  {'Particle':12s} {'MTFT (MeV)':>12s} {'PDG (MeV)':>12s} {'Error':>10s}")
print(f"  {'-'*46}")
for p in ["electron", "muon", "tau"]:
    m = pred["masses_MeV"][p]
    e = pred["experimental_MeV"][p]
    err = pred["errors_percent"][p]
    tag = "(anchor)" if p == "electron" else f"{err:.2f}%"
    print(f"  {p:12s} {m:12.3f} {e:12.3f} {tag:>10s}")

# ── 8. Jacobian stiffness matrix ─────────────────────────────
print(f"\nJACOBIAN STIFFNESS on J_0(143)")
print("-" * 60)
print(f"  lambda_3/lambda_2 = {JACOBIAN.ratio_32:.2f}  (delta_x = 7.24, match: 2.3%)")
print(f"  lambda_2/lambda_1 = {JACOBIAN.ratio_21:.2f}  (delta_F = 4.67, match: 9.8%)")
print(f"  Electron-muon coupling: {JACOBIAN.coupling_12:.3f} (decoupled)")
print(f"  Muon-tau coupling:      {JACOBIAN.coupling_23:.2f} (strongly mixed)")
print(f"  AL -> mass basis rotation: {JACOBIAN.rotation_deg} deg")
