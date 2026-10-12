"""TL-01 (2026-10-07): the gamma tilt is a rotation — it cancels between up and down quarks, and where it survives, V_ub measures the down quarks' share of V_cb.

1. The tilt is a rotation (EXACT).  Over each curve point the gamma bundle restricts to the flat torus bundle N_{c phi(x)} (EPS-01), the pattern of a Poincare
   pullback: its curvature mixes the curve and torus directions (odd under w -> -w, the parity breaking DD-01 needs).  On the jets at the jet point a mixed
   curvature makes the x- and w-jets non-orthogonal, so the physical jets are a rotated pair: rotation-type breaking, LV-01's pure lever.
2. Rotations of matter blocks cannot make the Cabibbo angle (EXACT; checked to 1e-17).  The quark doublet block is shared by up and down, so its tilt cancels in
   U_u^dagger U_d; right-handed tilts drop out of Y Y^dagger; with the parity-decoupled vertex V_us = V_ub = 0.  The tilt must act inside the vertex and differently
   for up and down — EPS-01's eps_d / eps_u = 2.3 split of the condensate's SU(2)_R components.
3. Where the tilt survives, V_ub measures the down quarks' share of V_cb (EXACT at leading order; exact diagonalisation agrees to 3%, the residual being cos u_d).  Sector-dependent
   aligned tilts give V_us = |u_d - u_u|, V_ub = V_us |s23_d|, hence c = V_ub / (V_us V_cb) = |s23_d| / V_cb.
4. The programme's numbers (NUMERICAL).  FL-01's shared right-handed tail and CD-01's dressing give s23_u / s23_d = (m_c/m_t) / (z m_s/m_b) = 0.19 / z, so
   c >= 1 / (1 + 0.19/z).  At DP-01's z = 0.25, c >= 0.565 — outside the pre-registered window [0.35, 0.45]: the test fails.  The window needs z <= 0.157 with nearly
   opposite 2–3 phases, the up quarks supplying at least 55% of V_cb.
5. Verdict (STRUCTURAL).  The data's c = 0.40 needs either LV-01's local arm (opposite arms with a twist) or stronger up-sector dressing (z <= 0.16, beyond what
   DP-01's towers already strain to give).  Both need the condensate's SU(2)_R split at the jet point — now the one computation every route points to.
CC-58 caveats CA-1 (the r-window reads V_us as a block tilt, which cancels; it maps the tilt's size only to order of magnitude)."""
import numpy as np

MC_MT, MS_MB = 0.0037, 0.01924
Z_DP01 = 0.25
WINDOW = (0.35, 0.45)


def _rot12(t):
    c, s = np.cos(t), np.sin(t); return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], complex)

def sector(u, b, bp, m11, phase23=0.0):
    """Parity-symmetric vertex (w decoupled up to a small (1,1) mass) tilted on the left by angle u."""
    return _rot12(u) @ np.array([[m11, 0, 0], [0, 0, b * np.exp(1j * phase23)], [0, bp, 1.0]], complex)

def ckm(Yu, Yd):
    Uu = np.linalg.svd(Yu)[0][:, ::-1]; Ud = np.linalg.svd(Yd)[0][:, ::-1]; return Uu.conj().T @ Ud

def left_angle23(Y):
    return abs(np.linalg.svd(Y)[0][:, ::-1][1, 2])

def block_tilts(t_q=0.2, t_u=0.1, t_d=0.3):
    """Matter-block rotations: shared Q tilt on the left, separate right-handed tilts. Returns |V_us|, |V_ub|."""
    Y0u = np.array([[2e-5, 0, 0], [0, 0, 0.02], [0, 0.3, 1.0]], complex); Y0d = np.array([[1e-3, 0, 0], [0, 0, 0.04], [0, 0.45, 1.0]], complex)
    Yu = _rot12(t_q).T @ Y0u @ _rot12(t_u); Yd = _rot12(t_q).T @ Y0d @ _rot12(t_d)
    V = ckm(Yu, Yd); return abs(V[0, 1]), abs(V[0, 2])

def down_share_check(ud=0.25, uu=0.03, bd=0.016, bu=0.025, phase=np.pi):
    Yd = sector(ud, bd, 0.45, 1e-3); Yu = sector(uu, bu, 0.3, 5e-5, phase); V = ckm(Yu, Yd)
    return abs(V[0, 2]) / (abs(V[0, 1]) * abs(V[1, 2])), left_angle23(Yd) / abs(V[1, 2])

def ratio_up_down(z):
    return MC_MT / (z * MS_MB)

def c_range(z):
    R = ratio_up_down(z); return 1 / (1 + R), (1 / abs(1 - R) if abs(1 - R) > 1e-12 else float("inf"))

def z_for_window(c_max=WINDOW[1]):
    return MC_MT / (MS_MB * (1 / c_max - 1))
