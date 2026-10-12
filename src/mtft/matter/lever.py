"""LV-01 (2026-10-07): the lever at the jet point — what V_ub measures, and a correction to FH-4.

The first family (the w-jet at P1) can reach the rest of the down sector by two routes once parity is broken there.
  the lever arm u = Y13 / Y23: its grip on the HEIGHT of the top's tail at P1, relative to the x-family's grip on the tail's SLOPE
                               (a misaligned D-term component times the tail's value-to-slope ratio);
  the local arm v = Y12 / m2:  its direct coupling to the x-family at P1, relative to the second family's mass.
The jet point is the fulcrum; the tail is the long arm.

1. The lever formulas (EXACT at leading order; checked against exact diagonalisation to ~10%).  With texture zeros at (1,1) and (2,2) and equal left and
   right arms, the down sector gives V_us = |u - v|, V_ub / V_cb = |u|, m_d / m_s = |v (2u - v)|, hence c = V_ub / (V_us V_cb) = |u| / |u - v|: c is the long
   arm's share of the Cabibbo angle.  A pure lever gives c = 1 and a massless first family; a pure local arm gives GST (m_d/m_s = V_us^2) and leaves V_ub to the
   up sector, c = theta_12^u / V_us = 0.20.
2. A correction (EXACT).  FH-4 attached the texture-zero relation |V_ub/V_cb| = theta_12^u to rotation-type breaking.  Rotation-type breaking is the pure lever
   (Y13 = eps Y23): c = 1, first family massless.  The texture-zero relation belongs to the pure local arm.  CC-57 supersedes FH-4.
3. What the data ask (NUMERICAL).  Down sector alone, c = 0.40 needs the arms opposite: u = 0.090, v = -0.135, which gives m_d/m_s = 0.042 against 0.050.
4. The CP test (ESTIMATE: up sector with its symmetric 1–2 angle sqrt(m_u/m_c) and one phase).  Fitting V_us, V_cb, V_ub, m_d/m_s and m_s/m_b exactly with real
   arms predicts the unitarity angle gamma = 17 deg (opposite arms) or 141 deg (same-sign), against the measured 65.7 deg: real arms fail.  A relative twist between
   the arms fits it — beta = -130 deg (lever share 0.26) or -89 deg (lever share 0.46).  The right-angle branch mirrors the measured CKM's near-right unitarity
   angle alpha rather than anything specific to the mechanism.
5. Where the test goes (STRUCTURAL).  c and gamma are set by the condensate's jets at the jet point — its first-order misalignment, its mixed second-order
   component and their relative phase — which no wave has computed.  The pre-registered window c in [0.35, 0.45] moves there, with gamma as a second target.
Provenance: fits by scipy.optimize.least_squares in refit(); the numbers below are frozen from it."""
import numpy as np

V_US, V_CB, V_UB = 0.2250, 0.0410, 0.0037
MD_MS, MS_MB, MU_MC, MC_MT = 0.0503, 0.01924, 0.0020, 0.0037
GAMMA_MEASURED = 65.7
REAL_ARM_FITS = ({"u": 0.1258, "v": -0.1351, "gamma": 16.8, "J": 9.62e-6}, {"u": 0.0364, "v": 0.2656, "gamma": 141.0, "J": 2.09e-5})
TWIST_FITS = ({"u": 0.0590, "v_abs": 0.1846, "beta": -129.7, "lever_share": 0.262}, {"u": 0.1038, "v_abs": 0.1827, "beta": -88.6, "lever_share": 0.461})


def y_down(b, bp, u, v, up=None, vp=None):
    up = u if up is None else up; vp = v if vp is None else vp; m2 = b * bp
    return np.array([[0, v * m2, u * b], [vp * m2, 0, b], [up * bp, bp, 1.0]], complex)

def down_only(u, v, b=0.02, bp=0.03):
    """Exact diagonalisation with the up sector diagonal: returns V_us, V_ub/V_cb, m1/m2, c."""
    U, s, _ = np.linalg.svd(y_down(b, bp, u, v)); U, s = U[:, ::-1], s[::-1]
    return abs(U[0, 1]), abs(U[0, 2]) / abs(U[1, 2]), s[0] / s[1], abs(U[0, 2]) / (abs(U[0, 1]) * abs(U[1, 2]))

def formulas(u, v):
    return abs(u - v), abs(u), abs(v * (2 * u - v)), abs(u) / abs(u - v)

def opposite_arms(c=V_UB / (V_US * V_CB), vus=V_US):
    u = c * vus; v = u - vus; return u, v, abs(v * (2 * u - v))

def up_sector_c():
    return np.sqrt(MU_MC) / V_US

def refit(seed=1, tries=300):
    """Re-run the real-arm fit (needs scipy); returns the distinct solutions as (u, v, gamma)."""
    from scipy.optimize import least_squares
    def ckm(phi, b, bp, u, v):
        A = np.exp(1j * phi) * np.sqrt(MU_MC) * MC_MT
        Yu = np.array([[0, A, 0], [A, MC_MT, 0], [0, 0, 1.0]], complex)
        Uu, su, _ = np.linalg.svd(Yu); Ud, sd, _ = np.linalg.svd(y_down(b, bp, u, v))
        return Uu[:, ::-1].conj().T @ Ud[:, ::-1], sd[::-1]
    def res(p):
        b, bp, u, v, phi = p; V, sd = ckm(phi, b, bp, u, v)
        return [abs(V[0, 1]) / V_US - 1, abs(V[1, 2]) / V_CB - 1, abs(V[0, 2]) / V_UB - 1, sd[0] / sd[1] / MD_MS - 1, sd[1] / sd[2] / MS_MB - 1]
    rng = np.random.default_rng(seed); out = set()
    for _ in range(tries):
        r = least_squares(res, [0.041 * rng.uniform(0.5, 1.5), 0.47 * rng.uniform(0.5, 1.5), rng.uniform(-0.4, 0.4), rng.uniform(-0.4, 0.4), rng.uniform(0, 6.28)],
                          bounds=([0.005, 0.05, -1, -1, -10], [0.2, 2, 1, 1, 10]))
        if np.max(np.abs(r.fun)) < 1e-8:
            b, bp, u, v, phi = r.x; V, _ = ckm(phi, b, bp, u, v)
            g = np.degrees(np.angle(-V[0, 0] * np.conj(V[0, 2]) / (V[1, 0] * np.conj(V[1, 2]))))
            out.add((round(abs(u), 3), round(np.sign(u * v), 0), round(abs(g), 0)))
    return sorted(out)
