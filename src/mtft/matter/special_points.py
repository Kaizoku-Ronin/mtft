"""SP-01 (2026-10-05): where the top can live — distances on X_0(143) in its canonical metric.

1. Distances from the blown-up cusp (NUMERICAL).  In the canonical (Arakelov) metric of area 1, from the Bergman density of the thirteen differentials (the cusp
   density 0.02434 reproduces EPS-01; truncation at 140 terms is below 3e-11 at these heights), the path lengths from the cusp are: 0.178 to the closest points
   over the 2-torsion points tau/2 and (1+tau)/2 (height 0.0335), 0.189 over the real 2-torsion point (height 0.0276), 0.107 to the Atkin–Lehner point i/sqrt(143),
   and 0.214 to the opposite cusp 0 (twice that, since W_143 is an isometry).  The region above height 0.0335 holds 15 per cent of the area.
2. The size of the curve (ESTIMATE).  The physical separation is sqrt(A_X) times these, plus the torus offset of the point's fibre.  V_cb = 0.041 with the (r/d)^3
   tail of TA-01 fixes the curve's area: A_X ~ 18 A_E (top blow-up A_E/3) to 34–36 A_E (0.55 A_E) for a top over a 2-torsion point on the graph of phi, and
   about 19–32 A_E for a top at another cusp.  With torus flux b = 0, D-flatness leaves A_X free, so the mixing measures it.
3. CC-37 revisited (STRUCTURAL).  In the canonical metric the other cusps are only slightly farther than the 2-torsion points: at A_X ~ 18 A_E the cusp-to-cusp
   mixing is 0.034–0.046, not <= 1/143.  CK-01's bound used the hyperbolic distance log 143 with exponential decay; the blown-up cusp is a smooth point of the
   physical metric, where tails fall as powers.  CK1-1 and CC-37 become caveats (CC-47): split assignments are not excluded by mixing alone.
4. A coincidence worth keeping (STRUCTURAL).  The points over the 2-torsion of E are where gamma-charged fermion zero modes localise (phi^-1(P_eta), NU-01): the
   top's candidate homes are the spin structure's points, or the Atkin–Lehner images of the jets' cusp.
OPEN: which home the dynamics picks; an independent handle on A_X (gauge couplings); the tail law beyond its leading scaling."""
import math
import numpy as np

PATHS = {"tau/2": (0.0335, 0.1781, 0.5629), "(1+tau)/2": (0.0335, 0.1781, 0.5629), "1/2": (0.0276, 0.1886, 0.4943)}
CUSP_TO_CUSP = 0.214
W143_POINT = 0.107


def density_setup(nmax=140, dps=30):
    from mtft.periods.channels import _setup
    a, Yi = _setup(nmax=nmax, dps=dps)
    A = np.array([[complex(a[n][i]) for i in range(13)] for n in range(len(a))])
    Y = np.array([[float(Yi[i, j]) for j in range(13)] for i in range(13)])
    return A, Y

def arakelov_density(tau, A, Y):
    q = np.exp(2j * np.pi * tau); w = (q ** np.arange(A.shape[0])) @ A
    return float(np.real(np.conj(w) @ Y @ w)) / 13

def path_length(x0, y0, A, Y, ytop=1.6, m=500):
    ys = np.concatenate([np.geomspace(y0, 0.2, m), np.linspace(0.2, ytop, m // 2)[1:]])
    return float(np.trapezoid([math.sqrt(max(arakelov_density(complex(x0, y), A, Y), 0.0)) for y in ys], ys))

def curve_area_for_Vcb(path, offset, eps2, V=0.041):
    r2 = math.sqrt(eps2 / math.pi); d = r2 / V ** (1 / 3)
    return (d ** 2 - offset ** 2) / path ** 2

def cusp_mixing(AX, eps2):
    return (math.sqrt(eps2 / math.pi) / (math.sqrt(AX) * CUSP_TO_CUSP)) ** 3

def report():
    out = {}
    for name, (h, L, off) in PATHS.items():
        for e in (1 / 3, 0.55):
            AX = curve_area_for_Vcb(L, off, e); out[(name, round(e, 3))] = (round(AX, 1), round(cusp_mixing(AX, e), 4))
    out[("cusp", 1 / 3)] = round(curve_area_for_Vcb(CUSP_TO_CUSP, 0.0, 1 / 3), 1); out[("cusp", 0.55)] = round(curve_area_for_Vcb(CUSP_TO_CUSP, 0.0, 0.55), 1)
    return out
