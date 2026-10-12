"""AG-01 (2026-10-07): the canonical Green's-function kernel of X0(143) — the pre-registered test of FA-01, run.

1. The kernel is conformal (EXACT).  Green's functions for different volume forms differ by f(P) + f(Q) + C, so every mixed second derivative is the same for
   all of them: d_P dbar_Q g = pi w(P)^T Y^-1 conj(w(Q)) and d_P d_Q g = B(P, Q) - pi w(P)^T Y^-1 w(Q), with B the fundamental bidifferential.  The tail kernel is a
   conformal invariant of X0(143); the canonical (Arakelov) metric enters only through the local normalisation sqrt(mu(P) mu(Q)).
2. The second-kind differential (EXACT construction, NUMERICAL checks).  No weight-2 eta quotient has a lone double pole at the cusp, and the smallest
   modular function with poles only there has order 60.  It is eta(t)eta(11t)/(eta(13t)eta(143t)) — orders (-6, 6, 6, -6) at (inf, 0, 1/11, 1/13) — times a
   weight-2 cusp form vanishing to order 5 at infinity and 7 at 1/13.  Its residue vanishes to 3e-14; its holomorphic part is fixed by Petersson orthogonality
   over the Ford domain of the Atkin–Lehner quotient (area 14 pi to 3e-6), whose integration reproduces Im Omega to 8e-6.  Eta factors are evaluated as
   products: the q-series of a weakly holomorphic form loses ~2 pi / y digits at height y and is unusable there.
3. The kernel (NUMERICAL).  K = |d_P d_Q g| / sqrt(mu mu) in area-one units from the blown-up cusp: 22.87 to the tau/2 preimage, 74.55 to the Fricke point,
   13.02 to the cusp 0, 1.73 and 1.27 to the cusps 1/11 and 1/13; a second construction (shifted by a holomorphic differential) agrees to <= 0.8%.  Against the
   flat 1/l^2 at the same canonical distance: 0.73, 0.85, 0.60 — the curvature correction FA-3 flagged is 15–40%.
4. The pre-registered test (NUMERICAL).  r = 2 pi r1 r2^3 K / V_cb (2D regime), V_cb in [0.030, 0.055], against the Cabibbo window 91–896 (CA-1):
   Fricke point 167–306, central 224 — PASS; tau/2 preimage 51–94, central 69 — touches the edge only; cusp 0: 29–53; cusps 1/11, 1/13: 3–7 — FAIL.
   Of the four arithmetic homes the CM point passes and the three rational ones fail.
5. Reading (STRUCTURAL).  One pass in four homes, at a point CA-3 had already flagged post hoc: modest evidence by the look-elsewhere count.  It fixes the top's
   home at the Fricke point and r near 170–300, consistent with V_us (CA-1) and NU-02's side.
Provenance: scripts/arakelov_kernel_x0143.py (needs cypari2); the numbers below are frozen from it."""
import math

R1, R2 = math.sqrt(0.225 / math.pi), math.sqrt(0.55 / math.pi)
VCB = (0.030, 0.041, 0.055)
CABIBBO_R = (90.8, 896.0)
MULTIPLIER = (1, 1, -1, -1)                       # eta exponents at (1, 11, 13, 143)
KERNEL = {"Fricke point": (74.5499, 74.5440), "tau/2 preimage": (22.8695, 22.8746), "cusp 0": (13.0238, 13.0239),
          "cusp 1/11": (1.7334, 1.7229), "cusp 1/13": (1.2714, 1.2819)}
FLAT_L = {"Fricke point": 0.107, "tau/2 preimage": 0.1781, "cusp 0": 0.214}
ARITHMETIC = ("Fricke point", "cusp 0", "cusp 1/11", "cusp 1/13")
VALIDATION = {"ford_area": 43.98218, "target_area": 14 * math.pi, "gram_rel_err": 8.2e-6, "mu_cusp": 0.0243397, "residue": 3.0e-14}


def multiplier_orders(r=MULTIPLIER):
    r1, r11, r13, r143 = r
    return ((r1 + 11*r11 + 13*r13 + 143*r143) / 24, (143*r1 + 13*r11 + 11*r13 + r143) / 24,
            (13*r1 + 143*r11 + r13 + 11*r143) / 24, (11*r1 + r11 + 143*r13 + 13*r143) / 24)

def kernel(home):
    a, b = KERNEL[home]; return (a + b) / 2

def spread(home):
    a, b = KERNEL[home]; return abs(a - b) / ((a + b) / 2)

def ratio_to_flat(home):
    return kernel(home) * FLAT_L[home] ** 2

def r_from_vcb(home, v):
    return 2 * math.pi * R1 * R2 ** 3 * kernel(home) / v

def r_range(home):
    return (r_from_vcb(home, VCB[2]), r_from_vcb(home, VCB[0]))

def passes(home):
    """Pre-registered criterion: r at the central V_cb lies in the Cabibbo window."""
    return CABIBBO_R[0] <= r_from_vcb(home, VCB[1]) <= CABIBBO_R[1]

def touches(home):
    lo, hi = r_range(home); return max(lo, CABIBBO_R[0]) <= min(hi, CABIBBO_R[1])

def report():
    return {h: {"K": kernel(h), "r_range": r_range(h), "r_central": r_from_vcb(h, VCB[1]), "passes": passes(h), "touches": touches(h)} for h in KERNEL}
