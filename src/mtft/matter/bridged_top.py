"""YB-01 (2026-10-08): the bridged top — the second redesign fails its pre-registered test.

Moving the left-handed top to the jet point (so that the right-handed down 2–3 angle can be of order one, DS-3) leaves its two halves at different points: the
top Yukawa becomes a bridge, the right-handed top's tail evaluated at the jet point.

1. The bridge (ESTIMATE).  With normalised constant jets of radii r1 (jet point) and r2 (top's home), the bridged-to-local ratio is 4 (r1/r2)^2 t(s), t the top's
   tail value at the jet point.  The image-summed cubic tail at the Fricke separation s = 1.62–1.79 (r = 230–280) gives t = 0.25–0.28 and a ratio of 0.41–0.45.
2. The test (pre-registered: y_t >= 0.5).  The local constant-jet top coupling is at most the marginal 0.70 (TP-01, LB-4), so y_t <= 0.29–0.32 with the flat tail
   and 0.17–0.27 with AG-01's curvature factors (0.60–0.85): FAIL, by a factor of 1.6 or more.
3. Verdict (STRUCTURAL).  Both redesigns PA-4 left open are closed: the B stack cannot carry an odd source (OM-2, topological) and a bridged top is too weak.
   The two-point design on X0(143) x 143a1 with this stack does not carry the quark flavour structure.  The arithmetic results — the canonical kernel (AG-01), the
   Fricke placement (FH-01), the concordance of r near 230–280 — stand on their own."""
import math

R1, R2 = math.sqrt(0.225 / math.pi), math.sqrt(0.55 / math.pi)
Y_LOCAL_MAX = 0.70
THRESHOLD = 0.5
CURVATURE = (0.60, 0.85)
SEPARATION = (1.62, 1.79)


def tail_value(s, offset=0j):
    from mtft.matter.bl_wilson_line import image_sum
    return image_sum(s, offset, (0, 0), p=3) * R2 ** 3

def jet_factor():
    return 4 * (R1 / R2) ** 2

def bridged_ratio(s):
    return jet_factor() * tail_value(s)

def y_bridged(s, y_local=Y_LOCAL_MAX, curvature=1.0):
    return y_local * bridged_ratio(s) * curvature

def outcome():
    flat = [y_bridged(s) for s in SEPARATION]
    curved = [y_bridged(s, curvature=c) for s, c in ((SEPARATION[1], CURVATURE[0]), (SEPARATION[0], CURVATURE[1]))]
    return {"flat": flat, "curved": curved, "passes": max(flat) >= THRESHOLD, "shortfall": THRESHOLD / max(flat)}
