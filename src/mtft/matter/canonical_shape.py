"""CS-01 (2026-10-05): is there a canonical prescription for the curve-to-torus area ratio?  Pre-registered test.

Pre-registration (fixed before any candidate was evaluated).  Target: V_cb's window A_X/A_E in [18.1, 36.4] (SP-2).  A candidate counts as a hit only if
(i) its definition never refers to the target, (ii) it is the natural normalisation within its family, and (iii) its whole convention range lies inside the window.
Every candidate is logged, hits and misses.  Chance level: for a log-uniform prior on [1, 1000] the window covers ln(36.4/18.1)/ln 1000 = 0.10 of the range, so
with K candidates the chance of at least one accidental hit is 1 - 0.9^K.
1. Pullback through the modular parametrisation (EXACT): phi^*(flat torus) gives the curve area deg(phi) A_E, so A_X/A_E = 4.  Miss.  (The Petersson norm,
   4 pi^2 ||f||^2 = deg(phi) vol(E), is the same metric: also 4.)
2. Hyperbolic area over the Neron lattice (NUMERICAL): X_0(143) has hyperbolic area (pi/3)(168) = 56 pi (no elliptic points, 4 cusps); the Neron lattice of 143a1 has
   area Omega+^2 Im tau = 3.9709.  Ratio 44.30, or 22.15 if the torus volume is taken as int |omega ^ omega-bar| (factor 2).  The convention range straddles the
   window: not a hit by rule (iii).
3. Canonical metrics on both factors (EXACT): Arakelov area 1 against the flat torus of area 1 gives 1.  Miss.
4. Dimensionless invariants of the curve, the fractal and Hecke family (EXACT): genus 13, |chi| = 24, hyperbolic area / 2 pi = 28, index 168, cusp widths 1, 11, 13,
   143, four cusps, Hecke degrees p + 1.  Three of them (24, 26 = 2g, 28) fall in the window, but none comes with a derivation that makes it an area ratio: they fail
   rule (i)-(ii) as prescriptions and are recorded as numerology.
Verdict (STRUCTURAL).  No canonical prescription survives.  Four families tested; at chance level 0.10 each, about one accidental hit was expected, and the near-misses
(22.15, and the genus-type counts 24-28) are what a window this wide collects from a genus-13 curve.  The shape stays a parameter that V_cb measures.
OPEN: a physical mechanism (flux quantisation of the Kaehler class, a non-perturbative potential) that would make one of these numbers an area ratio."""
import math

WINDOW = (18.1, 36.4)
OMEGA_PLUS = 1.9699231645720703
IM_TAU = 1.0232745926964613
DEG_PHI = 4
INDEX = 168


def chance_per_candidate(lo=1.0, hi=1000.0):
    return math.log(WINDOW[1] / WINDOW[0]) / math.log(hi / lo)

def candidates():
    vol = OMEGA_PLUS ** 2 * IM_TAU
    hyp = math.pi / 3 * INDEX
    return {"pullback_phi": (DEG_PHI, DEG_PHI), "hyperbolic_over_neron": (hyp / (2 * vol), hyp / vol), "canonical_metrics": (1.0, 1.0)}

def invariants():
    g = 13
    return {"genus": g, "|chi|": 2 * g - 2, "2g": 2 * g, "hyp_area_over_2pi": 2 * g - 2 + 4, "index": INDEX, "cusp_widths": (1, 11, 13, 143), "cusps": 4}

def hit(rng):
    lo, hi = rng; return WINDOW[0] <= lo and hi <= WINDOW[1]

def report():
    c = candidates(); inv = invariants()
    in_window = [k for k, v in inv.items() if isinstance(v, int) and WINDOW[0] <= v <= WINDOW[1]]
    return {"candidates": c, "hits": {k: hit(v) for k, v in c.items()}, "invariants_in_window": in_window, "chance": chance_per_candidate(),
            "expected_accidental": 1 - (1 - chance_per_candidate()) ** 4}
