"""YA-01 (2026-10-04): a character audit of the Yukawa amplitude for jet families — why the minimal background does not rescue the leading coupling, and
what sets the textures instead.

1. The overall rotation (EXACT, local model).  Near the blown-up point the overall U(1), (x, w) -> e^{i alpha}(x, w), is a symmetry in every background we use: the
   blow-up metric and exceptional flux are U(2)-invariant, curve-degree and torus fluxes are U(1) x U(1)-invariant, and the gamma-flux has character (1, -1).
   Characters as forms: each family jet +2 (degree-2 coboundary classes); a Phi-type Higgs at Taylor order m, m + 2 (its dx ^ dw legs); the condensate density
   from jets of degree <= 2, between -2 and 2.  An amplitude in which the Higgs enters holomorphically — tree level (YK-01) or PB-01's D-term mixing, whose
   coefficient is int D conj(phi_X) phi_H — has total character 4 + (m + 2) + d >= 4 > 0 and vanishes.  PB-01's route does not evade YK-01 by itself, and the
   minimal gamma-neutral background, whose Higgs is Phi-type, falls under the same rule.
2. The scales (EXACT given the backgrounds).  The rule holds only as far as the blow-up is small.  The exceptional area is A_E/3, a linear size 0.58 of the torus,
   so the torus rotation is broken by the lattice of 143a1 at the blow-up scale; the curve is larger (A_X > A_E) and its Bergman data nearly rotation-symmetric
   (XQ-3), so the curve rotation survives approximately.
3. Curve charge alone (EXACT under 2).  With only the curve rotation kept, (a + 1) + (l - k) + (4 - n) = 0 for Higgs curve order a >= 0 and density curve charge
   l - k >= -2, which allows only n = 3, 4: the leading texture is [[0,0,0],[0,0,Y3],[0,Y3,Y4]] — rank 2, one family massless — with Y3 from the Higgs's value at
   the cusp times the condensate's c2 c0* component and Y4 from its value and slope times the tilted components.  An untilted condensate gives nothing; a coherent
   tilt t gives m2/m3 ~ t^2.
4. Consequence (STRUCTURAL).  The jet families' Yukawas are not leading local amplitudes: they are set by how the global geometry — mostly the torus lattice at the
   blow-up scale — breaks the point's rotations, and they need the condensate's tilt.  PB-01's mechanism, EPS-01's Y ~ D_{q+q'} and HJ-01's leading texture tracked
   torus charge only and become caveats (CC-44).
OPEN: the lattice breaking at the blow-up scale (the torus lattice's Fourier components around the point) and whether it gives an O(1) top Yukawa."""
import math
import numpy as np

NU = np.array([1.0, math.sqrt(2.0), 1.0])
FAMILY_CHARACTER = 2
DENSITY_RANGE = (-2, 2)


def higgs_character(m):
    return m + 2

def holomorphic_amplitude_character(m, d):
    return 2 * FAMILY_CHARACTER + higgs_character(m) + d

def forbidden_holomorphic(mmax=12):
    return all(holomorphic_amplitude_character(m, d) > 0 for m in range(mmax + 1) for d in range(DENSITY_RANGE[0], DENSITY_RANGE[1] + 1))

def curve_allowed(amax=6):
    """(n, a, k - l) with (a + 1) + (l - k) + (4 - n) = 0, a >= 0, |l - k| <= 2."""
    return [(n, a, k_minus_l) for n in range(5) for a in range(amax + 1) for k_minus_l in range(-2, 3)
            if (a + 1) - k_minus_l + (4 - n) == 0]

def scale_ratio(exceptional_over_torus=1 / 3):
    return math.sqrt(exceptional_over_torus)

def coherent_amplitudes(t):
    c = np.array([1.0, math.sqrt(2.0) * t, t * t]); return c / np.linalg.norm(c)

def leading_texture(t, h0=1.0, h1=0.0):
    c = coherent_amplitudes(t)
    Y3 = h0 * c[2] * c[0]; Y4 = h0 * (c[1] * c[0] + c[2] * c[1]) + h1 * c[2] * c[0]
    Yn = [0.0, 0.0, 0.0, Y3, Y4]
    return np.array([[NU[q] * NU[r] * Yn[q + r] for r in range(3)] for q in range(3)])

def masses(Y):
    s = np.sort(np.linalg.svd(Y, compute_uv=False))[::-1]; return s / s[0]

def report():
    return {"holomorphic_forbidden": forbidden_holomorphic(), "curve_allowed": curve_allowed(), "scale_ratio": scale_ratio(),
            "texture_t0.1": masses(leading_texture(0.1)), "texture_untilted": float(np.abs(leading_texture(0.0)).max())}
