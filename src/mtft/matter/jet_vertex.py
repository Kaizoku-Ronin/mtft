"""HJ-01 (2026-10-04): the family vertex of the jet classes — selection rule, locality, the round blow-up, and where X_0(143) enters.

1. Selection rule (EXACT in the local model).  The blown-up point is rotation-symmetric in the curve direction near the point.  The jet classes come from
   holomorphic jets through an equivariant coboundary, so family q carries curve charge 2 - q (as a form).  A messenger whose curve factor is a holomorphic
   differential (the curve copies of PA-01/PB-01, whose sum is the Bergman kernel) carries curve charge beta + 1 >= 1 at Taylor order beta.  The vertex
   phi ^ psi ^ psi needs total charge 0, but (4 - n) + (beta + 1) > 0 always: these messengers do not couple at leading order, only through the breaking of
   the point's rotation symmetry.  Messengers that do couple must depend antiholomorphically on the curve, i.e. be massive along it.
2. Locality (STRUCTURAL).  Both vertices sit at the blown-up point, so the coupling messengers' propagator is dominated by the blow-up scale, where the rotation
   symmetry holds: charge-conserving and universal at leading order.  The global data of X_0(143) then enter through the light Higgs's own profile — its
   Taylor coefficients at the cusp, which are q-expansion coefficients when the profile is a cusp form — with the regular part of the curve Green's function
   as a correction.
3. The round blow-up gives no hierarchy (EXACT).  With the full SU(2) the vertex depends on the families only through their product (spin 2), and a coherent
   condensate supplies only the m = 0 component along its axis, so Y = R2 [[0,0,1],[0,2,0],[1,0,0]] in the axis frame: singular values 2 : 1 : 1 for every tilt
   (checked on random axes).  The hierarchy must come from the curve/torus anisotropy and the Higgs's Taylor data.
4. Leading texture (MODEL).  For a coherent condensate along x^2: Y_qq' ~ nu_q nu_q' h_(4-n) t_n — the Higgs profile's curve Taylor coefficient of order 4 - n
   (up to a fixed offset from the form legs) times its torus Taylor coefficient of order n, times local weights.  Illustrative curve log-convexities
   h4 h2 / h3^2: 0.25 for the 143a1 newform (a5 a3 / a4^2), 0.43 for the profile peaked at the cusp (Bergman first column).  Which light-Higgs combination
   takes the vev is set by the Higgs potential.
OPEN: the light Higgs's curve and torus profiles in the ratio-3 background; the torus Taylor data; the order offset; the condensate tilt's contribution."""
import math
import numpy as np

NU = np.array([1.0, math.sqrt(2.0), 1.0])


def curve_charge(n, beta):
    """Total curve charge of the vertex for a family pair of torus charge n and a holomorphic messenger of Taylor order beta."""
    return (4 - n) + (beta + 1)

def _mul(p, q):
    out = np.zeros(len(p) + len(q) - 1, complex)
    for i, x in enumerate(p):
        for j, y in enumerate(q): out[i + j] += x * y
    return out

def round_texture(alpha, beta):
    """Texture of a fully SU(2)-symmetric blow-up for a coherent condensate along the axis (alpha, beta)."""
    up = np.array([alpha, beta]); dn = np.array([-np.conj(beta), np.conj(alpha)])
    P = _mul(_mul(up, up), _mul(dn, dn))
    R = np.array([np.conj(P[n]) * math.factorial(4 - n) * math.factorial(n) / 24 for n in range(5)])
    return np.array([[NU[q] * NU[r] * R[q + r] for r in range(3)] for q in range(3)])

def singular_ratios(Y):
    s = np.sort(np.linalg.svd(Y, compute_uv=False))[::-1]; return s / s[-1]

def newform_coefficients(nmax=12):
    from mtft.periods.forms import _elliptic_an, CURVE_143A1
    return _elliptic_an(CURVE_143A1, 143, nmax)

def log_convexity(h):
    """h4 h2 / h3^2 for Taylor coefficients h_0..h_4 (h_a = q-coefficient a + 1)."""
    return h[4] * h[2] / h[3] ** 2

def profiles():
    a = newform_coefficients(); new = [a[k + 1] for k in range(5)]
    from mtft.matter.overlaps import cusp_bergman_jets
    b = cusp_bergman_jets()[:, 0]; berg = [b[k] for k in range(5)]
    return {"newform": (new, log_convexity(new)), "bergman_peaked": (berg, log_convexity(berg))}

def taylor_texture(h, t):
    """nu-Hankel texture with entries Y_n = h_(4-n) t_n, n = 0..4."""
    Yn = [h[4 - n] * t[n] for n in range(5)]
    return np.array([[NU[q] * NU[r] * Yn[q + r] for r in range(3)] for q in range(3)])

def report():
    rng = np.random.default_rng(0); axes = [rng.normal(size=2) + 1j * rng.normal(size=2) for _ in range(5)]
    ratios = [singular_ratios(round_texture(*(v / np.linalg.norm(v)))) for v in axes]
    return {"min_curve_charge": min(curve_charge(n, b) for n in range(5) for b in range(12)),
            "round_ratios": ratios, "profiles": profiles()}
