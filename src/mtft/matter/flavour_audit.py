"""FA-01 (2026-10-06): the flavour audit — one shared artefact found, four independent failures confirmed, and the placement fixed by arithmetic.

1. Scoreboard (STRUCTURAL).  Eight flavour observables, six continuous parameters introduced to meet them (eps_P, c, the P1–P2 separation, the dressing ratio z, the
   condensate ratio c1/c2, the tail ratios), two pre-registered tests (CS-01 no hit, TO-01 failed).  Nothing has yet passed a test it was not tuned to.
2. The tail model was single-image (NUMERICAL).  Every separation and every r taken from V_cb used DD-01's tail from one torus image.  On the flat torus the tail is an image
   sum; beyond the torus size it crosses over to 2D, V_cb -> 2 pi r1 r2^3 / s_X^2 instead of 3 r1 r2^3 / d^4.  Redone with image sums (V_cb window 0.055–0.030), V_cb's r rises by
   2.5–3.4: from 22–34 to 63–130 over the 2-torsion points and from 90–122 to 196–359 at the Fricke point, both overlapping the Cabibbo window 91–896 (CA-1); the opposite
   cusp (49–90) touches its edge.  The factor-3.3 gap of CA-2 was mostly the model, not the geometry.  Post hoc: nothing was tuned, but the correction was found after the failure.
3. The curve is not flat on that scale (EXACT, Gauss–Bonnet).  Genus 13 forces total curvature -48 pi, so the mean curvature radius is sqrt(A_X/48 pi) and a canonical
   path of length l (area-1 units) spans l sqrt(48 pi) = 12.28 l radii, at every A_X: 1.3 to the Fricke point, 2.2–2.3 to the 2-torsion preimages, 2.6 to the
   opposite cusp.  The flat-curve tail is uncontrolled; the right object is the Green's function of the canonical (Arakelov) metric.
4. Placement by arithmetic (EXACT, from the literature).  X0(143) has genus 13, so it has finitely many rational points (Faltings); by Mazur and Kenku no elliptic curve over
   Q has a rational cyclic 143-isogeny, so X0(143)(Q) is exactly its four cusps (all rational, N squarefree).  P1 sits on one; P2's curve coordinate has been a fitted
   continuous parameter.  Freezing P2 to arithmetic points (cusps; CM points such as the Fricke point) removes that parameter.
5. What survives the correction (STRUCTURAL).  TO-01's failure is a ratio of tails (unchanged by any common normalisation), the lepton failure is exact (CI-2, BL-1),
   b–tau normalisation and the Pati–Salam scale do not involve tails.  Four independent open failures remain.
OPEN: the Arakelov Green's function of X0(143) between the blown-up cusp and the candidate homes (closed forms exist for cusps of X0(N), N squarefree and prime to 6);
the pre-registered version of item 2."""
import math

R1, R2 = math.sqrt(0.225 / math.pi), math.sqrt(0.55 / math.pi)
CABIBBO_R = (90.8, 896.0)                       # CA-1, gamma-charge 1
VCB = (0.030, 0.041, 0.055)
CANONICAL = {"over 1/2": (0.189, "1/2"), "over tau/2": (0.178, "tau/2"), "Fricke": (0.107, None), "opposite cusp": (0.214, None)}
MAZUR_KENKU = set(range(1, 20)) | {21, 25, 27, 37, 43, 67, 163}

SCOREBOARD = [
    ("V_us = 0.225", "fit", "eps_P from the gamma-tilt; needs r = 91-896", "CA-1, DD-01"),
    ("V_cb = 0.041", "fit", "sets the P1-P2 separation and so r", "DD-4, SP-2, CA-2"),
    ("V_ub", "fit", "coefficient c = 0.40, consistent but not sharp", "CA-4"),
    ("m_s/m_b, m_c/m_t with V_cb", "failed (pre-registered)", "symmetric tails miss by 1.9; asymmetric tails or dressing needed", "TO-01, FL-01, CD-01"),
    ("lepton pattern 3.09, 0.30", "failed", "dressing, insertion and B-L Wilson line all closed", "QL-01, CI-01, BL-01"),
    ("y_tau/y_b = 1.40-1.53", "open", "needs a 30-35% threshold", "GC-01"),
    ("Pati-Salam scale R = 3.68", "open", "flux designs quantise R; thresholds need 17+ pairs", "GC-01, TH-01"),
    ("two measurements of r", "post-hoc pass", "single-image artefact; image sums overlap the Cabibbo window", "CA-2, FA-01"),
]


def _bl():
    from mtft.matter import bl_wilson_line as B
    return B

def _sE(label):
    return 0j if label is None else _bl().half_period(label)

def vcb_images(sx, label):
    B = _bl(); return 3 * R1 * R2 ** 3 * sx * B.image_sum(sx, _sE(label))

def vcb_2d(sx):
    return 2 * math.pi * R1 * R2 ** 3 / sx ** 2

def sx_for(v, label):
    lo, hi = 0.05, 6.0
    for _ in range(70):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if vcb_images(mid, label) > v else (lo, mid)
    return (lo + hi) / 2

def r_images(name, v):
    l, label = CANONICAL[name]; return (sx_for(v, label) / l) ** 2

def r_single(name, v):
    """DD-01's single-image reading: 3 r1 r2^3 / d^4 = v with d^2 = s_X^2 + |s_E|^2."""
    l, label = CANONICAL[name]
    d = (3 * R1 * R2 ** 3 / v) ** 0.25
    sE = 0.0 if label is None else _bl().nearest(label)
    return (math.sqrt(max(d * d - sE * sE, 0.0)) / l) ** 2

def curvature_radii(l):
    return l * math.sqrt(48 * math.pi)

def cusps_rational(N=143):
    import sympy
    ds = sympy.divisors(N); return len(ds), all(math.gcd(d, N // d) == 1 for d in ds)

def no_rational_isogeny(N=143):
    return N not in MAZUR_KENKU

def overlap(rng, window=CABIBBO_R):
    return max(rng[0], window[0]) <= min(rng[1], window[1])

def report():
    out = {}
    for name in CANONICAL:
        new = [r_images(name, v) for v in VCB]; old = [r_single(name, v) for v in VCB]
        out[name] = {"r_single (0.030, 0.041, 0.055)": old, "r_images": new, "overlaps_cabibbo": overlap((min(new), max(new))),
                     "curvature_radii": curvature_radii(CANONICAL[name][0])}
    return out
