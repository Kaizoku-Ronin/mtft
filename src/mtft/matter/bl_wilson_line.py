"""BL-01 (2026-10-06): a B–L Wilson line on 143a1 — the right shape, the wrong size, and an unstable vacuum.

A flat Z/2 holonomy on the colour stack (quark bifundamentals twisted by a character chi of the torus lattice, lepton ones untwisted) leaves every local entry alone
and multiplies each tail by an image-sum ratio R = |sum chi(lam) t(lam)| / sum t(lam).  Lepton tails then exceed quark tails by 1/R; CI-3 needs 1/R = 1.758 per chirality.
Tails: the curve-gradient of the cubic tail (r_2/d)^3 that DD-01 uses, summed over torus images (weight s_X d^-5), with P2 over a 2-torsion point.

1. Two of the three characters kill the quark tails (EXACT).  Over a 2-torsion point the images pair up at equal distance; if chi is -1 on the doubled offset every pair
   cancels and the quark 2–3 entries vanish identically.  One character per 2-torsion point survives, with R falling from 0.92 (over 1/2) or 0.79 (over tau/2) as the
   curve separation s_X grows.
2. The size is wrong (NUMERICAL).  R = 1/1.758 at s_X = 0.716 (over 1/2) or 0.499 (over tau/2); there the twisted quark V_cb is 0.133 or 0.184, against the pre-registered
   window [0.030, 0.055]: a factor 2.4–6.1 too large.  Moving out until V_cb = 0.055 to 0.030 makes the lepton tails 2.6 to 4.2 times the quark tails, so the muon factor
   comes out 6.6–18 instead of 3.09.  The Wilson line gives the lepton pattern or V_cb, never both.
3. The vacuum is a saddle (EXACT at one loop, matched spectra).  The flux-free blocks charged under the holonomy see V = -sum over odd Scherk–Schwarz windings of
   cos(theta.w) |w|^-6.  At every Z/2 point with chi distinct from the Scherk–Schwarz character V lies above theta = 0 and its Hessian has a negative eigenvalue (all six
   pairings; also for |w|^-4 and |w|^-8).  With chi equal to it, the leptoquark gauginos become periodic and massless.
4. By-product (NUMERICAL).  Over a 2-torsion point the untwisted tail is an image sum led by the equidistant pair: at s_X = 0.85 it is 3.1 (over 1/2) to 3.5
   (over tau/2) times DD-01's single-image estimate, which moves the separation V_cb wants outward by about 1.3.
OPEN: whether any flavour mechanism remains for the lepton pattern; the matter programme's flavour failures audited together."""
import math
import numpy as np

R1, R2 = math.sqrt(0.225 / math.pi), math.sqrt(0.55 / math.pi)
WINDOW = (0.030, 0.055)
_CACHE = {}


def tau_143a1():
    if "tau" not in _CACHE:
        import mpmath as mp
        t = mp.findroot(lambda t: mp.re(mp.kleinj(mp.mpc(0.5, t)) * 1728) + mp.mpf(262144) / 1859, 1.02)
        _CACHE["tau"] = complex(0.5, float(t))
    return _CACHE["tau"]

def _grid(N=50):
    key = ("grid", N)
    if key not in _CACHE:
        tau = tau_143a1(); s = 1 / math.sqrt(tau.imag)
        m, n = np.meshgrid(np.arange(-N, N + 1), np.arange(-N, N + 1), indexing="ij")
        _CACHE[key] = (m, n, (m + n * tau) * s, s)
    return _CACHE[key]

def half_period(label):
    tau = tau_143a1(); _, _, _, s = _grid()
    return {"1/2": 0.5, "tau/2": tau / 2, "(1+tau)/2": (1 + tau) / 2}[label] * s

def character(c):
    m, n, _, _ = _grid(); return (-1.0) ** (c[0] * m + c[1] * n)

def image_sum(sx, sE, c=(0, 0), p=5):
    _, _, lam, _ = _grid()
    return float(np.sum(character(c) * (sx ** 2 + np.abs(sE + lam) ** 2) ** (-p / 2)))

def ratio(sx, label, c):
    sE = half_period(label); return abs(image_sum(sx, sE, c)) / image_sum(sx, sE)

def nearest(label):
    _, _, lam, _ = _grid(); return float(np.min(np.abs(half_period(label) + lam)))

def vcb(sx, label, c=(0, 0)):
    return 3 * R1 * R2 ** 3 * sx * abs(image_sum(sx, half_period(label), c))

def _solve(f, target, lo=0.02, hi=3.0):
    for _ in range(70):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) > target else (lo, mid)
    return (lo + hi) / 2

def lepton_point(label, c, r_target=1 / math.sqrt(3.092)):
    sx = _solve(lambda x: ratio(x, label, c), r_target)
    return {"s_X": sx, "d": math.hypot(sx, nearest(label)), "vcb_twisted": vcb(sx, label, c), "vcb_untwisted": vcb(sx, label)}

def vcb_point(label, c, v):
    sx = _solve(lambda x: vcb(x, label, c), v)
    r = ratio(sx, label, c); return {"s_X": sx, "lepton_over_quark_tail": 1 / r, "muon_factor": 1 / r ** 2}

def potential(theta, ss, p=6):
    m, n, w, _ = _grid(); L = np.abs(w); L[m.shape[0] // 2, m.shape[1] // 2] = np.inf
    odd = character(ss) == -1.0
    return float(-np.sum(np.cos(theta[0] * m + theta[1] * n)[odd] * (L[odd] ** -p)))

def hessian_at(ss, c, p=6):
    m, n, w, _ = _grid(); L = np.abs(w); L[m.shape[0] // 2, m.shape[1] // 2] = np.inf
    odd = character(ss) == -1.0; f = (L ** -p) * character(c)
    H = np.array([[np.sum((m * m * f)[odd]), np.sum((m * n * f)[odd])], [np.sum((m * n * f)[odd]), np.sum((n * n * f)[odd])]])
    return np.linalg.eigvalsh(H)

def image_factor(sx, label):
    """Untwisted image sum over the single nearest image, at curve separation sx."""
    return image_sum(sx, half_period(label)) / (sx ** 2 + nearest(label) ** 2) ** -2.5

def report():
    pts = {k: lepton_point(*k) for k in (("1/2", (0, 1)), ("tau/2", (1, 0)))}
    edges = {k: [vcb_point(*k, v) for v in (WINDOW[1], 0.041, WINDOW[0])] for k in pts}
    saddles = {(ss, c): hessian_at(ss, c).tolist() for ss in ((1, 0), (0, 1), (1, 1)) for c in ((1, 0), (0, 1), (1, 1)) if ss != c}
    return {"lepton_points": pts, "vcb_points": edges, "saddles": saddles}
