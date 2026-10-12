"""CP-01 (2026-10-03): the condensate's jet alignment — spin-1 phases, the tree-level potential, and a Georgi–Jarlskog structure in the texture.

1. Spin-1 structure (EXACT).  The local blow-up of a point of C^2 is U(2)-symmetric; its degree-2 jets x^2, x w, w^2 carry the spin-1 representation of SU(2)
   (orthonormal basis x^2, sqrt(2) x w, w^2 for the U(2)-invariant inner product).  An SU(2)-invariant quartic on spin-1 has exactly two invariants,
   (c^dag c)^2 and |theta(c)|^2 with theta(c) = c_0^2 - 2 c_+ c_-, so a quartic potential has two phases: COHERENT (theta = 0: a perfect square (a x + b w)^2)
   and POLAR (|theta| = |c|^2: a product of orthogonal lines (a x + b w)(-b* x + a* w)).  Anisotropy between the curve and torus directions orients them.
2. The tree-level potential prefers polar (EXACT for the local quartic).  The D-term quartic is proportional to the integral of |Sigma|^4 over the region where the
   jets live; on the sphere S^3 of directions E|s|^4 = (1/5)(c^dag c)^2 - (1/15)|theta|^2 (checked to 3e-14 on random states); normalised by (E|s|^2)^2 it runs
   from 1.8 (coherent, the maximum) to 1.2 (polar, the minimum).  At tree level the alignment is polar.
3. What each phase does to the texture (EXACT within the Hankel class Y = [[1, D1, D2], [D1, D2, 0], [D2, 0, 0]]).  m2/m3 = |D2 - D1^2|, m1/m3 = |D2|^3 / m2.
   A pure polar alignment along x w gives D1 = D2 = 0 (rank one); polar along x^2 - w^2 gives no hierarchy; a slightly tilted coherent alignment (orthonormal
   basis) gives D2 = D1^2/2 — the quark-like hierarchy.  The quark fits sit near coherent (rho = D2/D1^2 = 0.59 for down quarks), so the observed vacuum needs the
   spin-dependent quartic to change sign relative to tree level: a definite target for the one-loop (Legendre-transform) effective potential.
4. A Georgi–Jarlskog structure (EXACT within the Hankel class; MECHANISM OPEN).  If leptons see the charge-2 coupling with the opposite sign (D2 -> -D2, same D1):
   m3 is unchanged (b–tau), m2 is multiplied and m1 divided by the same factor F = (1 + rho)/(1 - rho), so m_mu m_e = m_s m_d exactly; rho = 1/2 (orthonormal-coherent)
   gives 3 and 1/3 in the limit of small D1 (2.91 and 0.345 at D1 = 0.2).  Data at the b–tau scale: m_mu/m_s = 3.06, m_e/m_d = 0.295, product 0.90.  With rho = 0.59 from the down-quark fit the full texture
   predicts m_mu/m_s = 3.57 (17 per cent high) and m_e/m_d = 0.282 (4 per cent low).  The sign flip was spotted after the LEP-01 fits, so it is a hypothesis to derive, not evidence; a reflection
   1 - 2 P_lepton through the condensate's direction on lepton lines would supply it.
5. [RETRACTED, CC-40] Where each sector sits (MODEL).  Reading each LEP-01 equal-overlap fit as a condensate profile, |theta|/|c|^2 (0 coherent, 1 polar) is 0.006 for up quarks,
   0.033 for down quarks and 0.926 for charged leptons: quarks coherent, leptons polar.  A coherent profile has D2 = +D1^2/2 and a polar-leaning one D2 < 0, so the
   sign flip of item 4 is the same distinction.  The condensate points in the lepton direction of the 4, so leptons are the natural sector to inherit the
   tree-level (polar) alignment; what gives quarks a coherent profile is open.
6. Zero flips (EXACT; source of CC-40).  The texture depends on the profile only through its autocorrelation, i.e. through |p|^2 on the unit circle; reflecting a
   zero of the profile z -> 1/conj(z) leaves it unchanged.  The lepton fit (1, 3.26, -0.31) and the x^2-dominated profile obtained by flipping its outer zero give the
   same masses (the flip gives (1, 0.205, -0.028)), so item 5's coherent/polar reading is an artefact.  The invariants are D1 and rho = D2 D0 / D1^2, and two mass ratios fix rho only up to two branches
   per sector: up 0.453 or -4.80, down 0.577 or 3.73, leptons 0.302 or -0.765 (from (|rho| / |1 - rho|)^3 = m1 m3 / m2^2).  Item 4's sign flip is one branch choice.
7. Parity and the charge-1 component (EXACT in the model).  w -> -w flips the x w jet only.  In the Landau problem — twist anisotropy kappa(|c_0|^2 + 2|c_-|^2)
   against the tree-level quartic -|theta|^2 at fixed norm — the ground state is x^2 for kappa >= 2 and an x^2 / w^2 mixture below, always with c_0 = 0.  Then D1 = 0
   and the light masses are geometric, m1 m3 = m2^2.  The data give m1 m3 / m2^2 = 0.57 (up), 2.55 (down), 0.081 (leptons): the observed families need the charge-1
   component, i.e. w -> -w broken in the condensate's alignment.
8. The sharpened one-loop target (STRUCTURAL).  In the ratio-3 background the condensate entry carries no gamma-flux, so its potential is parity-even at tree level;
   parity-odd terms enter through loops of gamma-charged fields (gamma -> -gamma under w -> -w).  That parity-odd term sets D1 — the quantity to compute next.
CP-03 (2026-10-03): the parity-odd term as a transverse field.
9. Zeeman structure (EXACT).  A constant (1,1) field strength at the point is a Hermitian 2x2 matrix: trace (Kaehler, U(2)-invariant) plus a spin-1 vector.  By
   Wigner–Eckart the degree-2 jets couple to it only through <c| B.J |c> (one reduced coefficient).  Curve/torus flux differences give the longitudinal component
   (P-even); the mixed gamma-flux dx ^ dw-bar gives the transverse component (P-odd, since gamma -> -gamma under w -> -w).
10. Coherent vacuum (EXACT).  The minimum of -B.J is the spin-coherent state along B: a tilted perfect square, so rho = 1/2 exactly for every tilt; D1 = 0.202
    at tilt 0.29 rad (B_perp/B_z = 0.30).
11. Stability against the tree-level quartic (EXACT).  With E = -<J.B_hat> - lambda |theta|^2, the coherent state is the exact ground state for lambda <= 1/2;
    at 1/2 it destabilises towards the opposite pole and rho runs far from 1/2 at once (-5.3 at lambda = 0.55).  Small departures from 1/2 need fine tuning.
12. Consequence (EXACT within the Hankel class).  rho = 1/2 gives m1 m3 = m2^2, as D1 = 0 does: parity breaking is necessary (CP-8) but not sufficient.  The data's
    departures — rho = 0.453 (up), 0.577 (down) on the first branch, i.e. m1 m3/m2^2 = 0.57 and 2.55 — must come from unequal overlaps X_qq', which the
    equal-overlap texture sets to 1.  Leptons (0.081) need more (the Georgi–Jarlskog question).
13. Size (STRUCTURAL).  For a gamma-charged condensate the transverse field is tree-level (EPS-01's estimate is this tilt); for the gamma-neutral ratio-3
    condensate it is one-loop (CP-9), suppressed by a loop factor times a charge trace, so a tilt of 0.29 rad there would need a large trace.
OPEN: the one-loop sign of the spin-dependent quartic; the mechanism for the lepton sign; the order-one coefficients beyond the Hankel idealisation."""
import itertools
import math
import numpy as np


# ------------------------------------------------------------------ 1-2. spin-1 states as degree-2 polynomials on S^3
def sphere_moment(a, b):
    """E |x|^{2a} |w|^{2b} for (x, w) uniform on S^3 in C^2."""
    return math.factorial(a) * math.factorial(b) / math.factorial(a + b + 1)

def poly(c):
    """c = (c_plus, c_0, c_minus) in the orthonormal basis x^2, sqrt(2) x w, w^2 -> {(i, j): coefficient of x^i w^j}."""
    return {(2, 0): c[0], (1, 1): math.sqrt(2) * c[1], (0, 2): c[2]}

def _mul(p, q):
    out = {}
    for (i1, j1), a in p.items():
        for (i2, j2), b in q.items():
            out[(i1 + i2, j1 + j2)] = out.get((i1 + i2, j1 + j2), 0) + a * b
    return out

def norm2(c):
    return sum(abs(v) ** 2 * sphere_moment(i, j) for (i, j), v in poly(c).items())

def quartic(c):
    """E |s|^4 = E |s^2|^2 (cross terms average out over the torus of phases)."""
    sq = _mul(poly(c), poly(c))
    return sum(abs(v) ** 2 * sphere_moment(i, j) for (i, j), v in sq.items())

def theta(c):
    return c[1] ** 2 - 2 * c[0] * c[2]

def fit_invariants(n=40, seed=0):
    rng = np.random.default_rng(seed); A, y = [], []
    for _ in range(n):
        c = rng.normal(size=3) + 1j * rng.normal(size=3)
        A.append([np.vdot(c, c).real ** 2, abs(theta(c)) ** 2]); y.append(quartic(c))
    coef, res, *_ = np.linalg.lstsq(np.array(A), np.array(y), rcond=None)
    pred = np.array(A) @ coef
    return {"a": coef[0], "b": coef[1], "max_error": float(np.max(np.abs(pred - y)))}

def coherent(t):
    """(x + t w)^2 normalised, in the orthonormal basis."""
    c = np.array([1.0, math.sqrt(2) * t, t * t]); return c / math.sqrt(np.vdot(c, c).real)

def polar(t):
    """(x + t w)(-t x + w) normalised (t real)."""
    c = np.array([-t, (1 - t * t) / math.sqrt(2), t]); return c / math.sqrt(np.vdot(c, c).real)


# ------------------------------------------------------------------ 3-4. the Hankel texture and the sign flip
def hankel_masses(D1, D2):
    Y = np.array([[1.0, D1, D2], [D1, D2, 0.0], [D2, 0.0, 0.0]])
    m = np.sort(np.linalg.svd(Y, compute_uv=False))[::-1]
    return m[1] / m[0], m[2] / m[0]

def D_from_alignment(c):
    """Autocorrelations of the orthonormal-basis amplitudes, normalised by D0."""
    c = np.asarray(c, float); D = [sum(c[q + n] * c[q] for q in range(3 - n)) for n in range(3)]
    return D[1] / D[0], D[2] / D[0]

def sign_flip_factor(rho):
    return (1 + rho) / (1 - rho)

def down_rho():
    """rho = D2 / D1^2 of the LEP-01 equal-overlap down-quark fit."""
    from mtft.matter.leptons import fit_alignment, DOWN_MZ
    s, _ = fit_alignment(DOWN_MZ); D = [sum(s[q + n] * s[q] for q in range(3 - n)) for n in range(3)]
    d1, d2 = D[1] / D[0], D[2] / D[0]
    return d1, d2, d2 / d1 ** 2

def gj_test():
    from mtft.matter.leptons import gj_at_matching
    g = gj_at_matching(); d1, d2, rho = down_rho()
    m2q, m1q = hankel_masses(d1, d2); m2l, m1l = hankel_masses(d1, -d2)
    return {"rho_down": rho, "F": sign_flip_factor(rho), "pred_mmu_over_ms": m2l / m2q, "pred_me_over_md": m1l / m1q,
            "data_mmu_over_ms": g["mmu_over_ms"], "data_me_over_md": g["me_over_md"], "data_product": g["mmu_over_ms"] * g["me_over_md"]}

def sector_phases():
    """|theta| / |c|^2 for the LEP-01 equal-overlap fits (raw monomial coefficients -> orthonormal spin-1 amplitudes)."""
    from mtft.matter.leptons import fit_alignment, UP_MZ, DOWN_MZ, LEPTON_POLE
    out = {}
    for name, m in (("up", UP_MZ), ("down", DOWN_MZ), ("leptons", LEPTON_POLE)):
        s, _ = fit_alignment(m); c = (s[0], s[1] / math.sqrt(2), s[2])
        out[name] = abs(theta(c)) / sum(x * x for x in c)
    return out

def zero_flip(s):
    """Reflect the outer zero of the profile s0 + s1 z + s2 z^2 through the unit circle; returns the new coefficients (same autocorrelation up to scale)."""
    z = np.roots([s[2], s[1], s[0]])
    k = int(np.argmax(np.abs(z))); z[k] = 1 / np.conj(z[k])
    c = np.real_if_close(np.poly(z)) * 1.0                   # monic in z^2: (z - z1)(z - z2) -> coefficients [1, -(z1+z2), z1 z2]
    new = np.array([c[2], c[1], c[0]], dtype=complex)
    if abs(new[2]) > abs(new[0]): new = new[::-1]                     # reversal s0 <-> s2 also preserves the autocorrelation
    return np.real_if_close(new / new[0])

def raw_masses(s):
    from mtft.matter.leptons import mass_ratios
    return mass_ratios([float(np.real(x)) for x in s])[:2]

def rho_branches(m1, m2, m3):
    """The two values of rho allowed by two mass ratios: |rho / (1 - rho)|^3 = m1 m3 / m2^2."""
    k = (m1 * m3 / m2 ** 2) ** (1 / 3)
    first = k / (1 + k)                                                  # 0 < rho < 1
    second = -k / (1 - k) if k < 1 else k / (k - 1)                      # rho < 0 (k < 1) or rho > 1 (k > 1)
    return first, second

def landau_ground(kappa, starts=40, seed=0):
    """Minimise kappa(|c0|^2 + 2|c-|^2) - |theta|^2 on the unit sphere (real amplitudes)."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed); best = None
    f = lambda v: (lambda c: kappa * (c[1] ** 2 + 2 * c[2] ** 2) - theta(c) ** 2)(v / np.linalg.norm(v))
    for _ in range(starts):
        r = minimize(f, rng.normal(size=3), method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 20000})
        if best is None or r.fun < best.fun: best = r
    c = best.x / np.linalg.norm(best.x)
    return c * (np.sign(c[0]) if abs(c[0]) > 1e-9 else 1.0)

def geometric_test():
    from mtft.matter.leptons import UP_MZ, DOWN_MZ, LEPTON_POLE
    return {k: m[0] * m[2] / m[1] ** 2 for k, m in (("up", UP_MZ), ("down", DOWN_MZ), ("leptons", LEPTON_POLE))}

JZ = np.diag([1.0, 0.0, -1.0]); JX = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]]) / math.sqrt(2)

def zeeman_ground(beta, lam=0.0, starts=40, seed=1):
    """Minimise -<J_z> cos(beta) - <J_x> sin(beta) - lam |theta|^2 over unit real c = (c+, c0, c-)."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed); best = None
    def E(v):
        c = v / np.linalg.norm(v); return -(c @ JZ @ c) * math.cos(beta) - (c @ JX @ c) * math.sin(beta) - lam * theta(c) ** 2
    for _ in range(starts):
        r = minimize(E, rng.normal(size=3), method="Nelder-Mead", options={"xatol": 1e-13, "fatol": 1e-15, "maxiter": 40000})
        if best is None or r.fun < best.fun: best = r
    c = best.x / np.linalg.norm(best.x)
    return c * np.sign(c[0]), best.fun

def texture_invariants(c):
    """D1/D0 and rho = D2 D0 / D1^2 of the orthonormal autocorrelation."""
    c = np.asarray(c, float); D0 = float(c @ c); D1 = c[1] * c[0] + c[2] * c[1]; D2 = c[2] * c[0]
    return D1 / D0, (D2 * D0 / D1 ** 2 if abs(D1) > 1e-12 else float("inf"))

LAMBDA_CRITICAL = 0.5          # -1 + (2 - 4 lam) b^2 + ... : the coherent state destabilises towards the opposite pole at lam = 1/2

def rho_departures():
    from mtft.matter.leptons import UP_MZ, DOWN_MZ, LEPTON_POLE
    return {k: rho_branches(*m)[0] - 0.5 for k, m in (("up", UP_MZ), ("down", DOWN_MZ), ("leptons", LEPTON_POLE))}

def report():
    inv = fit_invariants()
    return {"invariants": inv, "coherent_quartic": quartic(coherent(0.3)), "polar_quartic": quartic(polar(0.3)),
            "texture_coherent": D_from_alignment(coherent(0.15)), "texture_polar_xw": D_from_alignment(polar(0.0)), "gj": gj_test(),
            "exact_gj_at_rho_half": (sign_flip_factor(0.5), 1 / sign_flip_factor(0.5))}
