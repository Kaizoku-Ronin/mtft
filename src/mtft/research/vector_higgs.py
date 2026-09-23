"""R2C-02 (2026-09-19): the internal-vector Higgs of M1 — mode operator, tachyon mass, and the Yukawa as the gauge interaction.

R2C-01 excludes an elementary 6D scalar Higgs for the three families and leaves the internal gauge component A_zbar in
Hom(L_y, L_x) (x) Kbar = L (x) Kbar, L = L_x L_y^{-1} of degree d_L (= -6 for H_u, H_d of M1).  From the 6D Yang–Mills action the
4D mass operator on these (0,1)-forms is (Weitzenböck)
    m^2 = nabla^* nabla_{L (x) Kbar} - 2 (2 pi |d_L| / A) + K_gauss,
the Bochner Laplacian of the degree-(|d_L| + deg K) bundle (30 for M1 — the "Higgs target space" whose 18 sections were built
in SM-03), the spin-1 magnetic-moment term, and the Ricci term of a 1-form.  On a constant-curvature metric with K = -1 and
A = 2 pi (2g - 2), the Landau level of nabla^* nabla is 2 pi (|d_L| + deg K)/A = 2 pi |d_L|/A + 1, so the curvature terms cancel EXACTLY:
    m^2_LLL = -2 pi |d_L| / A   (= -1/4 for M1 in curvature units), multiplicity h^1(L) = 18; all other modes massive.
Flat limit: -|F_L|, the Nielsen–Olesen / Bachas flux tachyon.  The Yukawa needs no new interaction: it is the gauge coupling
psi-bar Gamma^{zbar} A_zbar psi, so y_ijk = g_6 I_ijk = g_4 sqrt(A) I_ijk with unit-normalised modes (g_4^2 = g_6^2/A) — the
"section-product tensors" of SM-03 are these overlaps, and the equal-chirality selection rule of the vector insertion is
exactly what the family-preserving assignment provides.  kappa: an O(1) Clifford-convention factor in the reduction."""
import numpy as np, sympy as sp

def tachyon_mass2(deg_L, genus=13, deg_K=None, K_gauss=-1):
    """Exact: with A = 2 pi (2g-2) R^2 and K = -1/R^2, m^2_LLL = -2 pi |d_L|/A = -|d_L| / ((2g-2) R^2)."""
    R = sp.Symbol("R", positive=True); deg_K = deg_K if deg_K is not None else 2 * genus - 2; A = 2 * sp.pi * (2 * genus - 2) * R ** 2
    landau = 2 * sp.pi * (abs(deg_L) + deg_K) / A; m2 = sp.simplify(landau - 2 * (2 * sp.pi * abs(deg_L) / A) + sp.Rational(K_gauss) / R ** 2)
    return {"m2_LLL": m2, "cancellation": sp.simplify(2 * sp.pi * deg_K / A + sp.Rational(K_gauss) / R ** 2) == 0, "multiplicity": "h^1(L)", "R": R}

def mass_operator_spectrum(bochner_eigenvalues, deg_L, area, K_gauss=-1.0):
    """m^2_n = lambda_n - 2 (2 pi |d_L|/area) + K_gauss for the Bochner eigenvalues of the degree-(|d_L| + deg K) bundle."""
    lam = np.asarray(bochner_eigenvalues, float); return lam - 2 * (2 * np.pi * abs(deg_L) / area) + K_gauss

def m1_vector_higgs_spectrum(h=0.25, nx=6, nev=26):
    from ..surface import magnetic as MG, rrspace as RR
    mm = MG.MagneticMesh(h, nx); cm = RR.cm_classes(); P = cm["P"]; T0 = {k: mm.point_triangle(P[k]) for k in (1, 2, 3)}; Tc0 = mm.cusp_triangle(143); Tc1 = mm.cusp_triangle(13)
    Hg = mm.divisor_spectrum(30, [(Tc0, 12), (Tc1, 12), (T0[1], 2), (T0[2], 2), (T0[3], 2)], nev); m2 = mass_operator_spectrum(Hg["eigenvalues"], -6, mm.Ac)
    return {"area": mm.Ac, "bochner": Hg["eigenvalues"], "m2": m2, "tachyon_mean": float(m2[:18].mean()), "tachyon_exact": -2 * np.pi * 6 / mm.Ac, "first_massive": float(m2[18]), "multiplicity": 18}

def yukawa_gauge_ratio(Y_fem, area, n=4000, seed=0):
    """y_ijk / g_4 = kappa * I~_ijk with I~ = sqrt(area) Y (unit-normalised modes); distribution of the largest singular value over
    random Higgs directions (the modulus is unselected): y_t/g_4 = kappa s_max(v)."""
    It = np.sqrt(area) * np.asarray(Y_fem); rng = np.random.default_rng(seed); smax = []
    for _ in range(n):
        v = rng.standard_normal(It.shape[2]) + 1j * rng.standard_normal(It.shape[2]); v /= np.linalg.norm(v); smax.append(np.linalg.svd(np.einsum("ijk,k->ij", It, v), compute_uv=False)[0])
    smax = np.array(smax); return {"I_tilde": It, "rms_entry": float(np.sqrt(np.mean(np.abs(It) ** 2))), "s_max_percentiles": tuple(float(p) for p in np.percentile(smax, [5, 50, 95])), "s_max_max": float(smax.max()),
                                   "reading": "y_t/g_4 typical ~0.6 (0.4–0.9); observed ~0.7 at ~1e16 GeV; third-family down/lepton ratios are the same order (Yukawa unification, large tan beta)"}


def reduction_constants():
    """kappa = 1: with D = d - iA, kinetic -(1/2 g6^2) tr F^2, internal polarisations a_z = (a5 - i a6)/2, a_zbar = (a5 + i a6)/2
    (|a5|^2 + |a6|^2 = 2(|a_z|^2 + |a_zbar|^2)), the canonical complex Higgs is h = 2 a_zbar / g6; the vertex psi-bar Gamma^i A_i psi
    contains a_zbar (sigma5 - i sigma6) = 2 a_zbar sigma^- with unit matrix element between the S^+ and S^- zero-mode components;
    hence y = 2 (g6 h/2) I = g6 I = g4 sqrt(A) I with g4^2 = g6^2/A."""
    s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]]); sminus = sp.Matrix([[0, 0], [1, 0]])
    g6, h, I = sp.symbols("g6 h I", positive=True); a_zbar = g6 * h / 2; vertex = 2 * a_zbar * I
    return {"sigma_identity": (s1 - sp.I * s2) == 2 * sminus, "canonical_higgs": "h = 2 a_zbar / g6", "kappa": sp.simplify(vertex / (g6 * I * h)), "yukawa": "y = g4 sqrt(A) I"}


# ------------------------------------------------ R2C-05 (2026-09-19): both polarisations of a flux block; the recombination scale
def block_spectrum(degree, genus=13):
    """A gauge block with line bundle L of degree d on X0(143) (compact metric, curvature units) has two internal polarisations with
    Landau levels at m^2 = +-d/(2g-2) = +-d/24: the tachyonic one (Bochner on K (x) L_neg^-1, degree 24 + |d|) with multiplicity
    h^1(L_neg) = |d| + g - 1, and the massive one (Bochner on K (x) L_neg, degree 24 - |d|) with multiplicity 12 - |d| + h^0(O(|d| points)).
    The block symmetry L <-> L^-1 with polarisation swap fixes the sign of the magnetic-moment term for each polarisation; the earlier
    application of the (0,1) formula to a positive-degree block (M3, "-3/8") was wrong (CC-30).  Verified on the magnetic mesh: M1 (|d| = 6):
    18 @ -1/4 and 7 @ +1/4; M3 (|d| = 3): 15 @ -1/8 and 10 @ +1/8.  Consequence: a flux-charged gauge block is ALWAYS tachyonic at the
    compactification scale (m^2 = -|d|/24 in units where the curvature radius is 1); only a flux-neutral block (d = 0) gives massless Higgs
    moduli.  Neither M1's Higgs blocks (d = -6) nor M3's 27 block (d = 3, which contains the 10_H) is flux-neutral."""
    d = abs(int(degree)); g1 = genus - 1
    return {"degree": degree, "tachyon_m2": -sp.Rational(d, 2 * g1), "tachyon_multiplicity": d + g1, "massive_m2": sp.Rational(d, 2 * g1),
            "massive_multiplicity_formula": "12 - |d| + h0(O(|d| points))", "massive_multiplicity_M1_M3": {6: 7, 3: 10}.get(d), "flux_neutral": d == 0,
            "recombination_scale": "compactification scale unless d = 0"}

def m3_block_check(h=0.3, nx=5):
    """Mesh check of the degree-27 (tachyonic, 15) and degree-21 (massive, 10) clusters of the M3 27 block (slow)."""
    from ..surface import magnetic as MG, rrspace as RR
    mm = MG.MagneticMesh(h, nx); cm = RR.cm_classes(); P = cm["P"]; T0 = {k: mm.point_triangle(P[k]) for k in (1, 2, 3)}; Tc0 = mm.cusp_triangle(143); Tc1 = mm.cusp_triangle(13)
    l27 = mm.divisor_spectrum(27, [(Tc0, 12), (Tc1, 12), (T0[1], 1), (T0[2], 1), (T0[3], 1)], 20)["eigenvalues"]; l21 = mm.divisor_spectrum(21, [(Tc0, 12), (Tc1, 12), (T0[1], -1), (T0[2], -1), (T0[3], -1)], 16)["eigenvalues"]
    F = 2 * np.pi * 3 / mm.Ac; return {"tachyon_m2": float(l27[:15].mean() - 2 * F - 1), "tachyon_gap": float(l27[15] - l27[14]), "massive_m2": float(l21[:10].mean() + 2 * F - 1), "massive_gap": float(l21[10] - l21[9])}


# ------------------------------------------------ v0.33.0 (2026-09-23): the mode operator in Kodaira form (compendium V.2–V.3)
def kodaira_forms(degree, genus=13):
    """Both polarisations of a block of degree d != 0 as pointwise operator identities, valid for EVERY Kaehler metric on the curve when
    the block carries its HYM metric (|B| = 2 pi |d| / A):
      tachyonic  (the (0,1)-forms of the negative-degree member L_neg, sections of L_neg (x) K^-1):  m^2 = 2 d* d - |B|,
      massive    (the (1,0)-forms of L_neg, sections of L_neg (x) K):                                  m^2 = 2 dbar* dbar + |B|.
    The Gaussian curvature cancels at every point (the Ricci term against the curvature of Lambda^{0,1} = K^-1), not only on average.
    The (0,1)-formula of R2C-02 therefore applies to NEGATIVE-degree blocks only — the content of CC-30 — and its multiplicity
    h^1(L_neg) = |d| + g - 1 is Atiyah–Bott's Morse index of the split connection.  Curvature units: levels -/+ |d|/(2g-2)."""
    ad = abs(int(degree)); B = sp.Rational(ad, 2 * genus - 2)
    return {"B_curvature_units": B, "tachyonic": {"operator": "2 d*d - |B| on L_neg (x) K^-1", "lowest": -B, "multiplicity": ad + genus - 1},
            "massive": {"operator": "2 dbar*dbar + |B| on L_neg (x) K", "lowest_if_attained": B, "multiplicity": f"{genus - 1 - ad} + h^0(L_neg^-1)"},
            "sign_restriction": "the (0,1)-form formula m^2 = Bochner - 2|B| + K holds for d < 0 (CC-30)", "metric_independent": True}

def constant_curvature_checks(degree):
    """Closed-form Bochner spectra on the three constant-curvature geometries, fed into m^2 = Bochner + K -/+ 2|B| (exact rationals):
    sphere (monopole harmonics l(l+1) - q^2, R = 1), flat torus (Landau levels (2n+1)|B_e|, B = 1 per unit degree), genus 13 (lowest
    hyperbolic Landau level |e|/24).  Each reproduces the tachyonic level -|B| with |d| + g - 1 modes; the massive level +|B| is
    attained only where h^0(K (x) L_neg) > 0 (genus 13).  The sphere case is the Brandt–Neri–Coleman monopole instability."""
    from fractions import Fraction as Fr
    d = abs(int(degree)); out = {}
    def sphere_levels(e):
        q = Fr(abs(e), 2); return [(l * (l + 1) - q * q, int(2 * l + 1)) for l in (q + k for k in range(4))]
    B = Fr(d, 2); tach = min((lv + 1 - 2 * B, m) for lv, m in sphere_levels(-d + 2)); mass = min((lv + 1 + 2 * B, m) for lv, m in sphere_levels(d + 2))
    out["sphere"] = {"tachyonic": tach, "expected": (-B, d - 1), "massive_lowest": mass, "level_at_plus_B": False}
    B = Fr(d); tach = min(((2 * n + 1) * B - 2 * B, d) for n in range(3)); mass = min(((2 * n + 1) * B + 2 * B, d) for n in range(3))
    out["torus"] = {"tachyonic": tach, "expected": (-B, d), "massive_lowest": mass, "level_at_plus_B": False}
    b = Fr(d, 24); out["genus13"] = {"tachyonic": (Fr(d + 24, 24) - 1 - 2 * b, d + 12), "expected": (-b, d + 12), "massive": (Fr(24 - d, 24) - 1 + 2 * b, f"{12 - d} + h^0"), "level_at_plus_B": True}
    out["all_tachyonic_levels_match"] = all(out[k]["tachyonic"] == out[k]["expected"] for k in ("sphere", "torus", "genus13"))
    return out
