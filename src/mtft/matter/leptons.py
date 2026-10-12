"""LEP-01 (2026-10-02): the charged-lepton texture — the first place the construction touches an atom.

1. The Georgi–Jarlskog problem at the matching scale (NUMERICAL, SM one loop).  At the b–tau scale (1.1e6 GeV, PS-01) m_tau = m_b by construction, but
   m_mu / m_s = 3.06 and m_e / m_d = 0.295: the light families need factors ~3 and ~1/3 that an SU(4) relation cannot give.  No common rescaling of the
   hierarchy parameter helps: the second family needs eps_e / eps_d = 1.75, the first 0.74.
2. When Y_e = Y_d is exact (EXACT).  The charged-lepton and down-quark Yukawas use the same Higgs entry (L, B), the same heavy bidoublet modes and — when the flux
   keeps SU(4) (equal gamma-charges, curve degrees and exceptional charges on colour and A1) — identical jet representatives.  Every ingredient of the two vertices
   then coincides: Y_e = Y_d entry by entry, in every SU(4)-preserving background of FS-01.
3. The GJ-friendly order (STRUCTURAL).  An SU(4)-breaking gamma-flux (colour != A1) reshapes the lepton jets' harmonic representatives at the same order as the
   effects that generate the light families (eps), while the third family moves only at relative order eps^2.  So b–tau survives and the light-family ratios can
   move by order one — the pattern Georgi–Jarlskog needs.  Whether the factors come out near 3 and 1/3 needs the actual overlaps.
4. The equal-overlap texture (MODEL; corrected by CC-61).  In the idealised texture with equal overlap coefficients — Yukawa entries equal to the
   autocorrelations D_n = sum_q s_{q+n} s_q of the condensate's jet amplitudes s — each sector fits exactly, but not uniquely: the texture is even in s1 and
   blind to zero flips (CC-40), so every sector has several exact solutions (all_alignments).  On the rho branch nearest the coherent value 1/2 every sector
   has a charge-0-dominated fit — up (1, 0.082, 0.003), down (1, 0.221, 0.029), leptons (1, 0.339, 0.031) — and the leptons also fit with the charge-1 jet
   dominant, (1, 3.26, -0.31).  The texture does not decide which jet dominates; the earlier reading "leptons need the charge-1 jet" (LEP-4) was the optimizer's
   path and is retracted.  The coherent family s = (1, t, t^2) fails every sector (its two light masses stay within a factor two of each other).
5. Koide (NUMERICAL).  Q = 0.666664 for the charged leptons (pole masses), 0.746 for down and 0.888 for up quarks (running at m_Z).  With Y_e = Y_d an SU(4)
   texture would give leptons the down-quark Q, so Koide for leptons is one more face of item 1.  Nothing in the texture singles out Q = 2/3; a fit to two mass
   ratios satisfies it automatically.  Caution recorded: imposing s_1 s_2 = -1, noticed AFTER the lepton fit (-1.017), predicts m_e 5.5 per cent low and Q = 0.6672 —
   a family chosen after looking at the data, so it is numerology, not evidence.
6. Correction to RT-01 (CC-39).  Gauge invariance makes the gamma-charges of the three entries in a Yukawa vertex sum to zero, so the up and down vertices differ in
   two of their three charges (ratio-3 background: up (1, 0, 1), down (1, 2, 3)), and the Higgs entries also differ in curve degree (+1, -1).  The 'exactly 3 at
   leading order' of RT-1 is therefore an estimate: the ratio is near 3, not fixed at 3.  The plane-wave picture of FS-4 is a gauge artefact; the gauge-invariant
   effect is the reshaping of wavefunctions by the gamma field strength.
OPEN: the lepton dressing by the condensate (it sits in the lepton direction of the 4); the effective potential for the condensate's jet amplitudes s (the
Legendre-transform object that would fix s, eps and the textures)."""
import functools
import math
import numpy as np

LEPTON_POLE = (0.51099895, 105.6583755, 1776.93)        # MeV
DOWN_MZ = (0.0027, 0.055, 2.86)                            # GeV, MS-bar at m_Z (inputs of PS-01)
UP_MZ = (0.00127, 0.619, 171.0)

def koide(m):
    m = np.asarray(m, float); return float(m.sum() / np.sqrt(m).sum() ** 2)

def gj_at_matching():
    from mtft.matter.ps_breaking import b_tau_running
    bt = b_tau_running()
    r2 = lambda m: m[1] / m[2]; r1 = lambda m: m[0] / m[2]
    return {"mu_star_GeV": bt["mu_equal"], "mmu_over_ms": 1 / bt["ms_over_mmu_at_mu"], "me_over_md": 1 / bt["md_over_me_at_mu"],
            "eps_ratio_second_family": math.sqrt(r2(LEPTON_POLE) / r2(DOWN_MZ)), "eps_ratio_first_family": (r1(LEPTON_POLE) / r1(DOWN_MZ)) ** 0.25}

def vertex_charges(flux):
    """gamma-charges of (doublet, singlet, Higgs) entries in the down-quark and charged-lepton vertices of a FS-01 StackFlux (sum to zero by gauge invariance)."""
    G = flux.G
    down = (G("c") - G("L"), G("B") - G("c"), G("L") - G("B")); lep = (G("A1") - G("L"), G("B") - G("A1"), G("L") - G("B"))
    up = (G("c") - G("L"), G("A2") - G("c"), G("L") - G("A2"))
    return {"up": up, "down": down, "lepton": lep}

def su4_identity_holds(rows):
    """In every SU(4)-preserving scan background the lepton and down vertices have identical charge triples and curve-degree data."""
    from mtft.matter.flux_scan import StackFlux
    ok = True
    for r in rows:
        if not r["su4"]: continue
        f = StackFlux(r["a"], r["g"]); v = vertex_charges(f)
        ok &= v["down"] == v["lepton"] and f.A("c") == f.A("A1")
    return ok


# ------------------------------------------------------------------ the idealised autocorrelation texture
def autocorrelation_texture(s):
    s = np.asarray(s, float); D = [sum(s[q + n] * s[q] for q in range(3 - n)) for n in range(3)]
    return np.array([[D[a + b] if a + b <= 2 else 0.0 for b in range(3)] for a in range(3)])

def mass_ratios(s):
    m = np.sort(np.linalg.svd(autocorrelation_texture(s), compute_uv=False))[::-1]
    return m[1] / m[0], m[2] / m[0], m

def _rho(s):
    D = [sum(s[q + n] * s[q] for q in range(3 - n)) for n in range(3)]
    return D[2] * D[0] / D[1] ** 2

@functools.lru_cache(maxsize=None)
def all_alignments(masses, tol=1e-10):
    """Every exact equal-overlap solution s = (1, s1, s2) reached from a fixed grid of starts (s1 in [-4, 4], s2 in [-1, 1]).

    The texture is even in s1, so solutions are returned with s1 >= 0, deduplicated and sorted by s1.  Zero flips (CC-40) and the two rho branches (CP-6)
    make the solution set larger than one point in every sector; LEP-6 records the sets."""
    from scipy.optimize import least_squares
    r2t, r1t = masses[1] / masses[2], masses[0] / masses[2]

    def f(v):
        r2, r1, _ = mass_ratios([1, v[0], v[1]]); return [math.log(r2 / r2t), math.log(r1 / r1t)]

    found = []
    for a in np.linspace(-4.0, 4.0, 33):
        for b in np.linspace(-1.0, 1.0, 17):
            try:
                r = least_squares(f, (a, b), xtol=1e-14, ftol=1e-14)
            except ValueError:                                   # a start where a mass ratio vanishes
                continue
            if max(abs(x) for x in r.fun) < tol:
                s = (1.0, abs(float(r.x[0])), float(r.x[1]))
                if not any(abs(s[1] - t[1]) < 1e-6 * (1 + abs(t[1])) and abs(s[2] - t[2]) < 1e-6 * (1 + abs(t[2])) for t in found):
                    found.append(s)
    return tuple(sorted(found, key=lambda s: s[1]))

def fit_alignment(masses):
    """A canonical exact solution: the rho branch nearest the coherent value 1/2 (CP-11), then the smallest coefficients (the charge-0-dominated one).

    Returns (s, residual).  The representative is a convention, not a result: all_alignments gives the full set (CC-61)."""
    masses = tuple(float(x) for x in masses)
    sols = all_alignments(masses)
    if not sols:
        raise ValueError("no exact equal-overlap solution found")
    best = min(sols, key=lambda s: (round(abs(_rho(s) - 0.5), 9), max(abs(s[1]), abs(s[2]))))
    r2, r1, _ = mass_ratios(best)
    residual = max(abs(math.log(r2 / (masses[1] / masses[2]))), abs(math.log(r1 / (masses[0] / masses[2]))))
    return best, float(residual)

def coherent_light_mass_spread(ts=(0.05, 0.1, 0.2, 0.3, 0.5)):
    return {t: mass_ratios([1, t, t * t])[0] / mass_ratios([1, t, t * t])[1] for t in ts}

def numerology_example():
    """The family s = (1, t, -1/t), chosen after seeing the lepton fit: fit m_mu/m_tau, read off m_e/m_tau and Q."""
    from scipy.optimize import brentq
    r2t = LEPTON_POLE[1] / LEPTON_POLE[2]
    t = brentq(lambda t: mass_ratios([1, t, -1 / t])[0] - r2t, 2.5, 4.5)
    r2, r1, m = mass_ratios([1, t, -1 / t])
    return {"t": t, "me_over_mtau_pred": r1, "me_over_mtau_data": LEPTON_POLE[0] / LEPTON_POLE[2], "koide": koide(m)}

def report():
    return {"gj": gj_at_matching(), "koide": {"leptons": koide(LEPTON_POLE), "down": koide(DOWN_MZ), "up": koide(UP_MZ)},
            "fits": {k: fit_alignment(m)[0] for k, m in (("up", UP_MZ), ("down", DOWN_MZ), ("leptons", LEPTON_POLE))},
            "coherent_spread": coherent_light_mass_spread(), "numerology": numerology_example()}
