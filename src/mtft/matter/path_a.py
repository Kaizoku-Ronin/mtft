"""PA-01 (2026-10-02): Path A — the matter-untwisted completion of YK-01: its vector-like states, what could give them mass, and whether the background
is stable at one loop.

1. Content (EXACT; three-cusp background, M_L = O(P), P a cusp that is not blown up).  Per matter entry the fermions are
       R:     13 A-type  (H^1(X0(143), O(-P)) (x) H^0(143a1, O): constant along the torus)  +  12 Phi-type (H^0(K(-P)) (x) dz)
       Rbar:  12 A-type  (H^1(O(P)) (x) 1)  +  1 A-type (H^0(O(P)) (x) dzbar)  +  2 exceptional jets  +  13 Phi-type (H^0(K(P)) (x) dz)
   25 vector-like pairs and the net 3.
2. Tree-level masses (EXACT).  All 53 states per entry are zero modes.  A gauge-invariant mass for a pair (R, Rbar) needs a vev in a diagonal block — the
   stack-U(1) moduli: Phi-type moduli pair A-type with A-type, Wilson-line moduli pair Phi-type with A-type.  Nothing else in the spectrum can do it: the
   Higgs bidoublet links (4,2,1) to (4bar,1,2), and the Pati–Salam condensate links (4,2,1) to the 48 bidoublet jets only.
3. STABILITY (EXACT at one loop in the product approximation of Matter-05).  The untwisting Wilson line is a point theta_eta on the torus of flat connections
   of the relative U(1) T = (1 on the quark–lepton group, -1 on L and R).  Only the two matter blocks are charged under T, so its one-loop potential is
       V(theta) = -2 c sum_{w : eta(w) = -1} K(w) cos(q theta . w),   K(w) > 0,
   minimal at theta = 0 (matter twisted: jets) and MAXIMAL at theta_eta, where cos(q theta_eta . w) = -1 for every odd winding: the matter-untwisted point is
   a maximum, its Wilson-line modulus is tachyonic (Hessian = minus the matter contribution at the minimum), and the background rolls back to the jet vacuum.
   Nothing charged under T has the opposite sign: the adjoints and the bidoublet are T-neutral.  `potential_profile` evaluates V along the straight path
   0 -> theta_eta with K(w) ~ |w|^-8 (leading order in A_E / R_C^2).
4. Where the vector-like states go (EXACT).  At theta = 0 the product-cohomology modes see the twist eta on the torus: the lowest of them sit at the
   Scherk–Schwarz scale M = 2 pi dist(P_eta, O) / A_E — M sqrt(A_E) = 3.106 (real spin structure) or 3.537 (complex pair) — while the 3 net jets per entry stay
   massless.  So Path A's 25 pairs per entry are the heavy MESSENGERS of the stable vacuum, and the light families (jets) get Yukawas only by exchanging them
   between the Higgs and the Pati–Salam condensate: y_eff ~ y_1 y_2 v_R / M — Path B.  The two paths converge; the effective-Yukawa computation needs the
   couplings jet–messenger–Higgs and messenger–jet–condensate.
OPEN: the effective Yukawa matrices; stabilising the untwisted point instead (a rigid realisation of the twist — e.g. almost-commuting holonomies or an
orbifold projection of the T Wilson line) is a model-building alternative."""
import math
import numpy as np


def content():
    R = {"A-type H1(O(-P)) x 1": 13, "Phi-type H0(K(-P)) x dz": 12}
    Rbar = {"A-type H1(O(P)) x 1": 12, "A-type H0(O(P)) x dzbar": 1, "jets": 2, "Phi-type H0(K(P)) x dz": 13}
    from mtft.matter.yukawa_jets import untwisted_matter_entry
    e = untwisted_matter_entry()
    assert sum(R.values()) == e["R"] and sum(Rbar.values()) == e["Rbar"]
    return {"R": R, "Rbar": Rbar, "net": sum(Rbar.values()) - sum(R.values()), "pairs": sum(R.values())}

def mass_sources():
    return {"stack-U(1) Phi moduli": "pair A-type R with A-type Rbar", "stack-U(1) Wilson-line moduli": "pair Phi-type with A-type",
            "Higgs bidoublet": "links (4,2,1) to (4bar,1,2) (Yukawa, not vector-like mass)", "Pati–Salam condensate": "links (4,2,1) to bidoublet jets (48 states)"}

def T_charges():
    """Charges under T = diag(1 on the quark–lepton group, -1 on L and R) (traceless with ranks 4, 2, 2)."""
    q = {"4": 1, "L": -1, "R": -1}
    return {"adj": 0, "(4,2,1)": q["4"] - q["L"], "(4bar,1,2)": q["4"] - q["R"], "(1,2,2)": q["L"] - q["R"]}

def _lattice(eps, N=60):
    from mtft.research.theta_torus import tau_143a1
    tau = complex(tau_143a1()); s = math.sqrt(tau.imag); out = []
    for m in range(-N, N + 1):
        for n in range(-N, N + 1):
            if (m, n) == (0, 0): continue
            w = (m + n * tau) / s; odd = (eps[0] * m + eps[1] * n) % 2 == 1
            out.append((m, n, abs(w) ** 2, odd))
    return out

def potential_profile(eps=(0, 1), steps=11, N=60):
    """V(s) along theta = s theta_eta (s from 0 to 1), normalised by V(0): only odd windings, weights |w|^-8, phase pi s (eps . (m, n))."""
    L = _lattice(eps, N); V = []
    for k in range(steps):
        s = k / (steps - 1)
        V.append(-2 * sum(r2 ** -4 * math.cos(math.pi * s * (eps[0] * m + eps[1] * n)) for m, n, r2, odd in L if odd))
    V0 = V[0]
    return [(k / (steps - 1), v / abs(V0)) for k, v in enumerate(V)]

def hessian_sign(eps=(0, 1), N=60):
    """Second derivative of V along the path at s = 0 and at s = 1 (theta_eta)."""
    L = _lattice(eps, N)
    d2 = lambda s: 2 * sum(r2 ** -4 * (math.pi * (eps[0] * m + eps[1] * n)) ** 2 * math.cos(math.pi * s * (eps[0] * m + eps[1] * n)) for m, n, r2, odd in L if odd)
    return {"at_jet_vacuum": d2(0.0), "at_untwisted_point": d2(1.0)}

def messenger_masses():
    """Lowest eta-twisted mode on the unit-area 143a1: M sqrt(A_E) = 2 pi dist(P_eta, O)."""
    from mtft.research.theta_torus import tau_143a1
    tau = complex(tau_143a1()); s = math.sqrt(tau.imag)
    pts = {"real (P = 1/2)": 0.5, "complex pair (P = tau/2)": tau / 2}
    out = {}
    for k, P in pts.items():
        d = min(abs(P - (m + n * tau)) for m in (-1, 0, 1) for n in (-1, 0, 1)) / s
        out[k] = 2 * math.pi * d
    return out

def report():
    return {"content": content(), "T_charges": T_charges(), "profile": potential_profile(), "hessian": hessian_sign(), "messengers": messenger_masses()}
