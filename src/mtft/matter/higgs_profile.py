"""HP-01 (2026-10-04): the light Higgs of the ratio-3 background — where it lives, what type it is, and why it cannot carry the Yukawas.

1. Where (EXACT).  Ratio-3's Higgs entries carry gamma-charges Delta_u = 1 ((L, A2)) and Delta_d = 3 ((L, B)).  Over a curve point x an entry restricts to the
   torus as the degree-0 bundle P_{Delta phi(x)}, which has cohomology only when Delta phi(x) = O.  So R^0 pi_* vanishes and R^1 pi_* is a skyscraper on
   phi^-1(E[Delta]): 4 Delta^2 modes (deg phi = 4).  H_u: the four cusps (phi^-1(O); E(Q)_tors = 0), one of them the blown-up cusp.  H_d: 36 points over the
   nine 3-torsion points, four of them the cusps.  Near the blown-up cusp the H_u mode is organised by the graph w = phi(q), with
   phi(q) = sum a_n q^n / n = q - q^3/3 - q^4/2 - q^5/5 - 2q^7/7 + ... for 143a1 (a_2 = 0).
2. Type (EXACT).  The cohomology is H^1 only (chi = -4 Delta^2 - e(e-1)/2, H^0 = H^2 = 0): the modes are A-type, (0,1)-forms.  There is no Phi-type,
   global torus-constant Higgs.
3. No tree-level Yukawa (EXACT).  The families are A-type jets, and the 8D superpotential Tr Phi ^ [A, A] needs exactly one Phi-type field, so Q u^c H with all
   three A-type vanishes whatever the profiles; PB-01's D-term mixing also needs a Phi-type Higgs.  FS-01 item 3's linear plane-wave entry and RT-01's ratio-3
   mechanism therefore have no coupling (CC-43).
4. Where it works (EXACT).  The net families are the degree-2 coboundary jets of the e = -3 orientation, so the coupling Higgs carries exceptional charge +6.
   With a gamma-neutral bidoublet its Phi-type modes are sections of K (x) O(6E) (x) pi*L = pi*(K L) (x) O(7E), i.e. pullbacks of global forms: torus-constant
   times holomorphic sections of K_C (x) L_a on X_0(143).  Their Taylor data at the cusp are q-expansion coefficients — the HJ-01 structure.  Backgrounds with
   every Higgs entry global (FS-01's minimal family) qualify; ratio 3 does not.
OPEN: in a gamma-neutral background, which Phi-type Higgs combination is light and takes the vev, and its q-expansion at the cusp."""
from fractions import Fraction

RATIO3_HIGGS = {"H_u": {"block": "(L, A2)", "gamma": 1}, "H_d": {"block": "(L, B)", "gamma": 3}}
MODULAR_DEGREE = 4
CUSPS_OVER_O = 4


def localisation(delta):
    """Mode count and where the modes sit for a gamma-charge-delta entry."""
    n = MODULAR_DEGREE * delta * delta
    return {"modes": n, "torsion_points": delta * delta, "at_cusps": CUSPS_OVER_O, "elsewhere": n - CUSPS_OVER_O}

def phi_taylor(nmax=8):
    """Exact Taylor coefficients a_n / n of the modular parametrisation at the cusp."""
    from mtft.periods.forms import _elliptic_an, CURVE_143A1
    a = _elliptic_an(CURVE_143A1, 143, nmax)
    return [Fraction(a[n], n) for n in range(1, nmax + 1)]

def chi_entry(gamma, e):
    """chi of an entry with gamma-charge gamma and exceptional charge e, from the surface intersection form."""
    from mtft.surface.blowup import chi
    return chi((0, 0, gamma, e))

def superpotential_allowed(types):
    """Tr Phi ^ [A, A]: a cubic coupling needs exactly one Phi-type field."""
    return sum(1 for t in types if t == "Phi") == 1

def coupling_higgs_charge(e_family=-3):
    """Exceptional charge of the Higgs orientation that closes Q u^c H when both families sit in the e_family orientation."""
    return -2 * e_family

def report():
    return {"H_u": localisation(1), "H_d": localisation(3), "phi_taylor": [str(x) for x in phi_taylor()],
            "chi_u": chi_entry(1, 6), "chi_d": chi_entry(3, 6), "ratio3_tree_yukawa": superpotential_allowed(("A", "A", "A")),
            "global_tree_yukawa": superpotential_allowed(("A", "A", "Phi")), "coupling_higgs_charge": coupling_higgs_charge()}
