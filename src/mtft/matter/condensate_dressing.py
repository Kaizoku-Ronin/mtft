"""CD-01 (2026-10-05): the condensate dresses u^c and e^c but never d^c — and that carries the up–down asymmetry.

1. Gauge selection (EXACT).  The condensate points at nu^c: SU(4) lepton index, SU(2)_R upper component.  Its gauge–fermion mixing, sqrt(2) g Sigma* lambda psi, reaches
   a family fermion psi only if one broken generator carries the condensate's direction to psi's.  One step reaches u^c (a leptoquark generator: same SU(2)_R
   component, colour index), e^c (an SU(2)_R step: same SU(4) index) and nu^c itself (the broken diagonal U(1)).  d^c needs two steps and Q is another block: the
   dressing acts on u^c, e^c, nu^c only.
2. What the dressing does (STRUCTURAL).  The admixed gaugino modes have no Yukawa coupling to Q H, so each dressed family's entries are multiplied by its zero-mode
   weight Z_j.  theta_L^u = Y_23/Y_33 is untouched (Z_3 cancels), but m_c/m_t and theta_R^u pick up z = Z_2/Z_3: FL-2's shared theta_R fails for the up sector.
3. What the data now need (EXACT given the induced texture).  m_s/m_b = t_d theta_R, V_cb = t_d - t_u, m_c/m_t = t_u theta_R z give z = (m_c/m_t)/((m_s/m_b) - V_cb theta_R).
   z = 1 reproduces FL-2's 40 : 5.3 : 1.  At z = 0.25 the geometric tails are nearly symmetric — Q : d^c : u^c = 0.71 : 1 : 0.75 — with Q the shortest, as VS-3's
   stability ordering requires.  The dressing can carry the up–down asymmetry the geometry could not.
4. Direction and a bonus (STRUCTURAL).  z < 1 needs the charm's u^c dressed more than the top's; the condensate's jets sit at the blown-up cusp, beside the second
   family — the right direction.  The same dressing acts on e^c but not d^c: the lepton-specific dressing CP-4's Georgi–Jarlskog hypothesis needed.
OPEN: Z_2/Z_3 from the condensate and gaugino profiles at the two points; the e^c dressing and the Georgi–Jarlskog factors it gives."""
import itertools

COMPONENTS = [(s4, s2) for s4 in ("c1", "c2", "c3", "4") for s2 in ("up", "down")]
CONDENSATE = ("4", "up")
NAMES = {("4", "up"): "nu^c", ("4", "down"): "e^c"}


def name(comp):
    return NAMES.get(comp, "u^c" if comp[1] == "up" else "d^c")

def one_step(a, b):
    """Can a single generator of SU(4) x SU(2)_R x U(1) map component a to b? (changes at most one label)"""
    return (a[0] == b[0]) or (a[1] == b[1])

def dressed():
    return sorted({name(c) for c in COMPONENTS if one_step(CONDENSATE, c)})

def undressed():
    return sorted({name(c) for c in COMPONENTS if not one_step(CONDENSATE, c)})

def z_required(theta_R, mc_mt=None, ms_mb=None, V=0.041):
    if mc_mt is None:
        from mtft.matter.third_observable import sector_ratios
        mc_mt, ms_mb = sector_ratios()
    t_d = ms_mb / theta_R; t_u = t_d - V
    return mc_mt / (t_u * theta_R), (theta_R, t_d, t_u)

def report():
    z1, tails1 = z_required(0.381)
    zq, tailsq = z_required(0.117)
    return {"dressed": dressed(), "undressed": undressed(), "z_at_FL2": z1, "tails_FL2": tails1, "z_symmetric": zq, "tails_symmetric": tailsq}
